#include <iostream>
#include <fstream>
#include <math.h>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <map>
#include <numeric>
#include <functional>
#include <string.h>
#include <sstream>
#include <pthread.h>
#include <thread>
#include "zlib.h"
#include <sys/stat.h>

using namespace std;
//compile with openMP g++ linearBCV.cpp -O3 -std=c++11 -mavx2 -msse4.2 -pthread -fopenmp -lz -o linear11

//some global variables
int RAND_SAMPLES; //will all be set later (if needed)
int image_m; int image_n; int image_o; int image_d=1;
float SSD0=0.0; float SSD1=0.0; float SSD2=0.0; float distfx_global; float beta=1;
//float SIGMA=8.0;
int qc=1;

//struct for multi-threading of mind-calculation
struct mind_data{
	float* im1;
    float* d1;
    uint64_t* mindq;
    int qs;
    int ind_d1;
};

float quantile(float* values,int length,float quant1){
    float* values2=new float[length];
    for(int i=0;i<length;i++){
        values2[i]=values[i];
    }
    int quantind=length*min(max(quant1,0.01f),0.99f);
    //printf("quantind: %d/%d\n",quantind,length);
    nth_element(values2,values2+quantind,values2+length);
    float med1=values2[quantind];
    // delete values2;
    return med1;
}
bool RIGID=false;


struct parameters{
    float alpha; int levels; bool segment,affine,rigid;
    vector<int> grid_spacing; vector<int> search_radius;
    vector<int> quantisation;
    string fixed_file,moving_file,output_stem,moving_seg_file,affine_file,deformed_file;
};


#include "imageIOgzType.h"
#include "transformations.h"
#include "QRsolve.h"
#include "affineLTS1.h"
#include "MINDSSCbox.h"
#include "dataCostD.h"
#include "parseArguments.h"

int main (int argc, char * const argv[]) {
    
    //PARSE INPUT ARGUMENTS
    
    if(argc<4||argv[1][1]=='h'){
        cout<<"=============================================================\n";
        cout<<"Usage (required input arguments):\n";
        cout<<"./linearBCV -F fixed.nii.gz -M moving.nii.gz -O output\n";
        cout<<"optional parameters:\n";
        cout<<" -R <find rigid instead of affine transform> (default 0)\n";
        cout<<" -l <number of levels> (default 4)\n";
        cout<<" -G <grid spacing for each level> (default 7x6x5x4)\n";
        cout<<" -L <maximum search radius - each level> (default 5x4x3x2)\n";
        cout<<" -Q <quantisation of search step size> (default 4x3x2x1)\n";
        cout<<"=============================================================\n";
        return 1;
    }
    parameters args;
    //defaults
    args.grid_spacing={7,6,5,4};
    args.search_radius={5,4,3,2};
    args.quantisation={4,3,2,1};
    args.levels=4;
    parseCommandLine(args, argc, argv);
    
    size_t split_fixed=args.fixed_file.find_last_of("/\\");
    if(split_fixed==string::npos){
        split_fixed=-1;
    }
    size_t split_moving=args.moving_file.find_last_of("/\\");
    if(split_moving==string::npos){
        split_moving=-1;
    }
    
    
    if(args.fixed_file.substr(args.fixed_file.length()-2)!="gz"){
        cout<<"images must have nii.gz format\n";
        return -1;
    }
    if(args.moving_file.substr(args.moving_file.length()-2)!="gz"){
        cout<<"images must have nii.gz format\n";
        return -1;
    }
    
    cout<<"Starting linear reg. of "<<args.fixed_file.substr(split_fixed+1)<<" and "<<args.moving_file.substr(split_moving+1)<<"\n";
    cout<<"=============================================================\n";

    if(args.rigid){
        cout<<"Searching only for rigid (6 parameter) transform\n";
        RIGID=true;
    }

	
	qc=2;
    
	float alpha=1;
    //READ IMAGES and INITIALISE ARRAYS
    
    timeval time1,time2,time1a,time2a;
    
    RAND_SAMPLES=1; //fixed/efficient random sampling strategy
    
    float* im1; float* im1b;
    
    int M,N,O,P; //image dimensions
    
    //==ALWAYS ALLOCATE MEMORY FOR HEADER ===/
    char* header=new char[352];
    
    readNifti(args.fixed_file,im1b,header,M,N,O,P);
    image_m=M; image_n=N; image_o=O;

    readNifti(args.moving_file,im1,header,M,N,O,P);
    
    if(M!=image_m|N!=image_n|O!=image_o){
        cout<<"Inconsistent image sizes (must have same dimensions)\n";
        return -1;
    }
    
    int m=image_m; int n=image_n; int o=image_o; int sz=m*n*o;
    
    //assume we are working with CT scans (add 1024 HU)
    float thresholdF=-1024; float thresholdM=-1024;
    
    for(int i=0;i<sz;i++){
        im1b[i]-=thresholdF;
        im1[i]-=thresholdM;
    }
    
	float *warped1=new float[m*n*o];
    float *warped2=new float[m*n*o];
    
	int step1; int hw1; float quant1;
	
    vector<int> mind_step;
    for(int i=0;i<args.quantisation.size();i++){
        mind_step.push_back(floor(0.5f*(float)args.quantisation[i]+1.0f));
    }
    float* X=new float[4*4];
    
    float Xprev[16]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    for(int i=0;i<16;i++){
        X[i]=Xprev[i];
    }
    
    uint64_t* im1_mind=new uint64_t[m*n*o];
    uint64_t* im1b_mind=new uint64_t[m*n*o];
    uint64_t* warped_mind=new uint64_t[m*n*o];
	
	gettimeofday(&time1a, NULL);
    float timeDataSmooth=0;
	//==========================================================================================
	//==========================================================================================
	float* bench=new float[6*args.levels];
	for(int level=0;level<args.levels;level++){
        quant1=args.quantisation[level];
        step1=args.grid_spacing[level];
        hw1=args.search_radius[level];

        float Xinv[16]; float Ident[16]={1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
        qrsolve(Xinv,Xprev,Ident,4,4);
        
        warpAffine(warped2,im1b,Xinv,m,n,o);
        warpAffine(warped1,im1,Xprev,m,n,o);

        
        float prev=mind_step[max(level-1,0)];
        float curr=mind_step[level];
        
        float timeMIND=0; float timeSmooth=0; float timeData=0; float timeTrans=0;
		
        if(level==0|prev!=curr){
            gettimeofday(&time1, NULL);
            descriptor(im1_mind,im1,m,n,o,mind_step[level]);//max(min(quant1,2.0f),1.0f)
            descriptor(im1b_mind,im1b,m,n,o,mind_step[level]);
            gettimeofday(&time2, NULL);
            timeMIND+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
		}
		
		int len3=pow(hw1*2+1,3);
		int m1=m/step1; int n1=n/step1; int o1=o/step1; int sz1=m1*n1*o1;
        
        float* costall=new float[sz1*len3]; float* costall2=new float[sz1*len3];
		
        //cout<<"==========================================================\n";
		//cout<<"Level "<<level<<" grid="<<step1<<" hw="<<hw1<<" quant="<<quant1<<"\n";
		//cout<<"==========================================================\n";
		
        //FULL-REGISTRATION FORWARDS
        gettimeofday(&time1, NULL);
        gettimeofday(&time2, NULL);
		timeTrans+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"T"<<flush;
        gettimeofday(&time1, NULL);
		descriptor(warped_mind,warped1,m,n,o,mind_step[level]);

        gettimeofday(&time2, NULL);
		timeMIND+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"M"<<flush;
        gettimeofday(&time1, NULL);
        dataCostCL((unsigned long*)im1b_mind,(unsigned long*)warped_mind,costall,m,n,o,len3,step1,hw1,quant1,alpha,RAND_SAMPLES);
        gettimeofday(&time2, NULL);

		timeData+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"D"<<flush;
        gettimeofday(&time1, NULL);
        gettimeofday(&time2, NULL);
		timeSmooth+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"S"<<flush;
		
        //FULL-REGISTRATION BACKWARDS
        gettimeofday(&time1, NULL);
		gettimeofday(&time2, NULL);
		timeTrans+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"T"<<flush;
        gettimeofday(&time1, NULL);
		descriptor(warped_mind,warped2,m,n,o,mind_step[level]);

        gettimeofday(&time2, NULL);
		timeMIND+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"M"<<flush;
        gettimeofday(&time1, NULL);
        dataCostCL((unsigned long*)im1_mind,(unsigned long*)warped_mind,costall2,m,n,o,len3,step1,hw1,quant1,alpha,RAND_SAMPLES);
        gettimeofday(&time2, NULL);
		timeData+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"DS\n"<<flush;
        gettimeofday(&time1, NULL);
        estimateAffine2(X,Xprev,im1b,im1,costall,costall2,step1,quant1,hw1);

        gettimeofday(&time2, NULL);
		timeSmooth+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        //cout<<"S"<<flush;
		
        printf("t: MIND=%2.2f, data=%2.2f, affine=%2.2f, speed=%2.2e dof/s\n",timeMIND,timeData,timeSmooth,2.0*(float)sz1*(float)len3/(timeData+timeSmooth));
        
        gettimeofday(&time1, NULL);

        gettimeofday(&time2, NULL);
        float timeMapping=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
		
        
        for(int i=0;i<16;i++){
            Xprev[i]=X[i];
        }

        
        timeDataSmooth+=(timeSmooth+timeData+timeMIND+timeTrans);
        
        delete costall; delete costall2;
		
	}
    delete im1_mind;
    delete im1b_mind;
	//==========================================================================================
	//==========================================================================================
	
    gettimeofday(&time2a, NULL);
	float timeALL=time2a.tv_sec+time2a.tv_usec/1e6-(time1a.tv_sec+time1a.tv_usec/1e6);
    
    
    string outputfile;
    outputfile.append(args.output_stem);
    outputfile.append("_matrix.txt");

    ofstream matfile;
    matfile.open(outputfile);
    
    for(int i=0;i<4;i++){
        matfile<<X[i]<<"  "<<X[i+4]<<"  "<<X[i+8]<<"  "<<X[i+12]<<"\n";
        
        printf("%+4.3f | %+4.3f | %+4.3f | %+4.3f \n",X[i],X[i+4],X[i+8],X[i+12]);
    }
    matfile.close();
    
    // if SEGMENTATION of moving image is provided APPLY SAME TRANSFORM
    if(args.segment){
        short* seg2;
        readNifti(args.moving_seg_file,seg2,header,M,N,O,P);
        
        float* zero=new float[sz];
        fill(zero,zero+sz,0.0f);

        short* segw=new short[sz];
        fill(segw,segw+sz,0);
        
        warpAffineS(segw,seg2,X,zero,zero,zero);
        
        
        string outputseg;
        outputseg.append(args.output_stem);
        outputseg.append("_deformed_seg.nii.gz");
        
        
        
        gzWriteSegment(outputseg,segw,header,m,n,o,1);
    }
  
    
	cout<<"Finished. Total time: "<<timeALL<<" sec. ("<<timeDataSmooth<<" sec. for MIND+data+affine+trans)\n";
	
	
	return 0;
}
