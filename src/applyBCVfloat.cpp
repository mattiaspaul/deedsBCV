#include <iostream>
#include <fstream>
#include <math.h>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include <string.h>
#include <map>
#include <sstream>
#include <x86intrin.h>
#include <pthread.h>
#include <thread>
#include "zlib.h"
#include <sys/stat.h>

using namespace std;

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



struct parameters{
    float alpha; int levels=0; bool segment,affine,rigid;
    vector<int> grid_spacing; vector<int> search_radius;
    vector<int> quantisation;
    string fixed_file,moving_file,output_stem,moving_seg_file,affine_file,deformed_file;
};

#include "imageIOgzType.h"
#include "transformations.h"
//#include "primsMST.h"
//#include "regularisation.h"
//#include "MINDSSCbox.h"
#include "dataCostD.h"
#include "parseArguments.h"


int main (int argc, char * const argv[]) {
	//Initialise random variable

    
    //PARSE INPUT ARGUMENTS
    
    if(argc<4||argv[1][1]=='h'){
        cout<<"=============================================================\n";
        cout<<"Usage (required input arguments):\n";
        cout<<"./applyBCVfloat -M moving.nii.gz -O output -D deformed.nii.gz \n";
        cout<<"optional parameters:\n";
        cout<<" -A <affine_matrix.txt> \n";
        cout<<"=============================================================\n";
        return 1;
    }

    parameters args;
    parseCommandLine(args, argc, argv);
    
    size_t split_def=args.deformed_file.find_last_of("/\\");
    if(split_def==string::npos){
        split_def=-1;
    }
    size_t split_moving=args.moving_file.find_last_of("/\\");
    if(split_moving==string::npos){
        split_moving=-1;
    }
    
    
    if(args.deformed_file.substr(args.deformed_file.length()-2)!="gz"){
        cout<<"images must have nii.gz format\n";
        return -1;
    }
    if(args.moving_file.substr(args.moving_file.length()-2)!="gz"){
        cout<<"images must have nii.gz format\n";
        return -1;
    }
    
    cout<<"Transforming "<<args.moving_file.substr(split_moving+1)<<" into "<<args.deformed_file.substr(split_def+1)<<"\n";
    

    
    
    //READ IMAGES and INITIALISE ARRAYS
    
    timeval time1,time2,time1a,time2a;
    

    float* img2;
    int M,N,O,P; //image dimensions
    
    //==ALWAYS ALLOCATE MEMORY FOR HEADER ===/
    char* header=new char[352];
    
    readNifti(args.moving_file,img2,header,M,N,O,P);
    
    image_m=M; image_n=N; image_o=O;
    
    int m=image_m; int n=image_n; int o=image_o; int sz=m*n*o;
    
    
    
    //READ AFFINE MATRIX from linearBCV if provided (else start from identity)
    
    float* X=new float[16];
    
    if(args.affine){
        size_t split_affine=args.affine_file.find_last_of("/\\");
        if(split_affine==string::npos){
            split_affine=-1;
        }
        
        cout<<"Reading affine matrix file: "<<args.affine_file.substr(split_affine+1)<<"\n";
        ifstream matfile;
        matfile.open(args.affine_file);
        for(int i=0;i<4;i++){
            string line;
            getline(matfile,line);
            sscanf(line.c_str(),"%f  %f  %f  %f",&X[i],&X[i+4],&X[i+8],&X[i+12]);
        }
        matfile.close();
        
        
    }
    else{
        cout<<"Using identity transform.\n";
        fill(X,X+16,0.0f);
        X[0]=1.0f; X[1+4]=1.0f; X[2+8]=1.0f; X[3+12]=1.0f;
    }
    
    for(int i=0;i<4;i++){
        printf("%+4.3f | %+4.3f | %+4.3f | %+4.3f \n",X[i],X[i+4],X[i+8],X[i+12]);//X[i],X[i+4],X[i+8],X[i+12]);
        
    }

    
    
    
    string inputflow;
    inputflow.append(args.output_stem);
    inputflow.append("_displacements.dat");
    
    cout<<"Reading displacements from:\n"<<inputflow<<"\n";
    
    
    cout<<"=============================================================\n";
    
    vector<float> flow=readFile<float>(inputflow);
    
    int sz3=flow.size()/3;
    int grid_step=round(pow((float)sz/(float)sz3,0.3333333));

    cout<<"grid step "<<grid_step<<"\n";
    
    int step1; int hw1; float quant1;
    
    //set initial flow-fields to 0; i indicates backward (inverse) transform
    //u is in x-direction (2nd dimension), v in y-direction (1st dim) and w in z-direction (3rd dim)
    float* ux=new float[sz]; float* vx=new float[sz]; float* wx=new float[sz];
    for(int i=0;i<sz;i++){
        ux[i]=0.0; vx[i]=0.0; wx[i]=0.0;
    }
    
    int m1,n1,o1,sz1;
    m1=m/grid_step; n1=n/grid_step; o1=o/grid_step; sz1=m1*n1*o1;
    float* u1=new float[sz1]; float* v1=new float[sz1]; float* w1=new float[sz1];
   
    for(int i=0;i<sz1;i++){
        u1[i]=flow[i]; v1[i]=flow[i+sz1]; w1[i]=flow[i+sz1*2];
        
    }

    upsampleDeformationsCL(ux,vx,wx,u1,v1,w1,m,n,o,m1,n1,o1);


    
    float* warped=new float[sz];
    
    fill(warped,warped+sz,(short)0);

    warpAffine(warped,img2,img2,X,ux,vx,wx);

    
    
    gzWriteNifti(args.deformed_file,warped,header,m,n,o,1);

    
	
	
	return 0;
}
