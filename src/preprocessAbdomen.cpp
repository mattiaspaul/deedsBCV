#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <sys/time.h>
#include <math.h>
#include <inttypes.h>

using namespace std;


#include "zlib.h"
#include "imageIOgzType.h"

//#include "symmetricDiffeomorphic.h"

void matmult(float* A,float* B,float *C){
    
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            C[i+j*4]=0.0;
            for(int k=0;k<4;k++){
                C[i+j*4]+=A[i+k*4]*B[k+j*4];
                
            }
        }
    }
}
template <typename Type>
void warpAffine(float* warped,Type* input,float* X,int m,int n,int o,int m_orig,int n_orig,int o_orig){
    int m2=m_orig; int n2=n_orig; int o2=o_orig;
    float min_val=*min_element(input,input+m2*n2*o2);
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                
                float y1=(float)i*X[0]+(float)j*X[1]+(float)k*X[2]+(float)X[3];
                float x1=(float)i*X[4]+(float)j*X[5]+(float)k*X[6]+(float)X[7];
                float z1=(float)i*X[8]+(float)j*X[9]+(float)k*X[10]+(float)X[11];
                int x=floor(x1); int y=floor(y1);  int z=floor(z1);
                float dx=x1-x; float dy=y1-y; float dz=z1-z;
                
                if(y<0|y>=m2-1|x<0|x>=n2-1|z<0|z>=o2-1){
                    warped[i+j*m+k*m*n]=min_val;
                }
                else{
                    warped[i+j*m+k*m*n]=(1.0-dx)*(1.0-dy)*(1.0-dz)*(float)input[min(max(y,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
                    (1.0-dx)*dy*(1.0-dz)*(float)input[min(max(y+1,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
                    dx*(1.0-dy)*(1.0-dz)*(float)input[min(max(y,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
                    (1.0-dx)*(1.0-dy)*dz*(float)input[min(max(y,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2]+
                    dx*dy*(1.0-dz)*(float)input[min(max(y+1,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
                    (1.0-dx)*dy*dz*(float)input[min(max(y+1,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2]+
                    dx*(1.0-dy)*dz*(float)input[min(max(y,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2]+
                    dx*dy*dz*(float)input[min(max(y+1,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2];
                }
            }
        }
    }
    
    
}

template <typename Type>
void warpAffineL(short* warped,Type* input,float* X,int m,int n,int o,int m_orig,int n_orig,int o_orig){
    int m2=m_orig; int n2=n_orig; int o2=o_orig;
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                
                float y1=(float)i*X[0]+(float)j*X[1]+(float)k*X[2]+(float)X[3];
                float x1=(float)i*X[4]+(float)j*X[5]+(float)k*X[6]+(float)X[7];
                float z1=(float)i*X[8]+(float)j*X[9]+(float)k*X[10]+(float)X[11];
                int x=round(x1); int y=round(y1);  int z=round(z1);
                
                if(y<0|y>=m2-1|x<0|x>=n2-1|z<0|z>=o2-1){
                    warped[i+j*m+k*m*n]=0;
                }
                else{
                        warped[i+j*m+k*m*n]=input[y+x*m2+z*m2*n2];
                    
                }
            }
        }
    }
    
    
}

template <typename TypeReflect>

void reflect(TypeReflect*& vol,int m,int n,int o){
    int sz=m*n*o;
    //swap x and y-axes
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m/2;i++){
                swap(vol[i+j*m+k*m*n],vol[(m-1-i)+j*m+k*m*n]);
            }
        }
    }
    for(int k=0;k<o;k++){
        for(int j=0;j<n/2;j++){
            for(int i=0;i<m;i++){
                swap(vol[i+j*m+k*m*n],vol[i+(n-1-j)*m+k*m*n]);
            }
        }
    }
}
/*void reflect(short*& vol,int m,int n,int o){
    int sz=m*n*o; int sz2=sz/2;
    for(int i=0;i<sz2;i++){
        swap(vol[i],vol[sz-1-i]);
    }
}*/

int main(int argc, const char * argv[])
{
    
    if(argc<3){
        cout<<"=============================================================\n";
        cout<<"Usage (required input arguments either two or four):\n";
        cout<<" ./preprocess orig_moving_im.nii.gz moving_im.nii.gz orig_moving_seg.nii.gz moving_seg.nii.gz\n";
        cout<<"  ./preprocess orig_target_im.nii.gz target_im.nii.gz\n";
        cout<<"#pre-processing (requires segmentations to be uint8 and scans to be short)\n";
        return -1;
    }
    
    //int casenum=atoi(argv[5]);
    
    string filescan="temp0.nii.gz";
    
    short* im_orig;
    int M,N,O,K;
    char* header_orig=new char[352];
    
    readNifti(argv[1],im_orig,header_orig,M,N,O,K);
    
    printf("M=%d, N=%d, O=%d\n",M,N,O);
    float* vox; float voxel[3];
    vox=reinterpret_cast<float*>(header_orig+76);
    voxel[0]=(float)vox[1]; voxel[1]=(float)vox[2]; voxel[2]=(float)vox[3];
    
    int m=180; int n=140; int o=190; int sz=m*n*o;

    float* resampled=new float[sz];
    float voxelnew=2.2;//atof(argv[4]);
    float voxel2[3]={voxelnew,voxelnew,voxelnew};
    float affdiag[3]={voxel2[0]/voxel[0],voxel2[1]/voxel[1],voxel2[2]/voxel[2]};
    
    float translate[3]={(float)M/2.0f-((float)m/2.0f*affdiag[0]),(float)N/2.0f-((float)n/2.0f*affdiag[1]),
        (float)O/2.0f-((float)o/2.0f*affdiag[2])};
    
    //for CT only
    int transz=20;//atoi(argv[3]);
    translate[2]+=transz; cout<<"CT translate "<<transz<<" (z-axis)\n";
    
    float T[16]={affdiag[0],0,0,translate[0], 0,affdiag[1],0,translate[1], 0,0,affdiag[2],translate[2], 0,0,0,1};
    
    for(int i=0;i<4;i++){
        printf("%+4.3f | %+4.3f | %+4.3f | %+4.3f \n",T[i],T[i+4],T[i+8],T[i+12]);
    }
    //affine matrix
    //T0=diag([[2.5,2.5,2.5]./pixdim,1])
//    T(4,1:4)=(-[centre2,1]*T0+[centre1,2])
    
    //invT=diag([1./T0,1]); invT(4,1:3)=([centre2]-[centre1]./T0)
    
    warpAffine(resampled,im_orig,T,m,n,o,M,N,O);

    char* header=new char[352];
    copy(header_orig,header_orig+352,header);
    
    //set values for dimensions datatype and voxelsize
    vector<short> dim(8,1);
    dim[0]=o>1?3:2; dim[1]=m; dim[2]=n; dim[3]=o;
    vector<short> datatype(2,16); datatype[1]=32; //float-datatype 32 bit
    vector<float> pixdim(8,1); pixdim[0]=0;
    pixdim[1]=voxel2[0]; pixdim[2]=voxel2[1]; pixdim[3]=voxel2[2];
    
    memcpy(header+76,pixdim.data(),8*sizeof(float));
    memcpy(header+40,dim.data(),8*sizeof(short));
    memcpy(header+70,datatype.data(),2*sizeof(short));

    gzWriteNifti(argv[2],resampled,header,m,n,o,1);
    cout<<"intensity scan resampled\n";
    
    
    
        //   string segfile="temp";
      //  segfile.append(to_string(l));
      //  segfile.append(".nii.gz");
    if(argc==5){
        short* resampled_seg=new short[sz];
        for(int i=0;i<sz;i++){
            resampled_seg[i]=0;
        }
        
        char* segim1;
        readNifti(argv[3],segim1,header_orig,M,N,O,K);
        warpAffineL(resampled_seg,segim1,T,m,n,o,M,N,O);
        
        
        datatype[0]=4; datatype[1]=16; //short-datatype 16 bit
        memcpy(header+70,datatype.data(),2*sizeof(short));
        
        gzWriteSegment(argv[4],resampled_seg,header,m,n,o,1);
        cout<<"segmentation resampled\n";

    }
    
    return 0;
}


