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

void warpShort(short* segw,short* seg1,float* X,int m,int n,int o){
    int sz=m*n*o;
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                
                float y1=(float)i*X[0]+(float)j*X[1]+(float)k*X[2]+(float)X[3];
                float x1=(float)i*X[4]+(float)j*X[5]+(float)k*X[6]+(float)X[7];
                float z1=(float)i*X[8]+(float)j*X[9]+(float)k*X[10]+(float)X[11];
                int x=round(x1); int y=round(y1);  int z=round(z1);
                
                if((y)>=0&&(y)<m&&(x)>=0&&(x)<n&&(z)>=0&&(z)<o){
                    segw[i+j*m+k*m*n]=seg1[y+(x)*m+(z)*m*n];
                }
                else{
                    segw[i+j*m+k*m*n]=0;
                }
                
            }
        }
    }
}



void warpAffine(float* warped,float* input,float* X,int m,int n,int o){
    int m2=m; int n2=n; int o2=o;
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                
                float y1=(float)i*X[0]+(float)j*X[1]+(float)k*X[2]+(float)X[3];
                float x1=(float)i*X[4]+(float)j*X[5]+(float)k*X[6]+(float)X[7];
                float z1=(float)i*X[8]+(float)j*X[9]+(float)k*X[10]+(float)X[11];
                int x=floor(x1); int y=floor(y1);  int z=floor(z1);
                float dx=x1-x; float dy=y1-y; float dz=z1-z;
                
                if(y<0|y>=m2-1|x<0|x>=n2-1|z<0|z>=o2-1){
                    warped[i+j*m+k*m*n]=0.0;
                }
                else{
                    warped[i+j*m+k*m*n]=(1.0-dx)*(1.0-dy)*(1.0-dz)*input[min(max(y,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
                    (1.0-dx)*dy*(1.0-dz)*input[min(max(y+1,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
                    dx*(1.0-dy)*(1.0-dz)*input[min(max(y,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
                    (1.0-dx)*(1.0-dy)*dz*input[min(max(y,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2]+
                    dx*dy*(1.0-dz)*input[min(max(y+1,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
                    (1.0-dx)*dy*dz*input[min(max(y+1,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2]+
                    dx*(1.0-dy)*dz*input[min(max(y,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2]+
                    dx*dy*dz*input[min(max(y+1,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2];
                }
            }
        }
    }
    
    
}



void estimateAffine2(float* X,float* Xprev,float* im1,float* im2,float* costall,float* costall2,int step1,float quant1,int hw1){
    
    int m=image_m;
    int n=image_n;
    int o=image_o;
    
    int m1=m/step1; int n1=n/step1; int o1=o/step1; int sz1=m1*n1*o1;
    
    int len=hw1*2+1; int len3=len*len*len;
    int stephw=(step1-1)/2;
    float* xs=new float[len*len*len];
	float* ys=new float[len*len*len];
	float* zs=new float[len*len*len];
    
    
	for(int i=0;i<len;i++){
		for(int j=0;j<len;j++){
			for(int k=0;k<len;k++){
				xs[i+j*len+k*len*len]=(j-hw1)*quant1;
				ys[i+j*len+k*len*len]=(i-hw1)*quant1;
				zs[i+j*len+k*len*len]=(k-hw1)*quant1;
			}
		}
	}
    int nx[]={-1,1,0,0,0,0};
    int ny[]={0,0,-1,1,0,0};
    int nz[]={0,0,0,0,-1,1};
    // insert here: smoothing of datacost (linear b-spline)
    
    for(int i=0;i<sz1;i++){
        int z1=i/(m1*n1); int x1=(i-z1*m1*n1)/m1; int y1=i-z1*m1*n1-x1*m1;
        int i1=min(y1+1,m1-1)+x1*m1+z1*m1*n1;
        int i2=max(y1-1,0)+x1*m1+z1*m1*n1;
        for(int l=0;l<len3;l++){
            costall[i*len3+l]+=(0.5*costall[i1*len3+l]+0.5*costall[i2*len3+l]);
            costall2[i*len3+l]+=(0.5*costall2[i1*len3+l]+0.5*costall2[i2*len3+l]);
        }
    }
    for(int i=0;i<sz1;i++){
        int z1=i/(m1*n1); int x1=(i-z1*m1*n1)/m1; int y1=i-z1*m1*n1-x1*m1;
        int i1=y1+min(x1+1,n1-1)*m1+z1*m1*n1;
        int i2=y1+max(x1-1,0)*m1+z1*m1*n1;
        for(int l=0;l<len3;l++){
            costall[i*len3+l]+=(0.5*costall[i1*len3+l]+0.5*costall[i2*len3+l]);
            costall2[i*len3+l]+=(0.5*costall2[i1*len3+l]+0.5*costall2[i2*len3+l]);
        }
    }
    for(int i=0;i<sz1;i++){
        int z1=i/(m1*n1); int x1=(i-z1*m1*n1)/m1; int y1=i-z1*m1*n1-x1*m1;
        int i1=y1+x1*m1+min(z1+1,o1-1)*m1*n1;
        int i2=y1+x1*m1+max(z1-1,0)*m1*n1;
        for(int l=0;l<len3;l++){
            costall[i*len3+l]+=(0.5*costall[i1*len3+l]+0.5*costall[i2*len3+l]);
            costall2[i*len3+l]+=(0.5*costall2[i1*len3+l]+0.5*costall2[i2*len3+l]);
        }
    }
    
    
    float* minval=new float[sz1];
    int* minind=new int[sz1];

    float* minval2=new float[sz1];
    int* minind2=new int[sz1];

    int allcount=0;
    for(int i=0;i<sz1;i++){
        minval[i]=*min_element(costall+i*len3,costall+(i+1)*len3);
        minind[i]=min_element(costall+i*len3,costall+(i+1)*len3)-(costall+i*len3);

        minval2[i]=*min_element(costall2+i*len3,costall2+(i+1)*len3);
        minind2[i]=min_element(costall2+i*len3,costall2+(i+1)*len3)-(costall2+i*len3);

        int z1=i/(m1*n1); int x1=(i-z1*m1*n1)/m1; int y1=i-z1*m1*n1-x1*m1;
        int y2=y1*step1+stephw; int x2=x1*step1+stephw; int z2=z1*step1+stephw;
        
        float meanval1=0.0;
        float meanval2=0.0;
        for(int k1=0;k1<step1;k1++){
            for(int j1=0;j1<step1;j1++){
                for(int i1=0;i1<step1;i1++){
                    int ind1=y1*step1+i1+(x1*step1+j1)*m+(z1*step1+k1)*m*n;
                    meanval1+=im1[ind1];
                    meanval2+=im2[ind1];
                }
            }
        }
        meanval1/=(float)(step1*step1*step1);
        meanval2/=(float)(step1*step1*step1);
        
        ///removing (=assigning high value) to control points from background and/or boundaries might help ...
        //if(im1[y2+x2*m+z2*m*n]<10.0f)
        if(meanval1<5.0f){
            minval[i]=1e20;
        }
        if(meanval2<5.0f){
            minval2[i]=1e20;
        }
        
        if(z1<2|z1>o1-2|x1<2|x1>n1-2|y1<2|y1>m1-2){
            minval[i]=1e20;
            minval2[i]=1e20;
        }
        if(minval[i]<1e19)
            allcount++;
        if(minval2[i]<1e19)
            allcount++;
        
    }
    float median1=quantile(minval,sz1,0.5*allcount/(float)sz1/2.0f);
    float median2=quantile(minval2,sz1,0.5*allcount/(float)sz1/2.0f);
    
    //vector<int> validind;
    //vector<int> validind2;
    int valsz=0;
    
    for(int i=0;i<sz1;i++){
        if(minval[i]<median1){
            valsz++;
        }
        if(minval2[i]<median2){
            valsz++;
        }
    }
      printf("# points used: %d/%d, quantile: %f\n",valsz,allcount,0.5*allcount/(float)sz1/2.0f);
    
    float* pts1=new float[valsz*4];
    float* pts2=new float[valsz*4];

    int current=0;
    for(int i=0;i<sz1;i++){
        if(minval[i]<median1){
            int ind=minind[i];
            pts1[current+3*valsz]=1.0f; pts2[current+3*valsz]=1.0f;
            int k1=i/(m1*n1); int j1=(i-k1*(m1*n1))/m1; int i1=i-j1*m1-k1*m1*n1;
            int y1=i1*step1+stephw; int x1=j1*step1+stephw; int z1=k1*step1+stephw;
            pts1[current]=y1; pts1[current+valsz]=x1; pts1[current+valsz*2]=z1;
            pts2[current]=y1+ys[ind]; pts2[current+valsz]=x1+xs[ind]; pts2[current+valsz*2]=z1+zs[ind];
            current++;
        }
        if(minval2[i]<median2){
            int ind=minind2[i];
            pts1[current+3*valsz]=1.0f; pts2[current+3*valsz]=1.0f;
            int k1=i/(m1*n1); int j1=(i-k1*(m1*n1))/m1; int i1=i-j1*m1-k1*m1*n1;
            int y1=i1*step1+stephw; int x1=j1*step1+stephw; int z1=k1*step1+stephw;
            pts2[current]=y1; pts2[current+valsz]=x1; pts2[current+valsz*2]=z1;
            pts1[current]=y1+ys[ind]; pts1[current+valsz]=x1+xs[ind]; pts1[current+valsz*2]=z1+zs[ind];
            current++;
        }
    }

    float X1[16];
    affineRobust(X1,pts1,pts2,valsz);

    matmult(X1,Xprev,X);

  
    delete minval; delete minval2;
    delete minind; delete minind2;
    delete pts1; delete pts2;
    delete xs; delete ys; delete zs;
    
    
   // estimateAffine(X,Xprev,costall,validind,costall2,validind2,step1,quant1,hw1);
}
