
float norm(float* vector,int len){
    float n=0;
    for(int i=0;i<len;i++){
        n+=pow(vector[i],2);
    }
    return sqrt(n);
}

float dotprod(float* vector,float* vector2,int len){
    float dot=0;
    for(int i=0;i<len;i++){
        dot+=vector[i]*vector2[i];
    }
    return dot;
}

void qrsolve(float* X,float* A,float* b,int len,int len2){
    //first column
    float* Q=new float[4*len];
    float r11=norm(A,len);
    for(int i=0;i<len;i++){
        Q[i]=A[i]/r11;
    }
    float r12=dotprod(Q,A+len,len);
    float r13=dotprod(Q,A+len*2,len);
    float r14=dotprod(Q,A+len*3,len);
    //printf("r11: %f, r12: %f, r13: %f, r14: %f\n",r11,r12,r13,r14);
    
    //second column - misuse Q(:,2:4) for a2,a3,a4
    for(int i=0;i<len;i++){
        Q[i+len]=A[i+len]-Q[i]*r12;
        Q[i+len*2]=A[i+len*2]-Q[i]*r13;
        Q[i+len*3]=A[i+len*3]-Q[i]*r14;
    }
    float r22=norm(Q+len,len);
    for(int i=0;i<len;i++){
        Q[i+len]=Q[i+len]/r22;
    }
    float r23=dotprod(Q+len,A+len*2,len);
    float r24=dotprod(Q+len,A+len*3,len);
    //printf("r22: %f, r23: %f, r24: %f\n",r22,r23,r24);
    
    //third column
    for(int i=0;i<len;i++){
        Q[i+len*2]=A[i+len*2]-Q[i]*r13-Q[i+len]*r23;
        Q[i+len*3]=A[i+len*3]-Q[i]*r14-Q[i+len]*r24;
    }
    float r33=norm(Q+len*2,len);
    for(int i=0;i<len;i++){
        Q[i+len*2]=Q[i+len*2]/r33;
    }
    float r34=dotprod(Q+len*2,A+len*3,len);
    
    //fourth column
    for(int i=0;i<len;i++){
        Q[i+len*3]=A[i+len*3]-Q[i]*r14-Q[i+len]*r24-Q[i+2*len]*r34;
    }
    float r44=norm(Q+len*3,len);
    for(int i=0;i<len;i++){
        Q[i+len*3]=Q[i+len*3]/r44;
    }
    //printf("r33: %f, r34: %f, r44: %f\n",r33,r34,r44);
    for(int j=0;j<len2;j++){
        float d1=dotprod(Q,b+j*len,len);
        float d2=dotprod(Q+len,b+j*len,len);
        float d3=dotprod(Q+len*2,b+j*len,len);
        float d4=dotprod(Q+len*3,b+j*len,len);
        float x4=d4/r44;
        float x3=(d3-r34*x4)/r33;
        float x2=(d2-r23*x3-r24*x4)/r22;
        float x1=(d1-r12*x2-r13*x3-r14*x4)/r11;
        X[0+j*4]=x1; X[1+j*4]=x2; X[2+j*4]=x3; X[3+j*4]=x4;
        //printf("x1: %f, x2: %f, x3: %f, x4: %f\n",x1,x2,x3,x4);
    }
    
}

void jacobiSVD3(float* A,float* U,float* V){
    for(int i=0;i<9;i++){
        U[i]=A[i];
        V[i]=0;
    }
    float T[3];
    
    V[0]=1; V[4]=1; V[8]=1;
    for(int iter=0;iter<4;iter++){
        for(int j=1;j<3;j++){
            for(int i=0;i<j;i++){
                float alpha=dotprod(U+i*3,U+i*3,3);
                float beta=dotprod(U+j*3,U+j*3,3);
                float gamma=dotprod(U+i*3,U+j*3,3);
                float zeta=(beta-alpha)/(2.0*gamma);
                float t=((zeta>0)?1.0:-1.0)/(fabs(zeta)+sqrt(1.0+zeta*zeta));
                float c=1.0/sqrt(1+t*t);
                float s=c*t;
                T[0]=U[0+i*3]; T[1]=U[1+i*3]; T[2]=U[2+i*3];
                for(int k=0;k<3;k++){
                    U[k+i*3]=c*T[k]-s*U[k+j*3];
                }
                for(int k=0;k<3;k++){
                    U[k+j*3]=s*T[k]+c*U[k+j*3];
                }
                T[0]=V[0+i*3]; T[1]=V[1+i*3]; T[2]=V[2+i*3];
                for(int k=0;k<3;k++){
                    V[k+i*3]=c*T[k]-s*V[k+j*3];
                }
                for(int k=0;k<3;k++){
                    V[k+j*3]=s*T[k]+c*V[k+j*3];
                }
				
            }
        }
    }
    
    for(int j=0;j<3;j++){
        float singval=norm(U+j*3,3);
        for(int i=0;i<3;i++){
            U[i+j*3]/=singval;
        }
    }
    
	
    
}

void findRigid(float* RT,float* pts1t,float* pts2t,int len){
	
	//if needed transpose
	float* pts1=new float[len*4];
	float* pts2=new float[len*4];
	for(int i=0;i<4;i++){
		for(int k=0;k<len;k++){
			pts1[i+k*4]=pts1t[k+i*len];
			pts2[i+k*4]=pts2t[k+i*len];
		}
	}
	
		
	
	float* pts1m=new float[3];
	float* pts2m=new float[3];
	for(int i=0;i<3;i++){
		pts1m[i]=0.0; pts2m[i]=0.0;
	}
	for(int k=0;k<len;k++){
		pts1m[0]+=pts1[0+k*4]; pts1m[1]+=pts1[1+k*4]; pts1m[2]+=pts1[2+k*4];
		pts2m[0]+=pts2[0+k*4]; pts2m[1]+=pts2[1+k*4]; pts2m[2]+=pts2[2+k*4];
	}
	for(int i=0;i<3;i++){
		pts1m[i]/=(float)len; pts2m[i]/=(float)len;
	}
	
	float* K=new float[9];
	float* U=new float[9];
	float* V=new float[9];
	float* R=new float[16];
	for(int i=0;i<9;i++){
		K[i]=0.0; U[i]=0.0; V[i]=0.0;
	}

	for(int k=0;k<len;k++){
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				K[i+j*3]+=(pts1[j+k*4]-pts1m[j])*(pts2[i+k*4]-pts2m[i]);
			}
		}
	}
	
	
    jacobiSVD3(K,U,V);
    
    for(int i=0;i<16;i++){
        R[i]=0.0; RT[i]=0.0;
    }
    R[15]=1.0;
    RT[15]=1.0;
    
    for(int k=0;k<3;k++){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                R[j+i*4]+=U[j+k*3]*V[i+k*3];
            }
        }
    }
   // float* pts2_4=new float[len*4];
    float* Rpts2=new float[len*4];
    /*for(int i=0;i<len;i++){
        pts2_4[0+i*4]=pts2[0+i*3]; pts2_4[1+i*4]=pts2[1+i*3];
        pts2_4[2+i*4]=pts2[2+i*3]; pts2_4[3+i*4]=1.0f;
    }*/
    
    qrsolve(Rpts2,R,pts2,4,len);
    
	//    printf("Rpt: %f, %f, %f\n",Rpts2[0],Rpts2[len],Rpts2[2*len]);
	
    
    float* t1=new float[3];
	t1[0]=0; t1[1]=0; t1[2]=0;
    for(int i=0;i<len;i++){
        t1[0]+=(Rpts2[i*4]);//-pts1[i*3]);
        t1[1]+=(Rpts2[1+i*4]);//-pts1[1+i*3]);
        t1[2]+=(Rpts2[2+i*4]);//-pts1[2+i*3]);
    }
    for(int i=0;i<3;i++){
        t1[i]=t1[i]/(float)len-pts1m[i];
    }
	// printf("t1: %f, %f, %f\n",t1[0],t1[1],t1[2]);
	
    float t0[3]={0,0,0};
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            t0[i]+=R[i+j*4]*t1[j];
        }
    }
    
	// printf("t0: %f, %f, %f\n",t0[0],t0[1],t0[2]);
    
    //transpose
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            RT[i+j*4]+=R[j+i*4];
        }
        RT[3+i*4]=t0[i];
    }
	
	delete pts1;
	delete pts2;
    
   // delete pts2_4;
    delete Rpts2;
	
}



void affineRobust(float* RT,float* pts1,float* pts2,int len){
	
	if(RIGID)
        findRigid(RT,pts1,pts2,len);
    else
        qrsolve(RT,pts1,pts2,len,4);
	
	
	int lenh=len/2;
	float* err=new float[len];
	float* err2=new float[len];
	float* pts1b=new float[lenh*4];
	float* pts2b=new float[lenh*4];
	for(int i=0;i<lenh;i++){
		pts1b[i+0*lenh]=1; pts2b[i+0*lenh]=1;
		pts1b[i+1*lenh]=1; pts2b[i+1*lenh]=1;
		pts1b[i+2*lenh]=1; pts2b[i+2*lenh]=1;
		pts1b[i+3*lenh]=1; pts2b[i+3*lenh]=1;
	}
	
	for(int iter=0;iter<15;iter++){
		for(int l=0;l<len;l++){
			float x=0; float y=0; float z=0;
			for(int i=0;i<3;i++){
				y+=pts1[l+i*len]*RT[i];
				x+=pts1[l+i*len]*RT[i+4];
				z+=pts1[l+i*len]*RT[i+8];
			}
			y+=pts1[l+3*len]*RT[3]; x+=pts1[l+3*len]*RT[3+4]; z+=pts1[l+3*len]*RT[3+8];
			err[l]=pow(y-pts2[l],2)+pow(x-pts2[l+len],2)+pow(z-pts2[l+len*2],2);
			err2[l]=err[l];
		}
		nth_element(err,err+len/2,err+len);
		float median=err[len/2];
		int count=0;
		for(int i=0;i<len;i++){
			int l=min(count,lenh-1);
			if(err2[i]<=median){
				pts1b[l]=pts1[i];
				pts1b[l+lenh]=pts1[i+len];
				pts1b[l+2*lenh]=pts1[i+2*len];
				pts2b[l]=pts2[i];
				pts2b[l+lenh]=pts2[i+len];
				pts2b[l+2*lenh]=pts2[i+2*len];
				count++;
			}
		}
        
        if(RIGID)
            findRigid(RT,pts1b,pts2b,lenh);
        else
            qrsolve(RT,pts1b,pts2b,lenh,4);
        
	
        
		
		
	}	
	delete pts1b; delete pts2b; delete err; delete err2;
	
}