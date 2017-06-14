/*
 Mattias P. Heinrich
 Universitaet Luebeck, 2014
 */
template <typename Type>

void readNifti(string filestr,Type*& vol,char*& header,int& m,int& n,int& o,int& p){
    
    //convert input string into char* array
    char* filename=new char[filestr.size()+1];
    copy(filestr.begin(),filestr.end(),filename);
    filename[filestr.size()]='\0';
    //ifstream file(filename,ios::in|ios::binary);
    
    gzFile file=gzopen(filename,"rb");
    if(!file){
        printf("gzopen of '%s' failed.\n",filename);
        exit(-1);
    }
    
    
    gzread(file,header,352);
    
    short* dimensions;

    dimensions=reinterpret_cast<short*>(header+40);
    m=(int)dimensions[1];
    n=(int)dimensions[2];
    o=(int)dimensions[3];
    p=(int)dimensions[4];
    
    vol=new Type[m*n*o*p];

    //printf("Read image with dimensions %dx%dx%dx%d\n",m,n,o,p);
    short* datatype; //read datatype (and bitpix)
    datatype=reinterpret_cast<short*>(header+70);
    short bitpix=datatype[1];
    //
    char* filecharptr=new char[m*n*o*p*64];
    //raw empty datapointers of all 'types'
    unsigned char* ucharptr; float* floatptr; short* shortptr; int* intptr; double* doubleptr;
    //read input data values depending on datatype and convert to float
    //binary are read character by character, reinterpret_cast converts them
    //copy them into float* array afterwards
    switch(datatype[0]){
        case 2:
            gzread(file,filecharptr,m*n*o*p*sizeof(unsigned char));
            ucharptr=reinterpret_cast<unsigned char*>(filecharptr);
            for(int i=0;i<m*n*o*p;i++){
                vol[i]=ucharptr[i];
            }
            break;
        case 4:
            gzread(file,filecharptr,m*n*o*p*sizeof(short));
            shortptr=reinterpret_cast<short*>(filecharptr);
            for(int i=0;i<m*n*o*p;i++){
                vol[i]=shortptr[i];
            }
            break;
        case 8:
            gzread(file,filecharptr,m*n*o*p*sizeof(int));
            intptr=reinterpret_cast<int*>(filecharptr);
            for(int i=0;i<m*n*o*p;i++){
                vol[i]=intptr[i];
            }
            break;
        case 16:
            gzread(file,filecharptr,m*n*o*p*sizeof(float));
            floatptr=reinterpret_cast<float*>(filecharptr);
            for(int i=0;i<m*n*o*p;i++){
                vol[i]=floatptr[i];
            }
            break;
        case 64:
            gzread(file,filecharptr,m*n*o*p*sizeof(double));
            doubleptr=reinterpret_cast<double*>(filecharptr);
            for(int i=0;i<m*n*o*p;i++){
                vol[i]=doubleptr[i];
            }
            break;
        default:
            printf("Datatype %d not supported. Exiting.\n",datatype[0]);
            exit(1);
    }
    delete filecharptr;
    int sz=m*n*o*p;
    
    //gzread(file,vol,sz*sizeof(Type));
    //vol=reinterpret_cast<Type*>(vol);
    
    gzclose(file);
    
}


void writeNifti(string filestr,float* vol,char* header,int m,int n,int o,int k){
    //convert input string into char* array
    char* filename=new char[filestr.size()+1];
    copy(filestr.begin(),filestr.end(),filename);
    filename[filestr.size()]='\0';
    
    char* header2=new char[352];
    copy(header,header+352,header2);
    int dim=o>1?3:2;
    if(k>1){
        dim=4;
    }
    short* dimensions=new short[4];
    dimensions[0]=dim; dimensions[1]=m; dimensions[2]=n; dimensions[3]=o;
    dimensions[4]=k;
    char* dimchar=reinterpret_cast<char*>(dimensions);
    copy(dimchar,dimchar+8,header2+40);

    short* datatype=new short[2];
    datatype[0]=16; datatype[1]=32; //float datatype
    char* datachar=reinterpret_cast<char*>(datatype);
    copy(datachar,datachar+4,header2+70);
    
    //printf("Writing image with dimensions %dx%dx%d\n",m,n,o);

    ofstream file(filename,ios::out|ios::binary);
    //opens file for binary-output
	if(file.is_open()){
		file.write(header2,352);
		file.write(reinterpret_cast<char*>(vol),m*n*o*4);
		file.close();
		cout<<"File "<<filename<<" written.\n";
	}
    else{
        printf("File error. Could not write file.\n");
    }
}

void writeSegment(string filestr,short* vol,char* header,int m,int n,int o){
    //convert input string into char* array
    char* filename=new char[filestr.size()+1];
    copy(filestr.begin(),filestr.end(),filename);
    filename[filestr.size()]='\0';
    
    char* header2=new char[352];
    copy(header,header+352,header2);
    int dim=o>1?3:2;
    short* dimensions=new short[4];
    dimensions[0]=dim; dimensions[1]=m; dimensions[2]=n; dimensions[3]=o;
    char* dimchar=reinterpret_cast<char*>(dimensions);
    copy(dimchar,dimchar+8,header2+40);
    
    short* datatype=new short[2];
    datatype[0]=4; datatype[1]=512; //short datatype
    char* datachar=reinterpret_cast<char*>(datatype);
    copy(datachar,datachar+4,header2+70);
    
    //printf("Writing segmentation with dimensions %dx%dx%d\n",m,n,o);
    
    ofstream file(filename,ios::out|ios::binary);
    //opens file for binary-output
	if(file.is_open()){
		file.write(header2,352);
		file.write(reinterpret_cast<char*>(vol),m*n*o*2);
		file.close();
		cout<<"File "<<filename<<" written.\n";
	}
    else{
        printf("File error. Could not write file.\n");
    }
}


void gzWriteNifti(string filestr,float* data,char* header,int m,int n,int o,int p){
    
    int sz=m*n*o*p;
    //convert input string into char* array
    char* filename=new char[filestr.size()+1];
    copy(filestr.begin(),filestr.end(),filename);
    filename[filestr.size()]='\0';
    
    //do not change original header, define (new )dimensions and datatype
    char* header2=new char[352];
    copy(header,header+352,header2);
    int dim=o>1?3:2;
    if(p>1){ dim=4;}
    short* dimensions=new short[4];
    dimensions[0]=dim; dimensions[1]=m; dimensions[2]=n;
    dimensions[3]=o; dimensions[4]=p;
    char* dimchar=reinterpret_cast<char*>(dimensions);
    copy(dimchar,dimchar+10,header2+40);
    short* datatype=new short[2];
    datatype[0]=16; datatype[1]=32; //float datatype
    char* typechar=reinterpret_cast<char*>(datatype);
    copy(typechar,typechar+4,header2+70);
    
    //copy header and data into out
    int size=352+sz*sizeof(float);
    char* out2=new char[size];
    copy(header2,header2+352,out2);
    char* datachar=reinterpret_cast<char*>(data);
    copy(datachar,datachar+sz*sizeof(float),out2+352);
    
    //write gzipped file
    gzFile file2=gzopen(filename,"wb");
    if(!file2){
        printf("File error. Could not write file.\n");
        exit(1);
    }
    int err=gzwrite(file2,(unsigned char*)out2,size);
    gzclose(file2);
    
    delete out2; delete header2;
}

void gzWriteSegment(string filestr,short* data,char* header,int m,int n,int o,int p){
    int sz=m*n*o*p;
    //convert input string into char* array
    char* filename=new char[filestr.size()+1];
    copy(filestr.begin(),filestr.end(),filename);
    filename[filestr.size()]='\0';
    
    //do not change original header, define (new )dimensions and datatype
    char* header2=new char[352];
    copy(header,header+352,header2);
    int dim=o>1?3:2;
    if(p>1){ dim=4;}
    short* dimensions=new short[4];
    dimensions[0]=dim; dimensions[1]=m; dimensions[2]=n;
    dimensions[3]=o; dimensions[4]=p;
    char* dimchar=reinterpret_cast<char*>(dimensions);
    copy(dimchar,dimchar+10,header2+40);
    short* datatype=new short[2];
    datatype[0]=4; datatype[1]=16;// edit 08-16 - 512; //short datatype
    char* typechar=reinterpret_cast<char*>(datatype);
    copy(typechar,typechar+4,header2+70);
    
    //copy header and data into out
    int size=352+sz*sizeof(short);
    char* out2=new char[size];
    copy(header2,header2+352,out2);
    char* datachar=reinterpret_cast<char*>(data);
    copy(datachar,datachar+sz*sizeof(short),out2+352);
    
    //write gzipped file
    gzFile file2=gzopen(filename,"wb");
    if(!file2){
        printf("File error. Could not write file.\n");
        exit(1);
    }
    int err=gzwrite(file2,(unsigned char*)out2,size);
    gzclose(file2);
    
    delete out2; delete header2;
    
}


/*
void writeSegment(string filestr,short* vol,char* header,int m,int n,int o){
    //convert input string into char* array
    char* filename=new char[filestr.size()+1];
    copy(filestr.begin(),filestr.end(),filename);
    filename[filestr.size()]='\0';
    
    //create/change header part -> header2
    char* header2=new char[352];
    copy(header,header+352,header2);
    int dim=o>1?3:2;
    short* dimensions=new short[4];
    dimensions[0]=dim; dimensions[1]=m; dimensions[2]=n; dimensions[3]=o;
    char* dimchar=reinterpret_cast<char*>(dimensions);
    copy(dimchar,dimchar+8,header2+40);
    short* datatype=new short[2];
    datatype[0]=4; datatype[1]=512; //short datatype
    char* datatypechar=reinterpret_cast<char*>(datatype);
    copy(datatypechar,datatypechar+4,header2+70);
    
    int sz=m*n*o;
    
    printf("Writing segmentation with dimensions %dx%dx%d\n",m,n,o);
    
    int size=sz*sizeof(short)+352; //total size
    //unsigned long destLen = compressBound(size); // this is how you should estimate size

    //copy header and data into one char array
    char* data=new char[size];
    copy(header2,header2+352,data);
    char* datainchar=reinterpret_cast<char*>(vol);
    copy(datainchar,datainchar+sz*sizeof(short),data+352);
    
    
    printf("prepared header and data\n");
    
    gzFile file=gzopen(filename,"wb");
    if(!file){
        printf("gzopen of '%s' failed.\n",filename);
        exit(-1);
    }
    int err=gzwrite(file,(unsigned char*)data,size);
    gzclose(file);
    
    int sized=1000;
    char* outd=new char[sized];
    
    gzFile file2=gzopen("compressed2.nii.gz","wb");
    int err2=gzwrite(file2,(unsigned char*)outd,sized);
    gzclose(file2);
    
    //delete data;
    
}*/

void readFloat(char str2[], float*& pixels,int SZ){
	FILE * pFile;
	int i,j,offset;
	size_t result;
	
	pFile = fopen (str2, "rb" );
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	
	fseek(pFile,0,0);
	
	float* memory=new float[SZ];
	size_t lSize;
	float * buffer;
	lSize=SZ;
	buffer = (float*) malloc (sizeof(float)*lSize);
	if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}
	
	// copy the file into the buffer:
	result = fread (buffer,sizeof(float),lSize,pFile);
	if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
	
	for(i=0;i<(SZ);i++)
		pixels[i]=buffer[i];
	free (buffer);
	fclose (pFile);
}
void writeOutput(float* data,char* name,int length){
	ofstream ofs1(name,ofstream::binary);
	ofs1.write((char *) data,length*sizeof(float));
	ofs1.close();
}


void writeOutputI(int* data,char* name,int length){
	ofstream ofs1(name,ofstream::binary);
	ofs1.write((char *) data,length*sizeof(int));
	ofs1.close();
}

void writeOutputS(short* data,char* name,int length){
	
	ofstream ofs1(name,ofstream::binary);
	ofs1.write((char *) data,length*sizeof(short));
	ofs1.close();
	
	
	
}


