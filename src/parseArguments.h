parameters parseCommandLine(int argc, char * const argv[]){
    
    typedef pair<char,int> val;
    map<char,int> argin;
    argin.insert(val('F',0));
    argin.insert(val('M',1));
    argin.insert(val('O',2));
    argin.insert(val('a',3));
    argin.insert(val('l',4));
    argin.insert(val('G',5));
    argin.insert(val('L',6));
    argin.insert(val('Q',7));
    argin.insert(val('S',8));
    argin.insert(val('A',9));

    // parsing the input
    int requiredArgs=0;
    char* fixedfile=new char[200];
    char* movingfile=new char[200];
    char* outputstem=new char[200];
    char* movsegfile=new char[200];
    char* affinefile=new char[200];

    float alpha=1.6;
    int maxlevel=5;
    int num=maxlevel; int num2=maxlevel; int num3=maxlevel;
    bool s_set=false;
    int s_grid[10]={8,7,6,5,4,3,2,2,2,2};
    int s_search[10]={8,7,6,5,4,3,2,1,1};
    int s_quant[10]={5,4,3,2,1,1,1,1,1,1};
    char levelstr[]="%dx%dx%dx%dx%dx%dx%dx%dx%dx%d";
    
    bool symmetric=true;
    bool segment=false;
    bool affine=false;

    for(int k=1;k<argc;k++){
        if(argv[k][0]=='-'){
            if(argin.find(argv[k][1])==argin.end()){
                cout<<"Invalid option: "<<argv[k]<<" use -h for help\n";
                }
                
                switch(argin[argv[k][1]]){
                    case 0:
                        sprintf(fixedfile,"%s",argv[k+1]);
                        requiredArgs++;
                        break;
                    case 1:
                        sprintf(movingfile,"%s",argv[k+1]);
                        requiredArgs++;
                        break;
                    case 2:
                        sprintf(outputstem,"%s",argv[k+1]);
                        requiredArgs++;
                        break;
                    case 3:
                        alpha=atof(argv[k+1]);
                        break;
                    case 4:
                        maxlevel=atoi(argv[k+1]);
                        break;
                    case 5:
                        num=sscanf(argv[k+1],levelstr,&s_grid[0],&s_grid[1],&s_grid[2],&s_grid[3],&s_grid[4],&s_grid[5],&s_grid[6],&s_grid[7],&s_grid[8],&s_grid[9]);
                        s_set=true;
                        break;
                    case 6:
                        num2=sscanf(argv[k+1],levelstr,&s_search[0],&s_search[1],&s_search[2],&s_search[3],&s_search[4],&s_search[5],&s_search[6],&s_search[7],&s_search[8],&s_search[9]);
                        s_set=true;
                        break;
                    case 7:
                        num3=sscanf(argv[k+1],levelstr,&s_quant[0],&s_quant[1],&s_quant[2],&s_quant[3],&s_quant[4],&s_quant[5],&s_quant[6],&s_quant[7],&s_quant[8],&s_quant[9]);
                        s_set=true;
                        break;
                    case 8:
                        sprintf(movsegfile,"%s",argv[k+1]);
                        segment=true;
                        break;
                    case 9:
                        sprintf(affinefile,"%s",argv[k+1]);
                        affine=true;
                        break;
                    default:
                        cout<<"Invalid option: "<<argv[k]<<" use -h for help\n";
                        break;
                }
        }
    }
    
    
    if(requiredArgs!=3){
        cout<<"Missing argmuents, use -h for help.\n";
    }
    if(s_set){
        maxlevel=num;
    }
    
    if((num!=maxlevel)|(num2!=maxlevel)|(num3!=maxlevel)){
        cout<<"Max level and number of grid-spacing, search range\n or quantisation steps are not equal.\n";
        printf("maxlevel=%d, #grid=%d, #search=%d, #quant=%d\n",maxlevel,num,num2,num3);
    }
    
    //copy parameters into struct (and convert to strings/vectors)
    struct parameters args;
    string s1(fixedfile);
    args.fixed_file=s1;
    string s2(movingfile);
    args.moving_file=s2;
    string s3(outputstem);
    args.output_stem=s3;
    string s4(movsegfile);
    args.moving_seg_file=s4;
    string s5(affinefile);
    args.affine_file=s5;

    args.alpha=alpha;
    args.levels=maxlevel;
    args.segment=segment;
    args.affine=affine;
    
    for(int i=0;i<maxlevel;i++){
        args.grid_spacing.push_back(s_grid[i]);
        args.search_radius.push_back(s_search[i]);
        args.quantisation.push_back(s_quant[i]);
        
    }
    
    return args;
    
}




                //==========================================================================================
                //==========================================================================================
                //IMPORTANT SETTINGS FOR CONTROL POINT SPACING AND LABEL SPACE
                //int maxlevel=4;
                
                //int label_hw[]={6,5,4,2};//
                //half-width of search space L={±0,±1,..,label_hw}^3 * label_quant
                
                //int grid_step[]={8,6,4,2};//
                //spacing between control points in grid
                //d.o.f.: 2.8, 3.4, 4.3, 6.1(3.4)
                
                //float label_quant[]={3,2,1,1};//
                //quantisation of search space L, important: can only by integer or 0.5 so far!
                
                /*
                 float* mind_step=new float[maxlevel];
                 //[]={3,2,2,1,1};//{1,1,1,1,1};//
                 for(int i=0;i<maxlevel;i++){
                 mind_step[i]=ceil(label_quant[i]*0.49999);
                 }
                 
                 //optionally write-out warped Labels
                 if(segment){
                 for(int i=0;i<sz;i++){
                 ux[i]=round(ux[i]);
                 vx[i]=round(vx[i]);
                 wx[i]=round(wx[i]);
                 }
                 short *seg;
                 readNiftiShort(movingsegin,seg,M,N,O,header);
                 short* warpedseg=new short[sz];
                 warpImage(warpedseg,seg,ux,vx,wx);
                 
                 writeNiftiShort(output3,warpedseg,header,sz);
                 }
                 
                 writeOutput(flow1,output2,sz1*3);
                 */
                
