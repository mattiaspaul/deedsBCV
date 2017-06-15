/*
 Mattias P. Heinrich
 Universitaet Luebeck, 2014
 */
#include <iostream>     // std::cout
#include <algorithm>    // std::count
#include <fstream>       // filein
#include <numeric>   // std::accumulate
#include <vector>       // std::vector
#include "zlib.h"

using namespace std;

#include "imageIOgzType.h"


int main(int argc,char * const argv[]){
    
    if(argc<3){
        cout<<"Not enough input arguments. Example syntax:\n";
        cout<<"./diceMulti seg1.nii seg2.nii\n";
        exit(1);
    }
    
    string filein1=argv[1];
    string filein2=argv[2];
    
    char* header=new char[352];
    short* segment1;
    int m,n,o,p; //dimensions of image (will be set be readNifti)
    readNifti(filein1,segment1,header,m,n,o,p); //load segmentation1
    short* segment2;
    readNifti(filein2,segment2,header,m,n,o,p); //load segmentation2
    
    int sz=m*n*o;
    
    vector<float> dice1;

    bool* countvol=new bool[sz];
    
    vector<short> found_label(segment1,segment1+sz);
    sort(found_label.begin(),found_label.end());
    found_label.erase(unique(found_label.begin(),found_label.end()),found_label.end());
    
    int maxlabel=*max_element(segment1,segment1+sz);
    for(int i=1;i<found_label.size();i++){ //for all labels except background
        int label=found_label[i];
        for(int i=0;i<sz;i++)
            countvol[i]=segment1[i]==label&segment2[i]==label;
        float unionval=count(countvol,countvol+sz,true);
        for(int i=0;i<sz;i++)
            countvol[i]=segment1[i]==label;
        float card1=count(countvol,countvol+sz,true);
        for(int i=0;i<sz;i++)
            countvol[i]=segment2[i]==label;
        float card2=count(countvol,countvol+sz,true);
        dice1.push_back(2.0f*unionval/(card1+card2));
        
    }

    cout<<"Dice (for each label): ";
    for(int i=0;i<dice1.size();i++){
        printf(" %4.3f, ",dice1[i]);
    }
    float meandice=accumulate(dice1.begin(),dice1.end(),0.0f)/dice1.size();
    
    printf("\nAVERAGE DICE: %f\n",meandice);
    /*if(argc>3){
        ofstream outfile;
        outfile.open(argv[3], std::ios_base::app);
        outfile<<argv[4]<<" to "<<argv[5];//<<":\t"<<meandice<<"\n";
        for(int i=0;i<dice1.size();i++){
            outfile<<"\t"<<dice1[i];
        }
        outfile<<"\n";
    }*/
    delete countvol;
    return 0;
}

