#include <iostream>
#include "../headers/fitter.h"

using namespace std;
void logtest(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.12,0.14,0.2,0.1);
    Setting st;
    st.dot_size=0.8;
    st.markerstyle=20;
    st.color = kGreen;
    TGraph* graph = new TGraph;
    TF1* area = new TF1("area","0.5",216,264);
    //c1 -> SetLogy();
    for(int i=17;i<18;i++){
        TF1* gftest = new TF1("gftest","TMath::GammaDist(x*[0]*[1],[0],0,1)*[2]",0,10);
        gftest -> SetParameter(0,i);
        gftest -> SetParameter(1,1);
        gftest -> SetParameter(2,10);
        gftest -> Draw();
        //delete gftest;
    }
    
}