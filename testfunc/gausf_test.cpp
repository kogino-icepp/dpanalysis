#include <iostream>
#include "../headers/setting.h"
#include <TROOT.h>

Double_t gausf(double x,double p0,double mean,double sigma){
    return p0*TMath::Exp(-(x-mean)*(x-mean)/(2*sigma*sigma));
}
void gausf_test(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    Setting st;
    st.dot_size = 0.8;
    st.markerstyle = 20;
    st.color = kBlue;
    axrange ax = {-5,5,0,3,0,1};
    TF1* f_test = new TF1("f_test","gausf(x,[0],[1],[2])",-5,5);
    //TF1* f_test = new TF1("f_test","chiF_ndf2(x,1,[0],0.1)",0,10);
    f_test -> SetParameter(0,1);
    f_test -> SetParameter(1,0);
    f_test -> SetParameter(2,1);
    TGraph* graph = new TGraph(f_test);
    st.Graph(graph,ax);
    graph -> Draw("AP");
}