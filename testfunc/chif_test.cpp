#include <iostream>
#include "../headers/setting.h"
#include <TROOT.h>
vector<double> errorv = {0.0622,0.0646,0.0682,0.0809,0.0666,0.0865,0.0658,0.0612,0.0754,0.0618,0.0873,0.0634,0.0940,0.151,0.0859,0.102,0.0877,0.114,0.136,0.0898,0.0964,0.0799,0.108,0.0991};
Double_t chif(double x,double p0,double p1){
    return p0*TMath::Exp(-x/2)*pow(x,(p1/2)-1)/(pow(2,p1/2)*TMath::Gamma(p1/2));
}
Double_t chif_ndf(double x,double p0,double p1){
    return p1*chif(p1*x,p0,p1);
}
Double_t chiF(double x,double p0,double p1){
    return p0*TMath::Gamma(p1/2,x/2);
}
Double_t chiF_ndf(double x,double p0,double p1){
    return chiF(p1*x,p0,p1);
}
Double_t chiF_ndf2(double x,double p0,double p1,double bin){
    return (chiF_ndf((x+bin/2),p0,27)-chiF_ndf((x-bin/2),p0,27));
}
Double_t chif_free(double x,double p0,double k,double p1){
    //e^-x/2*x^(k/2-1)/2^(k/2)*Gamma(k/2) y=p1x
    return p0*TMath::Exp(-x/(2*p1))*pow(x,(k/2)-1)/(pow(2*p1,k/2)*TMath::Gamma(k/2));
}
Double_t chiF_free(double x,double p0,double k,double p1){
    return p0*TMath::Gamma(k/2,x/(2*p1));
}
Double_t chiF_freefit(double x,double p0,double p1,double k,double bin){
    return (chiF_free((x+bin/2),p0,k,p1)-chiF_free((x-bin/2),p0,k,p1));
}
//正規化されたchi^2/NDFにおけるp-valueを算出する
Double_t p_value(double x,double k){
    return 1-chiF_free(x,1,k,1/k);
}
/*
inc_gamma_c() = 1/Gamma int_x^inf
*/
void chif_test(){
    
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    Setting st;
    st.dot_size = 0.8;
    st.markerstyle = 20;
    st.color = kBlue;
    TF1* frui = new TF1("frui","chiF_free(x,26280,17,0.03)");
    //TF1* ffra = new TF1("ffra","chiF_freefit(x,26280,0,03,17,0.)")
    TGraph* graph = new TGraph();
    //ヒストグラムだから多分これで良い
    double bhaba = 0.1;
    
    TF1* f_test = new TF1("f_test","chiF_freefit(x,[0],[1],27,0.1)",0,10);
    vector<int> fpara = {1,2,3,4,27};
    vector<Color_t> color = {kBlack,kBlue,kGreen,kRed,kMagenta};
    
    
}