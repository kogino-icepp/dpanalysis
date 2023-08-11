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
    TGraph* graph = new TGraph;
    //ヒストグラムだから多分これで良い
    double bhaba = 0.1;
    /*TF1* test_func = new TF1("test_func","chiF_ndf2(x,[0],[1],[2])",0,10);
    test_func -> SetParameter(0,7000);
    test_func -> SetParameter(1,27);
    test_func -> SetParameter(2,0.1);
    for(int i=0;i<100;i++){
        //graph -> SetPoint(i,i*bhaba+0.05,(chif_ndf((i+1)*bhaba,1,)-chif_ndf((i)*bhaba,1,27))/0.1);
        graph -> SetPoint(i,i*bhaba+0.05,(chiF_ndf((i+1)*bhaba,1,27)-chiF_ndf((i)*bhaba,1,27)));
        cout << fixed;
        cout << setprecision(10) << endl;
        cout << "x : " << i*bhaba << " y: " << (chiF((i)*bhaba,1,2)) << endl;
    }
    axrange ax = {0,10,0,3000,0,1,"test;test;test"};
    st.Graph(graph,ax);
    c1 -> SetLogy();
    graph -> Draw("AL");
    test_func -> Draw("same");*/
    TF1* f_test = new TF1("f_test","chiF_freefit(x,1,[1],[0],0.1)",0,10);
    //TF1* f_test = new TF1("f_test","chiF_ndf2(x,1,[0],0.1)",0,10);
    vector<int> fpara = {1,2,3,4,27};
    vector<Color_t> color = {kBlack,kBlue,kGreen,kRed,kMagenta};
    rep(i,5){
        f_test -> SetLineColor(color[i]);
        f_test -> SetParameter(0,fpara[i]);
        double p2 = 1.0/fpara[i];
        f_test -> SetParameter(1,p2);
        TGraph* graph = new TGraph(f_test);
        if(i==0){
            axrange ax = {0,10,0,1};
            st.Graph(graph,ax);
            graph -> Draw();
        }
        else graph -> Draw("same");
    }
    //f_test -> Draw();
}