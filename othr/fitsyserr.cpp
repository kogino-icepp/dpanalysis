#include <iostream>
#include <random>
#include "../headers/fitter.h"
using namespace std;
const double c= 3.0*pow(10,8);
const double v0=220000.0;
const double vE=v0;
const double vc=v0;
const double dnu=0.000076296;//ビン幅[GHz]
const double dNu=0.0000885;//周波数分解能[GHz]
const double kb=1.38*pow(10,-23);
const double df=88.5*pow(10,3);
const double Tc=76;
const double Th=297;
double candpow[6] = {0,200,400,600,800,1000};
double candmass[3] = {220,240,260};
double v_conv(double f,double f0){
    double rtn=c*sqrt(1-((f0/f)*(f0/f)));
    return rtn;
}
double F_nu(double f,double f0){
    double rtn;
    double v=v_conv(f,f0);
    double p_kata=(v+vE)/v0;
    double m_kata=(v-vE)/v0;
    rtn = (vc/(2*sqrt(M_PI)*vE))*(exp(-(p_kata*p_kata))-exp(-(m_kata*m_kata)));
    rtn += 0.5*(erf(p_kata)+erf(m_kata));
    //rtn *= (c*f0*f0)/(f*f*f*sqrt(1-(f0/f)*(f0/f)));
    return rtn;
}
double F_sig2(double f,double f0,double P,double r){
    if(f+dnu*r<=f0)return 0;
    else if(f+r*dNu>f0 && f-(1-r)*dNu<=f0){
        return P*(F_nu(f+r*dNu,f0)-F_nu(f0,f0));
    }
    else if(f-dnu*(1-r)>f0){
        return P*(F_nu(f+r*dNu,f0)-F_nu(f-(1-r)*dNu,f0));
    }
    else return 0;
}
void fitsyserr(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);

    Setting st;
    
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_int_distribution<> rand100(0, 99);
    normal_distribution<double> normrand(0,100);
    //m_γ=220,240,260GHzにおいて周波数確度が-1~1kHzの不定性を持っている時のフィット精度を検証する
    double f0 = 220;//DM質量に対応する周波数[GHz]
    TH1* testhist = new TH1D("testhist","test;Freq[GHz];",20,f0-dnu,f0+19*dnu);
    //TH1* testhist2 = new TH1D("testhist2","test;Freq[GHz];",20,f0-dnu,f0+19*dnu);//実験用のヒストグラム、質量をちょっとずらす
    //testhist2 -> SetLineColor(kRed);
    rep(m,6){
        double P = 0;
        double Perr = 0;
        rep(ite,5000){
            TH1* fithist = new TH1D("fithist","",20,f0-dnu,f0+19*dnu);
            rep(bin,20){
                double freq = 220+(bin-1)*dnu;
                double perr = normrand(mt);
                fithist -> SetBinContent(bin,F_sig2(freq,220,1000,0.5)+perr);
                fithist -> SetBinError(bin,100);
            }
            delete fithist;
        }
    }
    /*rep(bin,20){
        double freq = 220+(bin-1)*dnu;
        
        double perr = normrand(mt);
        testhist -> SetBinContent(bin,F_sig2(freq,220,1000,0.5)+perr);
        testhist -> SetBinError(bin,100);
        //testhist2 -> SetBinContent(bin,F_sig2(freq,220-1*pow(10,-6),1000,0.5));
    }*/
    testhist -> SetFillColor(kBlue);
    st.Hist(testhist);
    testhist -> Draw("E");
    TF1* pfunc = new TF1("pfunc","F_sig2(x,220,[0],0.5)",f0-dnu,f0+19*dnu);
    testhist -> Fit(pfunc);
    //このフィットを5000回回せばええの？
}