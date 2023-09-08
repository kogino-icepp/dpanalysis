#include <iostream>
#include <queue>
#include "../headers/fitter.h"
const int nbin=32767;
const double c= 3.0*pow(10,8);
const double v0=220000.0;
const double vE=v0;
const double vc=v0;
const double dnu=0.000076296;//ビン幅
const double dNu=0.0000885;//周波数分解能
const double kb=1.38*pow(10,-23);
const double df=88.5*pow(10,3);
const double Tc=76;
const double Th=297;


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
double F_sig1(double f,double f0,double P,double r){
    double rtn = P*(F_nu(f+r*dNu,f0)-F_nu(f-r*dNu,f0));
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
void signal_test(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    double dfreq = 2.5/nbin;
    double fmin = 223-10*dfreq;
    double fmax = 223+15*dfreq;
    double Psig = 1.5*pow(10,-14);
    TH1D* hist = new TH1D("hist","test_signal;Freq[GHz];P[W/kHz]",26,fmin,fmax);
    TF1* sigf = new TF1("sigf","F_sig2(x,[0],[1],0.5)",222,224);
    sigf -> SetLineWidth(5);
    sigf -> FixParameter(0,223);
    sigf -> FixParameter(1,Psig);
    for(int i=0;i<30;i++){
        double freq = 223+(i-15)*dfreq;
        hist -> Fill(freq,F_sig2(freq,223,Psig,0.5));
    }
    st.Hist(hist);
    hist -> SetFillColor(kCyan);
    hist -> SetLineColor(kBlack);
    hist->SetMarkerStyle(1); // ビンのマーカースタイルを設定 (1 は普通の点)
    hist->SetMarkerColor(kBlack); // ビンのマーカーの色を黒に設定
    hist->SetLineColor(kBlack); // ビンの境界線の色を黒に設定
    hist -> SetLineWidth(2);
    hist -> Draw("HIST");
    sigf -> Draw("same");
}