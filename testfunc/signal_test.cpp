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
const double Aeff = 0.385;
const double rho = 0.3*1.602*pow(10,-2);//DMの密度(をSI単位系に直したもの)

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
    else if(f+r*dnu>f0 && f-(1-r)*dnu<=f0){
        return P*(F_nu(f+r*dnu,f0)-F_nu(f0,f0));
    }
    else if(f-dnu*(1-r)>f0){
        return P*(F_nu(f+r*dnu,f0)-F_nu(f-(1-r)*dnu,f0));
    }
    else return 0;
}
double normchi(int k,double x){
    return pow(k/2,(k/2))*pow(x,k/2-1)*TMath::Exp(-x*k/(2))/TMath::Gamma(k/2);
}
double CoupConst(double p,double dp){
    if(p<0)return sqrt(3*1.96*dp/(2*Aeff*rho));
    else return sqrt(3*(p+1.96*dp)/(2*Aeff*rho));
}
double CoupConst2(double p,double dp){
    if(p<0)return 4.5*pow(10,-14)*sqrt(dp*1.96*pow(10,23))*sqrt(1/Aeff);
    else return 4.5*pow(10,-14)*sqrt((dp*1.96+p)*pow(10,23))*sqrt(1/Aeff);
}
//もしかして理想的な分布から出てくるp値とフィットのp-valueって違う？？
double PValue(int k,double x){
    return 1-TMath::Gamma(k/2,(k/2)*x);
}
void CalculateCumulativeDistribution() {
    // ガウシアン関数のパラメータ
    double mean = 0.0;    // 平均
    double sigma = 1.0;   // 標準偏差

    // ガウシアン関数を定義
    TF1 *gaussian = new TF1("gaussian", "TMath::Gaus(x, [0], [1])", -10.0, 10.0);
    gaussian->SetParameters(mean, sigma);

    // 累積分布を計算
    TF1 *cumulative = new TF1("cumulative", "gaussian->Integral(-5, x)", -5.0, 5.0);

    
}


void signal_test(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    double dfreq = 2.5/nbin;
    double fmin = 223-10*dfreq;
    double fmax = 223+15*dfreq;
    double Psig = pow(10,-14);
    TH1D* hist = new TH1D("hist","test_signal;Freq[GHz];P[W/kHz]",26,fmin,fmax);
    TF1* sigf = new TF1("sigf","F_sig2(x,223,[0],0.5)",222,224);
    
    sigf -> SetLineWidth(5);
    TGraph* graph = new TGraph;
    for(int i=0;i<30;i++){
        double freq = 223+(i-15)*dfreq;
        graph -> SetPoint(i,freq,F_sig2(freq,223,Psig,0.5));
        hist -> Fill(freq,F_sig2(freq,223,Psig,0.5));
    }
    sigf -> SetParameter(0,Psig*1.015);
    /*cout << dfreq << endl;
    axrange ax = {fmin,fmax,0,5*pow(10,-15),0,1,"TestSignal;Freq[GHz];P[W/Hz]"};
    st.Graph(graph,ax);
    graph -> SetFillColor(kCyan);
    graph -> SetLineWidth(2);
    graph -> SetLineColor(kBlack);
    graph -> Draw("AB");
    sigf -> Draw("same");*/
    //TF1* testf = new TF1("testf","normchi(25,x)",0,10);
    //testf -> Draw();
    st.Hist(hist);
    hist -> SetLineWidth(5);
    hist -> SetLineColor(kBlack);
    hist -> SetFillColor(kCyan);
    hist -> Draw("HIST");
    sigf -> Draw("same");
    cout << PValue(26,3) << endl;
    cout << CoupConst2(-1,2*1.358*3.8*pow(10,-19)) << endl;
    //とりあえず比で誤魔化しているけどここは明らかにしないといけないパート
}