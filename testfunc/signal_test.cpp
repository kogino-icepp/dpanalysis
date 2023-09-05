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
}
void signal_test(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    TF1* sigf = new TF1("sigf","F_sig2(x,[0],[1],0.5)",223,224);
    sigf -> SetParameter(0,223);
    TGraph* graph = new TGraph;
    double dfreq = 2.5/nbin;
    for(int i=0;i<30;i++){
        double freq = 223+(i-15)*dfreq;
        graph -> SetPoint(i,freq,F_sig2(freq,223,10,0.5));
    }
    double fmin = 223-10*dfreq;
    double fmax = 223+15*dfreq;
    axrange axtest = {fmin,fmax,-1,5,0,1,"test_signal;Freq[GHz];Prec[K]"};
    st.Graph(graph,axtest);
    graph -> SetFillColor(kCyan);
    graph -> Draw("AB");

}