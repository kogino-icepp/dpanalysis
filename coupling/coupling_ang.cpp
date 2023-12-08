#include <iostream>
#include <TROOT.h>
using namespace std;
double cv=3.0*pow(10,11); 
int freq=240; 
double lam=cv/(freq*pow(10,9)); 
double w0=2.67;
double z0=M_PI*w0*w0/lam;
double k=2*M_PI/lam;
double west(double z){
    double w=w0*sqrt(1+pow(z/z0,2));
    return w;
}
double R(double z){
    double r=0;
    r=z+(1/z)*pow((M_PI*w0*w0/lam),2);
    return r;
}
double eff(double z){
    double e=0;
    e=sqrt(2/(M_PI*west(z)*west(z)));
    return e;
}
double xrot1(double x,double z,double theta){
    double rtn=0;
    rtn=-sin(theta)*z+cos(theta)*x;
    return rtn;
}
double zrot1(double x,double z,double theta){
    double rtn=0;
    rtn=cos(theta)*z+sin(theta)*x;
    return rtn;
}
complex<double> gauss(double x,double y, double z){
    complex<double> kata(-(x*x+y*y)/(west(z)*west(z)),-k*z-M_PI*(x*x+y*y)/(lam*R(z)));
    complex<double> g=exp(kata);
    return g;
}
//,double dx,double dz
complex<double> eta(double theta){
    complex<double> eta(0);
    for(int i=-350;i<351;i++){
        for(int j=-350;j<351;j++){
            double z=sqrt(1500*1500-i*i-j*j);
            //double inpro=0;
            //inpro=()/(1500*1500)
            eta+=eff(z)*eff(zrot1(i,z,theta))*gauss(i,j,z)*conj(gauss(xrot1(i,z,theta),j,zrot1(i,z,theta)));
        }
    }
    return eta;
}
void coupling_ang()
{
    const Int_t n=200;
    Double_t x[n]={0};
    Double_t y[n]={0};
    for(int i=-25;i<25;i++){
        x[i+50]=i*0.02;
        y[i+50]=pow(abs(eta(x[i+50])),4);
    }
    TGraph *graph=new TGraph(n,x,y);
    graph->SetMarkerColor(4);
    graph->SetMarkerStyle(4);
    graph->Draw("AP");
    graph->SetTitle("1theta;theta[rad];coupling");
}