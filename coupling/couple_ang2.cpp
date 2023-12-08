#include <iostream>
#include "../headers/fitter.h"
using namespace std;

double cv=3.0*pow(10,11); 
int freq=240; 
double lam=cv/(freq*pow(10,9)); 
double w0=2.67;
double z0=M_PI*w0*w0/lam;
double k=2*M_PI/lam;
double R0=1500;
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
complex<double> gauss(double x,double y, double z){
    complex<double> kata(-(x*x+y*y)/(west(z)*west(z)),-k*z-M_PI*(x*x+y*y)/(lam*R(z)));
    complex<double> g=exp(kata);
    return g;
}
double zrot(double z,double x,double theta){
    double rtn=0;
    rtn=cos(theta)*(z+R0*cos(theta)-R0)+sin(theta)*(x+R0*sin(theta));
    return rtn;
}
double xrot(double z,double x,double theta){
    double rtn=0;
    rtn=-sin(theta)*(z+R0*cos(theta)-R0)+cos(theta)*(x+R0*sin(theta));
    return rtn;
}
complex<double> eta(double theta){
    complex<double> eta(0);
    for(int i=-350;i<351;i++){
        for(int j=-350;j<351;j++){
            if(i*i+j*j<=350*350){
                double z=sqrt(1500*1500-i*i-j*j);
                double norm1=1500;
                double norm2=sqrt(xrot(z,i,theta)*xrot(z,i,theta)+zrot(z,i,theta)*zrot(z,i,theta)+j*j);
                double pd=(xrot(z,i,theta)*i+j*j+zrot(z,i,theta)*z)/(norm1*norm2);
                eta+=pd*eff(z)*eff(zrot(z,i,theta))*conj(gauss(i,j,z))*gauss(xrot(z,i,theta),j,zrot(z,i,theta));
            }            
        }
    }
    return eta;
}
complex<double> eta_2d(double dx,double theta,double dz){
    complex<double> eta(0);
    for(int i=-350;i<351;i++){
        for(int j=-350;j<351;j++){
            if(i*i+j*j<=350*350){
                double z=sqrt(1500*1500-i*i-j*j);
                double norm1=1500;
                double norm2=sqrt(xrot(z-dz,i-dx,theta)*xrot(z-dz,i-dx,theta)+zrot(z-dz,i-dx,theta)*zrot(z-dz,i-dx,theta)+j*j);
                double pd=(xrot(z-dz,i-dx,theta)*i+j*j+zrot(z-dz,i-dx,theta)*z)/(norm1*norm2);
                eta+=pd*eff(z)*eff(zrot(z-dz,i-dx,theta))*conj(gauss(i,j,z))*gauss(xrot(z-dz,i-dx,theta),j,zrot(z-dz,i-dx,theta));
            }
        }
    }
    return eta;
}

void couple_ang2(){
    Double_t xlo=-0.2;
    Double_t xhi=0.2;
    Double_t ylo=0;
    Double_t yhi=1;
    const Int_t n=200;
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,800,600);
    //TH1F *frame=gPad->DrawFrame(xlo,ylo,xhi,yhi);
    double wg=(2/pow(west(R0),2));
    //ウエストの違いがフェータルに効いてくるのでちゃんと調べないといけないのでは？
    double thetat=(lam/M_PI)*sqrt(wg);
    TGraph2D *graph2d = new TGraph2D;
    TGraph *graph1 = new TGraph;
    TGraph *graph2 = new TGraph;
    /*for(int i=0;i<n;i++){
        cout << i << endl;
        //double x=-0.01+0.0001*i;
        double dz=3*i-300;
        double y1=pow(abs(eta_2d(0,0,dz)),2);
        //double y2=exp(-2*x*x/pow(thetat,2));
        double y2=exp(-(dz*dz)/(2*58*58));
        graph1->SetPoint(i,dz,y1);
        graph2->SetPoint(i,dz,y2);
    }
    graph1->SetTitle("1z;z[mm];coupling");
    graph1->SetMarkerColor(kBlue);
    graph1->SetMarkerStyle(20);
    graph1->SetMarkerSize(0.7);
    graph2->SetLineColor(kRed);
    graph1->Draw("AP");
    graph2->Draw("CP");*/
    for(int i=0;i<n;i++){
        double theta = -1+0.01*i;
        for(int j=0;j<n;j++){
            double dx=0.2*j-20;
            double cp=pow(abs(eta_2d(dx,theta,0)),2);
            //cout << cp << endl;
            graph2d->SetPoint(1000*i+j,theta,dx,cp);
        }
        cout << i << endl;
    }
    axrange ax2d = {-1,1,-20,20,0,1,"Coupling;#theta[rad];x[mm];coupling"};
    graph2d->Draw("surf1");
}