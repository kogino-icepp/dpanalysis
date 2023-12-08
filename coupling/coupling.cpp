#include <iostream>
#include "../headers/fitter.h"
using namespace std;
double cv=3.0*pow(10,11); 
int freq=240; 
double lam=cv/(freq*pow(10,9)); 
double w0=2.67;
double z0=M_PI*w0*w0/lam;
double k=2*M_PI/lam;
double Rm=1500;
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
complex<double> gauss(double x,double y,double z){
    complex<double> kata(-(x*x+y*y)/(west(z)*west(z)),-k*z-M_PI*(x*x+y*y)/(lam*R(z)));
    complex<double> g=exp(kata);
    return g;
}
complex<double> etax(double dx,double dz){
    complex<double> eta(0);
    for(int x=-350;x<351;x++){
        for(int y=-350;y<351;y++){
            if(x*x+y+y<=350*350){
                double z=sqrt(1500*1500-x*x-y*y);
                double ab1=sqrt(x*x+y*y+z*z);
                double ab2=sqrt((x-dx)*(x-dx)+y*y+(z-dz)*(z-dz));
                eta+=eff(z)*eff(z-dz)*conj(gauss(x-dx,y,z-dz))*gauss(x,y,z)*((x/ab1)*((x-dx)/ab2)+(y*y)/(ab1*ab2)+(z*(z-dz))/(ab1*ab2));
            }
        }
    }
    return eta;
}
//ウエスト位置でカップリングを計算するように路線変更(ってどうやるの？)
//ミラーでビームが跳ね返った後
double draf(double r,double dr){
    double rtn=(-2/Rm)*r+dr;
    return rtn;
}
complex<double> etaw(double dx,double dz){
    complex<double> etaw(0);
    for(int x=0;x<100;x++){
        for(int y=0;y<100;y++){
            
        }
    }
}
void coupling()
{
    // キャンバスの準備
    int n=200;
    TGraph2D* graph = new TGraph2D;
    //データーの準備
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
  // グラフを作る
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            double x=0;
            double y=0;
            double z=0;
            x=0.1*i-10;
            y=0.1*j-10;
            z=abs(etax(x,y));
            cout <<x <<" " <<" " <<y<< " "<< z << endl;
            graph->SetPoint(1000*i+j,x,y,z);
        }
    }
    
    graph->SetTitle("Coupling;x[mm];z[mm];coupling");
    gStyle->SetPalette(1);
    axrange ax2d = {-10,10,-10,10,0,1,"Coupling;x[mm];z[mm];coupling"};
    st.Graph2D(graph,ax2d);
    graph->Draw( "surf1" );        // グラフを書く 
    
}