#include <iostream>
#include "../headers/fitter.h"
using namespace std;
# define PI 3.14159265359
double c = 3*pow(10,11);
double w0 = 2.67;//[mm]

double etaz(){
    return 0;
}
double z0(double freq){
    double lambda = c/freq;
    return PI*w0*w0/lambda;
}
double R(double z,double freq){
    return z*(1+(z0(freq)/z)*(z0(freq)/z));
}
double eta2(double R){
    return 4/(((w0/(0.00089*R))+((0.00089*R)/w0))*((w0/(0.00089*R))+((0.00089*R)/w0)));
}
double dBmtomW(double dbm){
    return pow(10,(dbm)/10);
}
double xdbm[7] = {-11.55,-11.25,-10.66,-10.12,-10.67,-11.42,-11.63};
double xpower[20] = {77.8,79.4,82.4,87.5,94.8,104.5,111.4,116.5,119.0,120.1,118.1,112.1,103.8,95.1,88.0,83.4,80.9,78.9,78.1,77.3};
double xpos[20] = {620,590,560,530,500,470,440,410,380,350,320,290,260,230,200,170,140,110,80,50};
double ypos[21] = {800,770,740,710,680,650,620,590,560,530,500,470,440,410,380,350,320,290,260,230,200};
double ypower[21] = {77.4,78.3,79.4,81.2,83.9,88.4,95.3,104.6,113.8,121.9,125.6,124.8,120.7,111.6,103.4,93.1,86.3,81.9,79.4,78.3,77.9};
double temp[3] = {-71.4,-197.2,23.5};
double xx[7] = {800,700,600,500,400,300,200};

void coupling(){
    Setting st;
    st.dot_size = 0.8;
    st.markerstyle = 20;
    st.color = kBlue;
    st.lcolor = kBlue;

    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    
    TGraphErrors* gxaline = new TGraphErrors;
    TGraphErrors* gyaline = new TGraphErrors;
    TGraphErrors* glinear = new TGraphErrors;
    glinear -> SetPoint(0,273.15+temp[0],1.45211);
    glinear -> SetPoint(1,273.15+temp[1],1.02329);
    glinear -> SetPoint(2,273.15+temp[2],1.78258);
    //errorもつける
    glinear -> SetPointError(0,1,1.45211*log(10)*0.1/10);
    glinear -> SetPointError(1,1,1.02329*log(10)*0.1/10);
    glinear -> SetPointError(2,1,1.78258*log(10)*0.1/10);
    axrange axl = {0,300,0,2,0,1,";Temperature[K];Power[uW]"};
    st.GraphErrors(glinear,axl);
    glinear -> Draw("AP");
    TF1* lf = new TF1("lf","[0]+[1]*x");
    glinear -> Fit(lf);
    // rep(bin,21){
    //     //gxaline -> SetPoint(bin,xpos[bin],xpower[bin]);
    //     //gxaline -> SetPointError(bin,1,xpower[bin]*log(10)/20);
    //     gyaline -> SetPoint(bin,ypos[bin],ypower[bin]);
    //     gyaline -> SetPointError(bin,1,ypower[bin]*log(10)/20);
    // }
    // axrange ax = {0,900,0,200,0,1,";pos[mm];Power[uW]"};
    // st.GraphErrors(gxaline,ax);
    // st.GraphErrors(gyaline,ax);
    // //gxaline -> Draw("AP");
    // gyaline -> Draw("AP");
    // TF1* fgaus = new TF1("fgaus","[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))+[3]");
    // fgaus -> SetParameter(3,0.05);
    // fgaus -> SetParameter(0,0.01);
    // fgaus -> SetParameter(1,500);
    // fgaus -> SetParameter(2,100);
    
    // gyaline -> Fit(fgaus);
    //各パラメータを変化させたときの理想的なカップリングロスを示してみる(二次元グラフとかがいいのかな)
}