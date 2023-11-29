#include <iostream>
#include "../headers/fitter.h"
using namespace std;
# define PI 3.14159265359
double c = 3*pow(10,11);//[mm/s]
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
//dzによるカップリング、ウエストサイズが違う場合を考えるときは補正かける
double etaz(double dz,double f){
    double lambda = c/(f*pow(10,9));//[mm]になるように調節(できているはず)
    return 4/(4+(lambda*dz/(PI*w0*w0))*(lambda*dz/(PI*w0*w0)));
}
//dzだけでなくdrも考慮した時のcoupling,上と同じくウエストサイズまでちゃんと考えるときは補正かける
double etar(double dz,double dr,double f){
    double lambda = c/(f*pow(10,9));
    double delta = sqrt(((4*pow(w0,4))+(lambda*dz/PI)*(lambda*dz/PI))/(2*w0*w0));
    return etaz(dz,f)*exp(-2*(dr/delta)*(dr/delta));
}
//角度がずれた時のカップリングの変化、
double etaarg(double dz,double dphi,double f){
    double lambda = c/(f*pow(10,9));
    double R = 1500;
    double deno = ((1/(w0*w0))+(1/(w0*w0)));
    double deltil = sqrt(((deno*deno)+(PI/(lambda*R))*(PI/(lambda*R)))/deno);
    return etaz(dz,f)*exp(-2*(dphi*dphi/(deltil*deltil)));
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
    TH2D* euceta = new TH2D("euceta","#eta (#Deltaz,#Delta r);#Delta z[mm];#Delta r[mm]",20,0,0.1,20,0,0.1);
    TH2D* argeta = new TH2D("argeta","#eta (#Deltaz,#Delta#theta);#Delta z[mm];#Delta #theta[rad]",20,0,0.025,20,0,0.025);

    rep(i,21){
        rep(j,21){
            euceta -> SetBinContent(i,j,etar(i/200.0,j/200.0,240));
            argeta -> SetBinContent(i,j,etaarg(0.00125*i,0.00125*j,240));
        }

    }
    st.Hist(euceta);
    st.Hist(argeta);
    euceta -> SetStats(0);
    euceta -> SetMinimum(0.997);
    euceta -> SetMaximum(1);
    euceta -> Draw("colz");
    argeta -> SetStats(0);
    //euceta -> SetMinimum(0.997);
    //euceta -> SetMaximum(1);
    argeta -> Draw("colz");
    /*st.GraphErrors(glinear,axl);
    glinear -> Draw("AP");
    TF1* f1 = new TF1("f1","[0]+[1]*x",0,300);
    glinear -> Fit(f1);
    //各軸からのずれによってどの程度カップリングロスがあるのかをTH2Dで表示してみる*/
    
}