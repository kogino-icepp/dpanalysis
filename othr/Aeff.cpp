#include <iostream>
#include <cmath>
#include "../headers/setting.h"
using namespace std;
const double R = 1500;//球面の曲率半径のこと
const double c = 3*pow(10,11);//光速[mm/s]
const double candfreq[3] = {218,240,264};
const double sigmax = 2.82;
const double sigmay = 2.71;//y方向のウエストサイズ[mm]
const double dtheta = M_PI/20000;//θの刻み値[rad]
const double dphi = M_PI/1000;//φの刻み値[rad]
const double theta_max = asin(344.75/1500);
const double dR = R*dtheta/(2*M_PI);
const double xkei = 131;
const double ykei = 137;
const double dxkei = 6;
const double dykei = 6;
const double dal = 0.001;
#define rep(i,n) for(int i=0;i<n;i++)
#define prep(i,n,m) for(int i=n;i<m;i++)
//ミラーの上で球面積分を行うプログラム
double step(double x){
    if(x>=0)return 1.0;
    else return 0.0;
}
double gauss_hosei(double x,double y,double height,double sigx,double sigy,double r_hatched){
    double r = sqrt(x*x+y*y);
    return height*exp(-(x/sigx)*(x/sigx)-(y/sigy)*(y/sigy));
}

void Aeff(){
    Setting st;
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,700);
    c1 -> SetMargin(0.15,0.15,0.2,0.1);
    double aeff = 0;//有効面積
    double reA = 0;
    double freq = candfreq[0]*pow(10,9);
    double wave_len = c/freq;
    double alpha_x = wave_len/(M_PI*sigmax);
    double alpha_y = wave_len/(M_PI*sigmay);
    double x_sigma = (alpha_x+dal)*(144.845+1315);
    double y_sigma = (alpha_y+dal)*(144.845+1315);
    TGraph2D* mapgraph = new TGraph2D;
    int mbin = 0;
    cout << gauss_hosei(0,0,1,200,200,0) << endl;;
    for(double theta=theta_max;theta>=0;theta-=dtheta){
        double r = 1500*sin(theta);
        double dz = R*(cos(theta)-cos(theta_max));
        double dr = dR*cos(theta);
        x_sigma += alpha_x*dz;
        y_sigma += alpha_y*dz;
        
        for(double phi=0;phi<=2*M_PI;phi+=dphi){
            double x_mesh = r*cos(phi);
            double y_mesh = r*sin(phi);
            reA += R*R*sin(theta)*dphi*dtheta;
            aeff += R*R*sin(theta)*dphi*dtheta*gauss_hosei(x_mesh,y_mesh,1,x_sigma,y_sigma,0);
            mapgraph -> SetPoint(mbin,x_mesh,y_mesh,100*gauss_hosei(x_mesh,y_mesh,1,x_sigma,y_sigma,0));
            mbin++;
            //cout << x_mesh << " " << y_mesh << " " << 100*gauss_hosei(x_mesh,y_mesh,1,x_sigma,y_sigma,0) << endl;
            //cout << gauss_hosei(x_mesh,y_mesh,1,csigmax,csigmay,1);
        }
        //cout << endl;
    }
    //cout << M_PI << endl;
    cout << reA << " : " << aeff << " -> " << aeff/reA << endl;
    axrange ax2d = {-350,350,-350,350,0,100,";X[mm];Y[mm];Power"};
    st.Graph2D(mapgraph,ax2d);
    mapgraph -> Draw("surf1");
    return;
}