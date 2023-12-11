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
double errdBmtomW(double dbm,double ddbm=0.01){
    //dBmの読み取り誤差が0.01dBmであると仮定、ダメなら適切に手で与える
    return ddbm*(dBmtomW(dbm)/10)*log(10);
}
double mWtodBm(double mw){
    return (10*log(mw))/log(10);
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
double thetadbm[12] = {-7.86,-7.97,-5.65,-5.71,-5.55,-5.48,-5.49,-5.88,-6.7,-5.48,-5.49,-5.48};
double phidbm[12] = {-8.19,-7.67,-8.05,-7.22,-5.85,-5.60,-7.97,-6.30,-5.78,-5.84,-6.8,-6.98};
double phi[12] = {0,1.0/6,1.0/3,2.0/3,5.0/6,1.0,-1.0/6,-1.0/3,-0.5,-2.0/3,-5.0/6,-1};
double theta[12] = {0,0.5,1,1.5,2,2.5,3,-0.5,-1,-1.5,-2,-2.5};

double xW1[7] = {1.097,0.993,0.866,0.830,0.866,0.93,0.89};
double xW2[8] = {66.0,65.0,65.5,65.3,66.2,67.3,68.4,70.1};
double xpos1[7] = {0,2,4,6,8,10,12};
double xpos2[8] = {6.4,6.0,5.5,5.0,4.5,7,7.5,8};
double yW1[11] = {0.824,0.728,0.705,0.723,0.695,0.687,0.692,0.671,0.669,0.707,0.819};
double yW2[9] = {0.655,0.657,0.652,0.658,0.66,0.671,0.683,0.72,0.757};
double ypos1[11] = {21,19,17,15,13,11,9,7,5,3,1};
double ypos2[9] = {7,7.5,8,8.5,9,6.5,6,5.5,5};

double xpower[20] = {77.8,79.4,82.4,87.5,94.8,104.5,111.4,116.5,119.0,120.1,118.1,112.1,103.8,95.1,88.0,83.4,80.9,78.9,78.1,77.3};
double xpos[20] = {620,590,560,530,500,470,440,410,380,350,320,290,260,230,200,170,140,110,80,50};
double ypos[21] = {800,770,740,710,680,650,620,590,560,530,500,470,440,410,380,350,320,290,260,230,200};
double ypower[21] = {77.4,78.3,79.4,81.2,83.9,88.4,95.3,104.6,113.8,121.9,125.6,124.8,120.7,111.6,103.4,93.1,86.3,81.9,79.4,78.3,77.9};
double temp[3] = {-71.4,-197.2,23.5};
double xx[7] = {800,700,600,500,400,300,200};

double xbpower[20] = {69.6,70.2,71.0,72.5,75.3,78.6,84.6,92.1,99.4,106.3,109.9,109.4,105.9,99.1,90.9,82.6,76.7,72.6,71.0,70.1};
double xbpos[20] = {800,770,740,710,680,650,620,590,560,530,500,470,440,410,380,350,320,290,260,230};
double ybpower[20] = {70.4,71.3,72.0,73.7,76.3,80.4,93.7,101.1,107.9,110.1,108.7,107.0,101.1,99.8,87.3,80.8,76.4,73.9,72.2};
double ybpos[20] = {50,80,110,140,170,200,230,260,290,320,350,380,410,440,470,500,530,560,590};
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
    TGraphErrors* gthetaaline = new TGraphErrors;
    TGraphErrors* gphialine = new TGraphErrors;
    
    TGraphErrors* gmirrorx = new TGraphErrors;
    TGraphErrors* gmirrory = new TGraphErrors;
    rep(i,11){
        gthetaaline ->SetPoint(i,theta[i],dBmtomW(thetadbm[i]));
        gphialine -> SetPoint(i,phi[i],dBmtomW(phidbm[i]));
        gphialine -> SetPointError(i,0.172,errdBmtomW(phidbm[i]));
        gxaline -> SetPoint(i,xpos1[i],xW1[i]);
        gxaline -> SetPointError(i,0.01,0.05);
        gyaline -> SetPoint(i,ypos1[i],yW1[i]);
        gyaline -> SetPointError(i,0.01,0.01);
    }
    rep(i,20){
        gmirrorx -> SetPoint(i,xbpos[i],xbpower[i]);
        double xerr = xbpower[i]*log(10)/100;
        gmirrorx -> SetPointError(i,0,xerr);
        gmirrory -> SetPoint(i,ybpos[i],ybpower[i]);
        cout << ybpos[i] << " " << ybpower[i] << endl;
        double yerr = ybpower[i]*log(10)/100;
        gmirrory -> SetPointError(i,0,yerr);
    }
    axrange ax = {0,22,0,1,0,1,";x[mm];Power[mW]"};
    TF1* fgaus = new TF1("fgaus","[0]*exp(-(x-[1])*(x-[1])/([2]*[2]))+[3]",0,22);
    fgaus -> SetParameter(3,60);
    fgaus -> SetParameter(2,100);
    fgaus -> SetParameter(0,1);
    fgaus -> SetParameter(1,500);
    
    st.GraphErrors(gyaline,ax);
    gyaline -> Draw("AP");
    //gyaline -> Fit(fgaus,"E");
    axrange axg = {0,1000,0,200,0,1,";x[mm];Power[#mu W]"};
    st.GraphErrors(gmirrorx,axg);
    st.GraphErrors(gmirrory,axg);
    gmirrory -> SetLineColor(kBlue);
    gmirrory -> Draw("AP");
    gmirrory -> Fit(fgaus);
}