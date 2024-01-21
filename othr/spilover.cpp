#include <iostream>
#include "../headers/fitter.h"
using namespace std;
double omegax = 157.085*sqrt(2);
double omegay = 164.195*sqrt(2);
double omegax0 = 2.82;//[mm]
double omegay0 = 2.71;//[mm]
double omegarec = 2.67;//[mm]
double dx = 0.05;//[mm]
double dy = 0.4;//[mm]
double dz = 10.03;//[mm]
double dtheta = 1.2*pow(10,-3);//[mrad]
double dphi = 0.5*pow(10,-3);//[mrad]
double c = 3*pow(10,11);//[mm/s]
//比で消える係数はガン無視で被積分部分のみを見る
double Epower(double x,double y){
    return exp(-2*x*x/(omegax*omegax)-2*y*y/(omegay*omegay));
}
double delta(double freq,double omegai){
    double lambda = c/(freq*pow(10,9));
    return sqrt((pow(omegarec*omegarec+omegai*omegai,2)+pow(lambda*dz/M_PI,2))/(omegarec*omegarec+omegai*omegai));
}
double R(double omega,double freq){
    double lamda = c/(freq*pow(10,9));
    double z0 = M_PI*omega*omega/lamda;
    return 1500*(1+(z0/1500)*(z0/1500));
}
double deltilt(double omegafoc,double freq){
    double lamda = c/(freq*pow(10,9));
    double z = (1/(omegarec*omegarec))+(1/(omegafoc*omegafoc));
    double z2 = (M_PI/lamda)*(M_PI/lamda)*pow(((1/R(omegarec,freq))-(1/R(omegafoc,freq))),2);
    return lamda*sqrt((z*z+z2)/z)/M_PI;
}
void spilover(){
    cout << 1-exp(-2*pow(dtheta/deltilt(omegax0,240),2))*exp(-2*pow(dphi/deltilt(omegay0,240),2)) << endl;
}