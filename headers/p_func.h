#include <iostream>
#include <math.h>
class p_func{
private:
    double c= 3.0*pow(10,8);
    double v0 = 220000.0;
    double vE = v0;
    double dNu=0.0000885;
    double vc=v0;
    double dnu=0.000076296;
public:
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
};