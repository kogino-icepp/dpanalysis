#include <iostream>
#include <random>
#include "setting.h"
#include <TMath.h>
#include <TROOT.h>
using namespace std;
#define rep(i,n) for(int i=0;i<n;i++)
const int dbin = 30;
class Fitter{
    public:
    void syoki_para(TGraphErrors* graph,TF1* f,int bin){
        //y=a(x-x0)(x-x2)+mx+nが(x1,y1)を通る
        //取得するパラメータはy=a(x-b)^2+cの方
        double x0 = graph -> GetPointX(bin);
        double y0 = graph -> GetPointY(bin);
        double x1 = graph -> GetPointX(bin+dbin/2);
        double y1 = graph -> GetPointY(bin+dbin/2);
        double x2 = graph -> GetPointX(bin+dbin-1);
        double y2 = graph -> GetPointY(bin+dbin-1);
        /*cout << x0 << " " << x1 << " " << x2 << endl;
        cout << y0 << " " << y1 << " " << y2 << endl;*/
        
        double m = (y2-y0)/(x2-x0);
        double n = (y0*x2-y2*x0)/(x2-x0);
        double a = (y1-m*x1-n)/((x1-x0)*(x1-x2));
        double b = (x0+x2-(m/a))/2;
        double c = a*x0*x2+n-a*b*b;
        /*cout << "a = " << a << endl;
        cout << "b = " << b << endl;
        cout << "c = " << c << endl;*/
        f -> SetParameter(0,a);
        f -> SetParameter(1,b);
        f -> SetParameter(2,c);
    }
    void syoki_para_d(TGraphErrors* graph,TF1* f,int bin){
        vector<double> ans(3);
        //y=a(x-x0)(x-x2)+mx+nが(x1,y1)を通る
        //取得するパラメータはy=a+bx+cx^の方
        double x0 = graph -> GetPointX(bin);
        double y0 = graph -> GetPointY(bin);
        double x1 = graph -> GetPointX(bin+dbin/2);
        double y1 = graph -> GetPointY(bin+dbin/2);
        double x2 = graph -> GetPointX(bin+dbin-1);
        double y2 = graph -> GetPointY(bin+dbin-1);
        
        double m = (y2-y0)/(x2-x0);
        double n = (y0*x2-y2*x0)/(x2-x0);
        double a = (y1-m*x1-n)/((x1-x0)*(x1-x2));
        double b = m-a*(x0+x2);
        double c = a*x0*x2+n;
        f -> SetParameter(0,c);
        f -> SetParameter(1,b);
        f -> SetParameter(2,a);
    }
    //初期値を、ランダムな3点のデータ点を取ることで与える
    void syoki_rand1(TGraphErrors* graph,TF1* f){
        random_device rnd;
        mt19937 mt(rnd());
        uniform_int_distribution<> rand_bin(0,dbin-1);
        int p0,p1,p2;
        p0 = rand_bin(mt);
        while(p0==p1) p1 = rand_bin(mt);

        while(p2==p0 || p2==p1) p2 = rand_bin(mt);

        double x0 = graph -> GetPointX(p0);
        double y0 = graph -> GetPointY(p0);
        double x1 = graph -> GetPointX(p1);
        double y1 = graph -> GetPointY(p1);
        double x2 = graph -> GetPointX(p2);
        double y2 = graph -> GetPointY(p2);
        
        double m = (y2-y0)/(x2-x0);
        double n = (y0*x2-y2*x0)/(x2-x0);
        double a = (y1-m*x1-n)/((x1-x0)*(x1-x2));
        double b = (x0+x2-(m/a))/2;
        double c = a*x0*x2+n-a*b*b;
        
        f -> SetParameter(0,a);
        f -> SetParameter(1,b);
        f -> SetParameter(2,c);
    }
    /*
    初期値をランダムな3つのパラメータを振ることで最適化を目指す
    課題：どの程度探す？(乱数で回数調整or一通り走査)
    */
    void syoki_rand2(TF1*f){
        random_device rnd;
        mt19937 mt(rnd());
        uniform_real_distribution<> rand0(-1,1);
        uniform_real_distribution<> rand1(-2,2);
        uniform_real_distribution<> rand2(0,1);
        f -> SetParameter(0,rand0(mt));
        f -> SetParameter(1,rand1(mt));
        f -> SetParameter(2,rand2(mt));
    }
    void rand_fit(TGraphErrors*graph,TF1*f,int ite,int fite,double fm,double fM){
        double chimin = 100.0*27;
        vector<double> good_para(3);
        rep(i,ite){
            syoki_rand1(graph,f);
            //syoki_rand2(f);
            rep(j,fite)graph -> Fit(f,"MQ","",fm,fM);
            double chi2 = f -> GetChisquare();
            double ndf = f -> GetNDF();
            chi2 /= ndf;
            cout << chi2 << endl;
            if(chi2 < chimin){
                chimin = chi2;
                rep(k,3)good_para[k] = f -> GetParameter(k);
            }
        }
        rep(i,3) f -> SetParameter(i,good_para[i]);
    }
    Double_t MyFunction(double x,double p0,double p1){
        return p0*TMath::Exp(-x/2)*pow(x,(p1/2)-1)/(pow(2,p1/2)*TMath::Gamma(p1/2));
    }
    double ChiValue(TGraphErrors* graph,TF1* f,int boffset){
        int npoints = graph -> GetN();
        double chi2 = 0.0;
        double x,y,ex,ey;
        rep(i,npoints){
            graph -> GetPoint(i+boffset,x,y);
            ex = graph -> GetErrorX(i+boffset);
            ey = graph -> GetErrorY(i+boffset);

            double expectedY = f -> Eval(x);
            double residual = (y-expectedY)/ey;

            chi2 += residual*residual;
        }
        return chi2;
    }
    /*void getmnstat(TMinuit*minu){
        Double_t fmin,fedm,errdef;
        Int_t npari,nparx,istat;
        minu -> mnstat(fmin,fedm,errdef,npari,nparx,istat);
        cout << "fmin = " << fmin << endl;

    }*/
    //要件定義:パラメータ空間において二次元グラフを作成する
    void pfield(TGraphErrors*graph,TGraph2D* graph2,TF1* f,double offset){
        int bin = 0;
        for(double i=-0.1;i<1.1;i+=0.1){
            for(double j=0;j<100;j+=1){
                f -> SetParameter(0,i);
                f -> SetParameter(1,j);
                f -> SetParameter(2,offset);
                double chi2 = ChiValue(graph,f,0);
                int ndf = graph -> GetN();
                ndf -= 3;
                cout << i << " " << j << " " << chi2/ndf << endl;
                graph2 -> SetPoint(bin,i,j,chi2/ndf);
                bin++;
            }
        }
    }
};