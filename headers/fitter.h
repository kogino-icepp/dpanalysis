#include <iostream>
#include <random>
#include "p_func.h"
#include <TMath.h>
#include <TROOT.h>
using namespace std;
#define rep(i,n) for(int i=0;i<n;i++)
#define prep(i,m,n) for(int i=m;i<n;i++)
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
        cout << "first" << endl;
        random_device rnd;
        mt19937 mt(rnd());
        uniform_int_distribution<> rand_bin(0,dbin-1);
        cout << "second" << endl;
        int p0,p1,p2;
        p0 = rand_bin(mt);
        cout << "third" << endl;
        while(p0==p1) p1 = rand_bin(mt);
        while(p2==p0 || p2==p1) p2 = rand_bin(mt);
        cout << "forth" << endl;
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
        cout << "third" << endl;
        f -> SetParameter(0,a);
        f -> SetParameter(1,b);
        f -> SetParameter(2,c);
        cout << a << " " << b << " " << c << endl;
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
    void rand_fit(TGraphErrors* graph,TF1* f,int ite,int fite,double fm,double fM,double &res){
        syoki_para(graph,f,0);
        rep(i,fite)graph -> Fit(f,"MQ","",fm,fM);
        double chi2 = f -> GetChisquare();
        double ndf = f -> GetNDF();
        double chimin = chi2/ndf;
        double chimin0 = chimin;
        vector<double> good_para(3);
        rep(i,3)good_para[i] = f -> GetParameter(i);
        
        rep(i,ite){
            syoki_rand1(graph,f);
            cout << i << endl;
            //syoki_rand2(f);
            rep(j,fite)graph -> Fit(f,"MQN","",fm,fM);
            chi2 = f -> GetChisquare();
            ndf = f -> GetNDF();
            chi2 /= ndf;
            //cout << chi2 << endl;
            if(chi2 < chimin){
                chimin = chi2;
                rep(k,3)good_para[k] = f -> GetParameter(k);
            }
        }
        rep(i,3) f -> SetParameter(i,good_para[i]);
        res = chimin;
    }
    void rand_fit_test(TGraphErrors* graph,TF1* f,int fite,double fm,double fM,double chi_test,int &num){
        syoki_para(graph,f,0);
        double chi2,ndf;
        rep(i,fite){
            graph -> Fit(f,"MQN","",fm,fM);
            chi2 = f -> GetChisquare();
            ndf = f -> GetNDF();
            
        }
        //cout << chi2/ndf << endl;
        chi2 = f -> GetChisquare();
        ndf = f -> GetNDF();
        double chimin = chi2/ndf;
        if(chi_test!=chimin){
            num++;
            cout << "Yes" << endl;
            cout << chi_test << " : " << chimin << endl;
        }
        
    }
    Double_t MyFunction(double x,double p0,double p1){
        return p0*TMath::Exp(-x/2)*pow(x,(p1/2)-1)/(pow(2,p1/2)*TMath::Gamma(p1/2));
    }
    //きちんと元の縮尺に戻してからヒストに入れてホワイトノイズの散らばりを見る
    void FillHist(TF1* f,TGraphErrors* graph,TH1D* hist,double yscale,double &dym){
        double x,y,dy,y1;
        double dymax = -100;
        rep(i,dbin){
            x = graph -> GetPointX(i);
            y = graph -> GetPointY(i);
            y1 = f -> Eval(x);
            dy = y1 - y;
            dy *= yscale;
            hist -> Fill(dy);
            if(abs(dy)>=dymax){
                //cout << dy << endl;
                dymax = abs(dy);
            }
        }
        if(dym<dymax)dym = dymax;
        //cout << "dymax = " << dymax << endl;

    }
    void make_scale(TGraphErrors* graph,TGraph* mgraph,int sbin,double &yscale){
        //走査範囲のレンジ調査
        double xmin,xmax,ymin,ymax;
        ymin = 10000;
        ymax = -200;
        double x1 = mgraph -> GetPointX(sbin);
        double x2 = mgraph -> GetPointX(sbin+dbin-1);
        xmin = min(x1,x2);
        xmax = max(x1,x2);
        double x,y;
        prep(bin,sbin,sbin+dbin){
            y = mgraph -> GetPointY(bin);
            if(y<ymin)ymin = y;
            if(y>ymax)ymax = y;
        }
        //レンジに合わせてセットポイントを変える
        //x -> (x-x0)/xrange, y ->  (y-ymin)/yrange ??
        
        double xrange = xmax-xmin;
        double yrange = ymax-ymin;
        prep(bin,sbin,sbin+dbin){
            x = mgraph -> GetPointX(bin);
            y = mgraph -> GetPointY(bin);
            x = (x-xmin)/xrange;
            y = (y-ymin)/yrange;
            graph -> SetPoint(bin-sbin,x,y);
            //cout << x << " " << y << endl;
            graph -> SetPointError(bin-sbin,0,0.1/yrange);
        }
        yscale = yrange;
    }
    vector<double> fparameters(TF1* f,int pnum){
        vector<double> para;
        rep(i,pnum){
            double p;
            p = f -> GetParameter(i);
            para.push_back(p);
        }
        return para;
    }
};