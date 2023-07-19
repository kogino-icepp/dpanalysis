#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <format>
#include "setting.h"
#include <TROOT.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TStyle.h>
using namespace std;
#define rep(i,n) for(int i=0;i<n;i++)
#define prep(i,m,n) for(int i=m;i<n;i++)
#define fore(i,v) for(auto& i:v)

//特定のデータだけ引っ張ってくる用
vector<int> lsbo = {1,3,5,13,15,17};
vector<int> lsbi = {2,4,6,14,16,18};
vector<int> usbi = {7,9,11,19,21,23};
vector<int> usbo = {8,10,12,20,22,24};
vector<int> lcolor = {4,2,3,6,7,5};
vector<int> sbin = {0,512,-512,-1024};
//共通の物理量
const double dcsigma = 0.000487;
const double dhsigma = 0.000457;
const double ddcsigma = 0.000508;
const double ddhsigma = 0.000496;
const double ddmsigma = 0.0004986;
const double dclim = 6.2;
const double dhlim = 6.2;
const double ddclim = 5;
const double ddhlim = 5;
const double ddmlim = 10;
const int nbin=32767;
const int ssb = 2621;//探索すべきビンの最初
const int sfb = 28835;//探索すべきビンの最後
const int sb = 2585;//切り出してくるビンの最初
const int fb = 28865;//切り出してくるビンの最後
const double gcsigma = 0.00072;
const double ghsigma = 0.00112;
const double calpha = 6.5;
const int total_num=(fb-sb)*2*4;
const double c= 3.0*pow(10,8);
const double v0=220000.0;
const double vE=v0;
const double vc=v0;
const double dnu=0.000076296;//ビン幅
const double dNu=0.0000885;//周波数分解能
const int dbin = 30;//フィッティングするビンの幅
const double kb=1.38*pow(10,-23);
const double df=88.5*pow(10,3);
const double Tc=76;
const double Th=297;

bool toge_hantei(vector<bool>ctoge, vector<bool>htoge, int bin, queue<int>& que){
    prep(i,bin,bin+dbin){
        if(ctoge[i] || htoge[i]){
            que.push(i);
            return true;
        }
    }
    return false;
}
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
int both = 0;
int tsafe = 0;
int tout = 0;
void toge_scan(bool (&hantei)[nbin],double input[nbin],double sigma,double limit){
    double binput = input[sb-1];
    double bbinput = input[sb-2];
    double ddinput;
    for(int bin=sb;bin<fb;bin++){
        ddinput = input[bin] + bbinput - 2*binput;
        ddinput /= input[bin]*sigma;
        
        bbinput = binput;
        binput = input[bin];
        if(abs(ddinput) > limit)hantei[bin] = true;
    }
}
void toge_value(double input[nbin],double (&output)[nbin],double sigma){
    double binput = input[sb-1];
    double bbinput = input[sb-2];
    double ddinput;
    for(int bin=sb;bin<fb;bin++){
        ddinput = input[bin] + bbinput - 2*binput;
        ddinput /= input[bin]*sigma;
        output[bin] = ddinput;
        bbinput = binput;
        binput = input[bin];
    }
}
vector<double> caliblate(double hot[nbin],double cold[nbin], double mirror[nbin]){
    long double gain,psys;
    vector<double> prec(nbin);
    rep(bin,nbin){
        gain=(hot[bin]-cold[bin])/(2*kb*(Th-Tc)*df);
        psys=(cold[bin]/gain)-2*kb*Tc*df;
        prec[bin]=((mirror[bin]/gain)-psys)/(2*kb*df);
        
    }
    return prec;
}
/*
要件定義
1. スケール関数-f
2. キャリブったグラフのレンジ分fを戻す(xはわかりそう、yはログ残しとかなきゃ)
3. 戻した後のデータを合流させる
*/
void rescale_graph(TGraphErrors* sgraph,TGraphErrors* resgraph,TF1* f,int& pos,double yrange,int bin,double freq[nbin],TH1D* hist){
    double p0 = f -> GetParameter(0);
    double p1 = f -> GetParameter(1);
    double p2 = f -> GetParameter(2);
    double x,y;
    rep(i,dbin){
        x = sgraph -> GetPointX(i);
        y = sgraph -> GetPointY(i);
        y -= p0*(x-p1)*(x-p1)+p2;
        y *= yrange;
        resgraph -> SetPoint(pos,freq[bin+i],y);
        pos++;
        hist -> Fill(y);
    }
}
void make_scale(TGraphErrors* graph,TGraph* mgraph,int sbin,double &yoffset,double &yran){
    //走査範囲のレンジ調査
    double xmin,xmax,ymin,ymax;
    ymin = 10000;
    ymax = -200;
    double x1 = mgraph -> GetPointX(sbin);
    double x2 = mgraph -> GetPointX(sbin+dbin);
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
    yran = yrange;
    yoffset = ymin;
    prep(bin,sbin,sbin+dbin){
        x = mgraph -> GetPointX(bin);
        y = mgraph -> GetPointY(bin);
        x = (x-xmin)/xrange;
        y = (y-ymin)/yrange;
        graph -> SetPoint(bin-sbin,x,y);
        //cout << x << " " << y << endl;
        graph -> SetPointError(bin-sbin,0,0.1/yrange);
    }
}
//パラメータ数を可変にしたときのフィッティング関数 できればフィット後のステータスも知りたい

/*
要件定義
1. coldとhotでアウトなものは無条件に抜く
2. mirrorに関しては4回分データを走査し、
*/
void peak_search2(){
    Setting setting;
    setting.dot_size = 0.8;
    setting.markerstyle = 20;
    setting.color = kGreen;
    //お絵描きの設定など
    Double_t xlo=215;
    Double_t xhi=264;
    Double_t ylo=0;
    Double_t yhi=100;
 
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    TH1F *frame=gPad->DrawFrame(xlo,ylo,xhi,yhi);
    filesystem::path path=filesystem::current_path();
    string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
    string savedir = "/Users/oginokyousuke/data/baseline_change/";
    
    filesystem::current_path(dir);

    TGraph* g_ddcm = new TGraph;
    TGraph* g_ddhm = new TGraph;
    int dcnum = 0;
    int dhnum = 0;
    
    queue<double> ratio;
    double rsum = 0;
    for(int i=1;i<2;i++){
        double good_num = 0;
        double bad_num = 0;
        cout << i << endl;
        //int lonum = 2;
        //if(i!=lsbo[lonum] && i!=lsbi[lonum] && i!=usbi[lonum] && i!=usbo[lonum])continue;
        double ifmin = 213.9+2*i;
        double ifmax = 216.1+2*i;
        string cdir=dir+"band"+to_string(i);
        filesystem::current_path(cdir);
        queue<pair<int,double>> que1;
        queue<pair<int,double>> que2;
        double prec_av[34303][8];
        double freq_av[34303];
        rep(bin,34303){
            freq_av[bin] = -1;
            rep(k,8)prec_av[bin][k] = -10000;//ありえない数字を入れておけばとりあえずオッケー
        }
        //事前に配列として持っておけば良いのでは？
        bool c_toge1[4][nbin],c_toge2[4][nbin],h_toge1[4][nbin],h_toge2[4][nbin],m_toge1[4][nbin],m_toge2[4][nbin];

        //フィットの改善を表すグラフ
        TGraph* gchange1 = new TGraph;
        TGraph* gchange2 = new TGraph;
        int gcbin1 = 0;
        int gcbin2 = 0;
        TH1D* chi_hist = new TH1D("chi_hist","chi_hist;chi2/Ndf;Count",100,0,10);

        for(int j=0;j<1;j++){
            double Freq1[nbin],cold1[nbin],hot1[nbin],mirror1[nbin],Freq2[nbin],cold2[nbin],hot2[nbin],mirror2[nbin];
            //ビンのシフトをここでいじる
            for(int data=0;data<576;data++){
                string file_name="test"+to_string(data)+".root";
                if(filesystem::exists(file_name)){
                    TFile* file = new TFile(file_name.c_str());
                    TTree* tree = (TTree*)file->Get("tree");
                    int time,xfft,shift,body,num;
                    long double freq,power;
                    tree->SetBranchAddress("time",&time);
                    tree->SetBranchAddress("XFFT",&xfft);
                    tree->SetBranchAddress("shift",&shift);
                    tree->SetBranchAddress("freq",&freq);
                    tree->SetBranchAddress("power",&power);
                    tree->SetBranchAddress("num",&num);
                    tree->SetBranchAddress("body",&body);
                    tree->GetEntry(0);
                    if(shift==sbin.at(j)){
                        //cout << file_name << endl;
                        if(num==1 && body==0){
                            for(int bin=0;bin<nbin;bin++){
                                tree -> GetEntry(bin);
                                Freq1[bin]=freq;
                                cold1[bin]=power;
                            }
                        }
                        else if(num==1 && body==1){
                            for(int bin=0;bin<nbin;bin++){
                                tree->GetEntry(bin);
                                hot1[bin]=power;
                            }
                        }
                        else if(num==1 && body==2){
                            for(int bin=0;bin<nbin;bin++){
                                tree->GetEntry(bin);
                                mirror1[bin]=power;
                            }
                        }
                        else if(num==0 && body==0){
                            for(int bin=0;bin<nbin;bin++){
                                tree->GetEntry(bin);
                                cold2[bin]=power;
                                Freq2[bin]=freq;
                            }
                        }
                        else if(num==0 && body==1){
                            for(int bin=0;bin<nbin;bin++){
                                tree->GetEntry(bin);
                                hot2[bin]=power;
                            }
                        }
                        else if(num==0 && body==2){
                            for(int bin=0;bin<nbin;bin++){
                                tree->GetEntry(bin);
                                mirror2[bin]=power;
                            }
                        }
                    }
                }
            }
            
            vector<double> prec1 = caliblate(hot1,cold1,mirror1);
            vector<double> prec2 = caliblate(hot2,cold2,mirror2);
            TGraphErrors* pgraph1 = new TGraphErrors;
            TGraphErrors* pgraph2 = new TGraphErrors;

            prep(bin,sb,fb){
                //cout << bin << " " << Freq1[bin] << " " << prec1[bin] << endl;
                pgraph1 -> SetPoint(bin-sb,Freq1[bin],prec1[bin]);
                pgraph1 -> SetPointError(bin-sb,0,0.1);
                pgraph2 -> SetPoint(bin-sb,Freq2[bin],prec2[bin]);
                pgraph2 -> SetPointError(bin-sb,0,0.1);
            }
            double ddc1[nbin],ddc2[nbin],ddh1[nbin],ddh2[nbin],ddm1[nbin],ddm2[nbin];
            toge_scan(c_toge1[j],cold1,ddcsigma,ddclim);
            toge_scan(c_toge2[j],cold2,ddcsigma,ddclim);
            toge_scan(h_toge1[j],hot1,ddhsigma,ddhlim);
            toge_scan(h_toge2[j],hot2,ddhsigma,ddhlim);
            
            setting.GraphErrors(pgraph1,ifmin,ifmax,0,100);
            pgraph1 -> SetTitle("Prec;Freq[GHz];Prec[K]");
            pgraph1 -> Draw("AP");
            
            double fm,fM;
            
            //キャリブった後のデータをグラフに貼る+ヒストに詰める→ヒストの結果を元にerrorをグラフに追加で載せる
            //一回目のデータの処理
            TGraphErrors* respgraph = new TGraphErrors;
            int pos = 0;
            double yoffset,yran;
            double yrange;
            TH1D* white_hist = new TH1D("white_hist","white_hist;white_noise[K];Count",100,-1,1);
            for(int bin=sb;bin<fb;bin+=dbin){
                //スケール
                TGraphErrors* spgraph = new TGraphErrors;
                //void make_scale(TGraphErrors* graph,TGraph* mgraph,int sbin,double &yoffset,double &yran)
                make_scale(spgraph,pgraph1,bin-sb,yoffset,yran);
                //棘のパートはパス
                bool hantei = false;
                prep(k,bin,bin+dbin){
                    if(c_toge1[j][k] || h_toge1[j][k]){
                        hantei = true;
                        break;
                    }
                }
                if(hantei)continue;
                
                TF1* f2 = new TF1("f2","[0]*(x-[1])*(x-[1])+[2]",0,1);
                syoki_para(spgraph,f2,0);
                rep(k,5)spgraph -> Fit(f2,"MQ","",0,1);
                //void rescale_graph(TGraphErrors* sgraph,TGraphErrors* resgraph,TF1* f,int& pos,double yrange,int bin,double freq[nbin])
                rescale_graph(spgraph,respgraph,f2,pos,yran,bin,Freq1,white_hist);
            }
            TF1* gaus1 = new TF1("gaus","[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",-5,5);
            gaus1 -> SetParameter(0,4000);
            gaus1 -> SetParameter(1,0.001);
            gaus1 -> SetParameter(2,0.1);
            white_hist -> Fit(gaus1);
            double error = gaus1 -> GetParameter(2);
            cout << fixed;
            cout << "error is " << setprecision(7) <<error << endl; 
            rep(bin,pos)respgraph -> SetPointError(bin,0,error);
            setting.GraphErrors(respgraph,ifmin,ifmax,-5,5);
            respgraph -> SetTitle("white_noise;Freq[GHz];white_noise[K]");
            respgraph -> Draw("AP");
            //二回目のデータの処理
            /*for(int bin=sb;bin<fb;bin+=dbin){
                TGraphErrors* spgraph = new TGraphErrors;
                TGraphErrors* respgraph = new TGraphErrors;
                make_scale(spgraph,pgraph2,bin-sb);
                bool hantei = false;
                prep(k,bin,bin+dbin){
                    if(c_toge2[j][k] || h_toge2[j][k]){
                        hantei = true;
                        break;
                    }
                }
                if(hantei)continue;
                
                TF1* f2 = new TF1("f2","[0]*(x-[1])*(x-[1])+[2]",0,1);
                syoki_para(spgraph,f2,0);
                rep(k,5)spgraph -> Fit(f2,"MQ","",0,1);
                
            }*/

        }

        /*filesystem::current_path(savedir);
        string gname = "chi_hist"+to_string(i)+".ps";
        c1 -> SaveAs(gname.c_str());
        filesystem::current_path(cdir);*/
    }
    
}