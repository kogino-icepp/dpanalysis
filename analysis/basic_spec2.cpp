#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <format>
#include <tuple>
#include "../headers/fitter.h"

using namespace std;
#define rep(i,n) for(int i=0;i<n;i++)
#define prep(i,m,n) for(int i=m;i<n;i++)
#define fore(i,v) for(auto& i:v)
double PValue(int k,double x){
    return 1-k*TMath::Gamma(k/2,x/(2*k));
}
//特定のデータだけ引っ張ってくる用
vector<int> lsbo = {1,3,5,13,15,17};
vector<int> lsbi = {2,4,6,14,16,18};
vector<int> usbi = {7,9,11,19,21,23};
vector<int> usbo = {8,10,12,20,22,24};
vector<int> lcolor = {4,2,3,6,7,5};
vector<int> sbin = {0,512,-512,-1024};
vector<double> errorv = {0.063,0.066,0.067,0.078,0.066,0.087,0.067,0.062,0.073,0.061,0.089,0.061,0.091,0.142,0.082,0.097,0.087,0.107,0.131,0.087,0.093,0.077,0.103,0.096};
//共通の物理量
const double dcsigma = 0.000487;
const double dhsigma = 0.000457;
const double ddcsigma = 0.000508;
const double ddhsigma = 0.000496;
const double ddmsigma = 0.0004986;
const double dclim = 6.2;
const double dhlim = 6.2;
const double ddclim = 10;
const double ddhlim = 10;
const double ddmlim = 10;
const int nbin=32767;
const int ssb = 2621;//探索すべきビンの最初
const int sfb = 28835;//探索すべきビンの最後
const int sb = 2585;//切り出してくるビンの最初
const int cb = 15725;//探索範囲の中点,ここを境にデータの振幅が変わっているデータがある
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

int both = 0;
int tsafe = 0;
int tout = 0;
void toge_scan(bool (&hantei)[nbin],double input[nbin],double sigma,double limit,TGraph* graph,double freq[nbin]){
    double binput = input[sb-1];
    double bbinput = input[sb-2];
    double ddinput;
    int num = 0;
    for(int bin=sb;bin<fb;bin++){
        ddinput = input[bin] + bbinput - 2*binput;
        ddinput /= input[bin]*sigma;
        graph -> SetPoint(num,freq[bin],ddinput);
        num++;
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
Double_t MyFunction(double x,double p0,double p1){
    return p1*p0*TMath::Exp(-p1*x/2)*pow(p1*x,(p1/2)-1)/(pow(2,p1/2)*TMath::Exp(TMath::LnGamma(p1/2)));
}

//パラメータ数を可変にしたときのフィッティング関数 できればフィット後のステータスも知りたい
/*
要件定義
1. coldとhotでアウトなものは無条件に抜く
2. mirrorに関しては4回分データを走査し、
*/
//棘の始まりと終わりをここに記録しておく、一旦正直な区間を書く、for文にすれば定義も手間じゃない
vector<pair<int,int>> toge_region = {{6626,6629},{9175,9177},{9757,9760},{15730,15730},{16384,16386},{21437,21437},{22356,22356},{22746,22746},{22753,22754},{22877,22878},{22881,22882},{22887,22890},{22896,22898},{22999,23002},{23005,23023},{23114,23119},{23125,23127},{23130,23136},{23141,23146},{23222,23225},{23269,23270},{23277,23277},{23294,23297},{23405,23406},{24576,24578},{26143,26144},{nbin,nbin}};

void basic_spec2(){
    //種々のヘッダー関数を用意
    Setting st;
    st.dot_size = 0.8;
    st.markerstyle = 20;
    st.color = kBlue;
    st.lcolor = kBlue;
    Fitter ft;
    CheckData cda;
    //お絵描きの設定など
    Double_t xlo=215;
    Double_t xhi=264;
    Double_t ylo=0;
    Double_t yhi=100;

    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    /*c1 -> Divide(1,2);
    c1 -> cd(1);*/
    TH1F *frame=gPad->DrawFrame(xlo,ylo,xhi,yhi);
    filesystem::path path=filesystem::current_path();
    string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
    string saveexe = "/Users/oginokyousuke/data/search_exe/";

    filesystem::current_path(dir);

    TGraph* g_ddcm = new TGraph;
    TGraph* g_ddhm = new TGraph;
    int dcnum = 0;
    int dhnum = 0;

    //普遍的な関数(ChiSqure/Ndf)
    
    TF1* chi_f = new TF1("chi_f","MyFunction(x,[0],[1])",0,10);
    
    queue<double> ratio;
    double rsum = 0;
    axrange axscale = {0,1,0,1,0,0,"After Scale;xscale;yscale"};
    vector<pair<double,double>> pairsigma;
    
    for(int i=2;i<3;i++){
        double dym = 0;
        int outnum = 0;
        
        vector<tuple<int,int,double>> revresult;
        double good_num = 0;
        double bad_num = 0;
        cout << i << endl;
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
        rep(j,4){
            rep(k,nbin){
                c_toge1[j][k] = false;
                c_toge2[j][k] = false;
                h_toge1[j][k] = false;
                h_toge2[j][k] = false;
                m_toge1[j][k] = false;
                m_toge2[j][k] = false;
            }
        }
        //フィットの改善を表すグラフ
        TGraph* gchange1 = new TGraph;
        TGraph* gchange2 = new TGraph;
        int gcbin1 = 0;
        int gcbin2 = 0;
        //要件定義: カイ二乗分布でフィットするのがなぜかうまく行かない理由を探る
        for(int j=1;j<2;j++){
            filesystem::current_path(cdir);
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
            
            double ddc1[nbin],ddc2[nbin],ddh1[nbin],ddh2[nbin],ddm1[nbin],ddm2[nbin];
            TGraph* ddcold1 = new TGraph;
            TGraph* ddcold2 = new TGraph;
            TGraph* ddhot1 = new TGraph;
            TGraph* ddhot2 = new TGraph;
            TGraph* gcold1 = new TGraph;
            TGraph* gcold2 = new TGraph;
            TGraph* ggraph1 = new TGraph;
            TGraph* ggraph2 = new TGraph;
            TGraph* sysgraph1 = new TGraph;
            TGraph* sysgraph2 = new TGraph;
            TGraphErrors* pgraph1 = new TGraphErrors;
            TGraphErrors* pgraph2 = new TGraphErrors;
            TGraph* togecold1 = new TGraph;
            TGraph* togecold2 = new TGraph;
            TGraph* ghot1 = new TGraph;
            TGraph* gmirror1 = new TGraph;
            toge_scan(c_toge1[j],cold1,ddcsigma,ddclim,ddcold1,Freq1);
            toge_scan(c_toge2[j],cold2,ddcsigma,ddclim,ddcold2,Freq2);
            toge_scan(h_toge1[j],hot1,ddhsigma,ddhlim,ddhot1,Freq1);
            toge_scan(h_toge2[j],hot2,ddhsigma,ddhlim,ddhot2,Freq2);
            TH1D* ddhist = new TH1D("ddhist","ddhist;ddcold/cold;Count",100,-10,10);
            TH1D* dhist = new TH1D("dhist","dhist;dcold/cold;Count",100,-10,10);
            
            int tnum=0;
            double gain1,gain2,psys1,psys2;
            rep(bin,nbin){
                gcold1 -> SetPoint(bin,Freq1[bin],cold1[bin]);
                ghot1 -> SetPoint(bin,Freq1[bin],hot1[bin]);
                gmirror1 -> SetPoint(bin,Freq1[bin],mirror1[bin]);
                gain1=(hot1[bin]-cold1[bin])/(2*kb*(Th-Tc)*df);
                psys1=(cold1[bin]/gain1)-2*kb*Tc*df;
                gain2=(hot2[bin]-cold2[bin])/(2*kb*(Th-Tc)*df);
                psys2=(cold2[bin]/gain2)-2*kb*Tc*df;
                sysgraph1 -> SetPoint(bin,Freq1[bin],psys1);
                sysgraph2 -> SetPoint(bin,Freq2[bin],psys2);
                ggraph1 -> SetPoint(bin,Freq1[bin],gain1);
                ggraph2 -> SetPoint(bin,Freq2[bin],gain2);
            }
            axrange axraw = {Freq1[23753],Freq1[23694],0,1*pow(10,16),0,1,"Mirror;Freq[GHz];Mirror[a.u.]"};
            axrange axg = {Freq1[23753],Freq1[23694],0,1*pow(10,30),0,1,"Gain;Freq[GHz];Gain[a.u.]"};
            axrange axsys = {Freq1[23753],Freq1[23694],0,1*pow(10,0),0,1,"Psys;Freq[GHz];Psys[a.u.]"};
            axrange axtoge = {Freq1[23753],Freq1[23694],-10,10,0,1,"ddhot;Freq[GHz];ddhot[a.u.]"};
            st.Graph(gmirror1,axraw);
            st.Graph(ghot1,axraw);
            st.Graph(gcold1,axraw);
            st.Graph(ggraph1,axg);
            st.Graph(sysgraph1,axsys);
            st.Graph(togecold1,axraw);
            st.Graph(ddcold1,axtoge);
            st.Graph(ddhot1,axtoge);
            gmirror1 -> SetMarkerColor(kGreen);
            ghot1 -> SetLineWidth(5);
            gmirror1 -> Draw("AP");
            ggraph1 -> SetMarkerColor(kBlack);
            ddhot1 -> SetLineColor(kRed);
            ddhot1 -> Draw("AL");
            //ggraph1 -> Draw("AP");
            sysgraph1 -> SetMarkerColor(kBlack);
            //sysgraph1 -> Draw("AP");
            //togecold1 -> Draw("P");
            /*
            TGraph* dmdg1 = new TGraph;
            TGraph* dmdg2 = new TGraph;
            //評価の仕方がわからないのでとりあえず差分取ってフロアで割る
            TGraph* g1g2 = new TGraph;
            double gain1,psys1,gain2,psys2;
            //とりあえずコンシステントかどうかだけは眺めてみる
            double fgain = 100;
            int num = 0;
            vector<double> prec(nbin);
            rep(bin,nbin){
                dmdg1 -> SetPoint(bin,Freq1[bin],mirror1[bin]/gain1);
                dmdg2 -> SetPoint(bin,Freq2[bin],mirror1[bin]/fgain);
                if(!c_toge1[j][bin-1] && !h_toge1[j][bin-1]){
                    pgraph1 -> SetPoint(num,Freq1[bin],(mirror1[bin]/gain1-psys1)/(2*kb*df));
                    num++;
                }
                prec[bin] = (mirror1[bin]/gain1-psys1)/(2*kb*df);
                //pgraph1 -> SetPoint(bin,Freq1[bin],(mirror1[bin]/gain2-psys2)/(2*kb*df));
                pgraph2 -> SetPoint(bin,Freq2[bin],(mirror1[bin]/gain2-psys2)/(2*kb*df));
                //g1g2 -> SetPoint(bin,Freq1[bin],(gain1-gain2)/(gain1+gain2));

                //if(h_toge1[j][bin])cout << "hot out " << bin << endl;
                //prec[bin]=((mirror[bin]/gain)-psys)/(2*kb*df);
                //pgraph ->SetPoint(bin,Freq1[bin],prec1);
            }
            axrange axg= {Freq1[23028],Freq1[23057],0,pow(10,31),0,1,"Gain;Freq[GHz];Gain[a.u.]"};
            axrange axsys= {Freq1[23028],Freq1[23057],0,300,0,1,"Psys;Freq[GHz];Psys[K]"};
            axrange axdmdg = {Freq1[23255],Freq1[23284],0,pow(10,-15),0,1,"mirror/gain;Freq[GHz];mirror/Gain[a.u.]"};
            axrange axp= {ifmin,ifmax,40,70,0,1,";Freq[GHz];Prec[K]"};
            axrange axg1g2 = {Freq1[23255],Freq1[23284],-1,1,0,1,"(Gain1-Gain2)/(Gain1+Gain2);Freq[GHz];(g1-g2)/(g1+g2)"};
            st.Graph(ggraph1,axg);
            st.Graph(pgraph2,axp);
            st.Graph(pgraph1,axp);
            pgraph1 -> SetLineColor(kGreen);
            pgraph2 -> SetMarkerColor(kRed);
            ggraph2 -> SetMarkerColor(kRed);
            pgraph1 -> SetLineWidth(4);
            pgraph1 -> Draw("AC");
            int sn = 0;
            bool first = true;
            for(auto v:toge_region){
                TGraph* togegraph = new TGraph;
                prep(bin,v.first,v.second+1){
                    togegraph -> SetPoint(bin-v.first,Freq1[bin],prec[bin]);
                }
                st.Graph(togegraph,axp);
                togegraph -> SetLineColor(kRed);
                togegraph -> SetLineWidth(5);
                togegraph -> Draw("L");
            }
            //pgraph2 -> Draw("P");
            //dmdg2 -> Draw("P");
            //pgraph2 -> Draw("P");
            //axg.title = "Gain1;Freq[GHz];Gain[a.u]";
            
            //c1 -> SetLogy();
            /*int testbin = 3035;
            for(int bin=testbin;bin<testbin+1;bin+=dbin){
                bool hantei = false;
                
                double yscale = 100000;
                TGraphErrors* spgraph = new TGraphErrors;
                ft.make_scale(spgraph,pgraph1,bin-sb,yscale);
                double res1,res2;
                TF1* f2 = new TF1("f2","[0]*(x-[1])*(x-[1])+[2]",0,1);
                st.GraphErrors(spgraph,axscale);
                spgraph -> Draw("AP");
                ft.all_fit(spgraph,f2,5,res2);
                f2 -> Draw("same");
                
            }*/
        }
    }
}