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
const double ddclim = 5;
const double ddhlim = 5;
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
void toge_scan(bool (&hantei)[nbin],double input[nbin],double sigma,double limit,TGraph* graph){
    double binput = input[sb-1];
    double bbinput = input[sb-2];
    double ddinput;
    int num = 0;
    for(int bin=sb;bin<fb;bin++){
        ddinput = input[bin] + bbinput - 2*binput;
        ddinput /= input[bin]*sigma;
        graph -> SetPoint(num,num,ddinput);
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
void make5tree(){
    //種々のヘッダー関数を用意
    Setting st;
    st.dot_size = 0.8;
    st.markerstyle = 20;
    st.color = kGreen;
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
    
    for(int i=5;i<6;i++){
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
        for(int j=1;j<4;j++){
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
            TGraphErrors* pgraph1 = new TGraphErrors;
            TGraphErrors* pgraph2 = new TGraphErrors;
            
            prep(bin,sb,fb+30){
                pgraph1 -> SetPoint(bin-sb,Freq1[bin],prec1[bin]);
                pgraph1 -> SetPointError(bin-sb,0,0.1);
                pgraph2 -> SetPoint(bin-sb,Freq2[bin],prec2[bin]);
                pgraph2 -> SetPointError(bin-sb,0,0.1);
            }
            double ddc1[nbin],ddc2[nbin],ddh1[nbin],ddh2[nbin],ddm1[nbin],ddm2[nbin];
            TGraph* ddcold1 = new TGraph;
            TGraph* ddcold2 = new TGraph;
            TGraph* ddhot1 = new TGraph;
            TGraph* ddhot2 = new TGraph;
            toge_scan(c_toge1[j],cold1,ddcsigma,ddclim,ddcold1);
            toge_scan(c_toge2[j],cold2,ddcsigma,ddclim,ddcold2);
            toge_scan(h_toge1[j],hot1,ddhsigma,ddhlim,ddhot1);
            toge_scan(h_toge2[j],hot2,ddhsigma,ddhlim,ddhot2);
            /*axrange axdd = {0,nbin,-10,10,0,1,"ddhot1;bin;ddhot[a.u]"};
            st.Graph(ddhot1,axdd);
            ddhot1 -> Draw("AP");
            prep(bin,sb,fb){
                if(c_toge1[j][bin])cout << "cbin1 : " << bin << endl;
                if(c_toge2[j][bin])cout << "cbin2 : " << bin << endl;
                if(h_toge1[j][bin])cout << "hbin1 : " << bin << endl;
                if(h_toge2[j][bin])cout << "hbin2 : " << bin << endl;
                
            }*/
            //後々のためにフィットで得られた諸データをrootファイルで保存したい
            filesystem::current_path(saveexe);
            string tfname = "all_basekai"+to_string(j)+".root";
            TFile *fout = new TFile(tfname.c_str(),"recreate");
            
            prep(offset,0,30){
                cout << "offset : " << offset << endl;
                string treeName = "tree" + to_string(offset);
                TTree *tree = new TTree(treeName.c_str(), ("tree"+to_string(offset)).c_str());
                double aF,bF,cF,chiF,freqF;
                int binF;
                //chiはndfで割った後の値
                tree -> Branch("a",&aF,"a/D");
                tree -> Branch("b",&bF,"b/D");
                tree -> Branch("c",&cF,"c/D");
                tree -> Branch("chi",&chiF,"chi/D");
                tree -> Branch("bin",&binF,"bin/I");
                tree -> Branch("freq",&freqF,"freq/D");
                for(int bin=sb;bin<fb;bin+=dbin){
                    //cout << bin << endl;
                    axrange axtest = {0,1,0,1,0,1,"test;xscale;yscale"};
                    bin += offset;
                    bool hantei = false;
                    prep(k,bin,bin+dbin){
                        if(c_toge1[j][k] || h_toge1[j][k]){
                            hantei = true;
                            break;
                        }
                    }
                    if(hantei){
                        bin -= offset;
                        continue;
                    }
                    double yscale = 100000;
                    TGraphErrors* spgraph = new TGraphErrors;
                    ft.make_scale(spgraph,pgraph1,bin-sb,yscale);
                    TGraphErrors* spgraphk = new TGraphErrors;
                    rep(spbin,dbin){
                        if(spbin<10 || spbin>=20){
                            double x = spgraph -> GetPointX(spbin);
                            double y = spgraph -> GetPointY(spbin);
                            double ye = spgraph -> GetErrorY(spbin);
                            spgraphk -> SetPoint(spbin,x,y);
                            spgraphk -> SetPointError(spbin,0,ye);
                        }
                    }
                    double res1,res2;
                    TF1* f1 = new TF1("f1","[0]*(x-[1])*(x-[1])+[2]",0,1);
                    //TF1* f2 = new TF1("f2","[0]*(x-[1])*(x-[1])+[2]",0,1);
                    ft.allfit2(spgraphk,f1,5,res1);
                    //ft.exfit(spgraph,f2,res2);
                    //cout << res1 << " : " << res2 << endl;
                    st.GraphErrors(spgraphk,axtest);
                    spgraphk -> Draw("AP");
                    
                    f1 -> Draw("same");
                    //ft.all_fit(spgraph,f2,5,res2);
                    aF = f1 -> GetParameter(0);
                    bF = f1 -> GetParameter(1);
                    cF = f1 -> GetParameter(2);
                    chiF = res1;
                    binF = bin;
                    freqF = Freq1[bin];
                    //cout << "parameter's are" << endl;
                    //cout << aF << " " << bF << " " << cF << endl;
                    //cout << res1 << " " << freqF << endl;
                    tree -> Fill();
                    bin-=offset;
                }
                
            }
            fout -> Write();
            fout -> Close();
            
            
            double ysmax = -10;
            //TH1D* ketah = new TH1D("hetah","ketah;PrecisionKeta;Count",21,0,20);
            //一回目のデータの処理
            int offbin = 0;
            
            
            //tree1->SaveAs(tfname1.c_str());
            /*tree1 -> Write();
            fout1 -> Close();*/

            /*string tfname2 = "fit_res"+to_string(i)+"_"+to_string(j)+"_2.root";
            TFile *fout2 = new TFile(tfname2.c_str(),"recreate");
            TTree *tree2 = new TTree("tree2","tree2");
            double aF2,bF2,cF2,chiF2,freqF2;
            int binF2;
            tree2 -> Branch("a",&aF2,"a/D");
            tree2 -> Branch("b",&bF2,"b/D");
            tree2 -> Branch("c",&cF2,"c/D");
            tree2 -> Branch("chi",&chiF2,"chi/D");
            tree2 -> Branch("bin",&binF2,"bin/I");
            tree2 -> Branch("freq",&freqF2,"freq/D");
            //double ysmax = -10;
            //TH1D* ketah = new TH1D("hetah","ketah;PrecisionKeta;Count",21,0,20);
            //2回目のデータの処理
            for(int bin=sb;bin<fb;bin+=dbin){
                cout << 2 << " " << bin << endl;
                bool hantei = false;
                prep(k,bin,bin+dbin){
                    if(c_toge2[j][k] || h_toge2[j][k]){
                        hantei = true;
                        break;
                    }
                }
                if(hantei){
                    continue;
                }
                double yscale = 100000;
                TGraphErrors* spgraph = new TGraphErrors;
                ft.make_scale(spgraph,pgraph2,bin-sb,yscale);
                double res2;
                TF1* f2 = new TF1("f2","[0]*(x-[1])*(x-[1])+[2]",0,1);
                
                ft.all_fit(spgraph,f2,5,res2);
                aF2 = f2 -> GetParameter(0);
                bF2 = f2 -> GetParameter(1);
                cF2 = f2 -> GetParameter(2);
                chiF2 = res2;
                binF2 = bin;
                freqF2 = Freq2[bin];
                tree2 -> Fill();
                chi_hist -> Fill(res2);
                cout << fixed;
                cout << freqF2 << endl;
                cout << aF2 << " : " << bF2 << " : " << cF2 << endl;
            }
            //tree2->SaveAs(tfname2.c_str());
            tree2 -> Write();
            fout2 -> Close();*/
        }
        
    }
}