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


vector<int> lsbo = {1,3,5,13,15,17};
vector<int> lsbi = {2,4,6,14,16,18};
vector<int> usbi = {7,9,11,19,21,23};
vector<int> usbo = {8,10,12,20,22,24};
vector<int> lcolor = {4,2,3,6,7,5};
vector<int> sbin = {0,512,-512,-1024};

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
const double DINF=1e9;
const double DeltaP = 0;
const double Aeff = 0.385;

string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
string rdir="/Users/oginokyousuke/data/root_file/";

void sec_deriv(double input[nbin],TH1D*hist,double (&vec)[nbin]){
    double binput = input[sb-1];
    double bbinput = input[sb-2];
    double ddinput;
    for(int bin=sb;bin<fb;bin++){
        ddinput = input[bin] + bbinput - 2*binput;
        ddinput /= input[bin];
        hist -> Fill(ddinput);
        vec[bin] = ddinput;
        bbinput = binput;
        binput = input[bin];
    }
}
void make_togemap(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    Setting st;
    
    for(int i=15;i<16;i++){
        string cdir=dir+"band"+to_string(i);
        filesystem::current_path(cdir);
        //いつもの感じで生データ取り出し、ddの処理までやる
        bool chantei[nbin][8];
        bool hhantei[nbin][8];
        bool mhantei[nbin][8];
        for(int j=0;j<4;j++){
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
            //それぞれの配列に対してddvalueを計算→treeに格納する
            double vddcold[2][nbin],vddhot[2][nbin],vddmirror[2][nbin];
            
            TH1D* cddhist = new TH1D("cddhist","cddhist;ddcold/cold;Count",100,-0.005,0.005);
            TH1D* hddhist = new TH1D("hddhist","hddhist;ddhot/hot;Count",100,-0.005,0.005);
            TH1D* mddhist = new TH1D("mddhist","mddhist;ddmirror/mirror;Count",100,-0.005,0.005);
            // sec_deriv(cold1,cold2,vddcold,cddhist);
            // sec_deriv(hot1,hot2,vddhot,hddhist);
            // sec_deriv(mirror1,mirror2,vddmirror,mddhist);
            sec_deriv(cold1,cddhist,vddcold[0]);
            sec_deriv(cold2,cddhist,vddcold[1]);
            sec_deriv(hot1,hddhist,vddhot[0]);
            sec_deriv(hot2,hddhist,vddhot[1]);
            sec_deriv(mirror1,mddhist,vddmirror[0]);
            sec_deriv(mirror2,mddhist,vddmirror[1]);
            double csigma,hsigma,msigma;
            TF1* fgaus = new TF1("fgaus","gaus",-0.005,0.005);
            st.Hist(cddhist),st.Hist(hddhist),st.Hist(mddhist);
            //c1 -> SetLogy();
            cddhist -> Draw();
            cddhist -> Fit(fgaus,"Q");
            csigma = fgaus -> GetParameter("Sigma");
            //hddhist -> Draw();
            //hddhist -> Fit(fgaus,"Q");
            //hsigma = fgaus -> GetParameter("Sigma");
            //mddhist -> Draw();
            //mddhist -> Fit(fgaus,"Q");
            //msigma = fgaus -> GetParameter("Sigma");
            //cout << csigma << " " << hsigma << " " << msigma << endl;
            
            rep(bin,nbin){
                vddcold[0][bin] /= csigma,vddcold[1][bin] /= csigma;
                vddhot[0][bin] /= hsigma, vddhot[1][bin] /= hsigma;
                vddmirror[0][bin] /= msigma, vddmirror[1][bin] /= msigma;
                if(abs(vddcold[0][bin])>5)chantei[bin][2*j] = true;
                if(abs(vddcold[1][bin])>5)chantei[bin][2*j+1] = true;
                if(abs(vddhot[0][bin])>5)hhantei[bin][2*j] = true;
                if(abs(vddhot[1][bin])>5)hhantei[bin][2*j+1] = true;
                if(abs(vddmirror[0][bin])>5)mhantei[bin][2*j] = true;
                if(abs(vddmirror[1][bin])>5)mhantei[bin][2*j+1] = true;
            }
            //rootファイルは諦めた、とりあえず間に合わせでも比較する。
            //使い回しでいけるところと
            // int tbin,ite,type,toge;
            // double ddvalue;
        }
        //なんとか判定の配列自体は作ったのでこれを比較してグルーピングする
        int toge_type[nbin][8];
        int not_toge = 0;
        int true_toge = 0;
        int detect_error = 0;
        TH1D* type_hist = new TH1D("type_hist","type;type;Count",9,-1,8);
        prep(bin,sb,fb){
            bool hantei = false;
            int x = 0;
            rep(j,8){
                if(chantei[bin][j] && hhantei[bin][j] && mhantei[bin][j])toge_type[bin][j] = 0,true_toge++;//全てで正常に検出できている、これが理想
                else if(!chantei[bin][j] && !hhantei[bin][j] && !mhantei[bin][j])toge_type[bin][j] = 7,not_toge++;//一番安全、そもそも棘がない
                else if(chantei[bin][j] && hhantei[bin][j] && !mhantei[bin][j])toge_type[bin][j] = 1,detect_error++;//1~6ではどこかで検出できていない、この状態がまずい
                else if(chantei[bin][j] && !hhantei[bin][j] && mhantei[bin][j])toge_type[bin][j] = 2,detect_error++;
                else if(!chantei[bin][j] && hhantei[bin][j] && mhantei[bin][j])toge_type[bin][j] = 3,detect_error++;
                else if(!chantei[bin][j] && !hhantei[bin][j] && mhantei[bin][j])toge_type[bin][j] = 4,detect_error++;
                else if(!chantei[bin][j] && hhantei[bin][j] && !mhantei[bin][j])toge_type[bin][j] = 5,detect_error++;
                else if(chantei[bin][j] && !hhantei[bin][j] && !mhantei[bin][j])toge_type[bin][j] = 6,detect_error++;
                type_hist -> Fill(toge_type[bin][j]);
                
            }
            // cout << x << endl;
            //if(x == 0)true_toge++;
            //else if(x==56)not_toge++;
            //else detect_error++;
            //そもそも2測定とものデータなのかどうかを統合する？
        }
        st.Hist(type_hist);
        type_hist -> Draw();
        cout << not_toge << " " << true_toge << " " << detect_error << endl;
    }
}