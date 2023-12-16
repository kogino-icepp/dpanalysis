//スプリアス検出コード並びにピークサーチのコードが実データ内でしっかり機能しているかどうかを、テストシグナルを用いて検証する
#include <iostream>
#include "../headers/fitter.h"
#include <random>
using namespace std;
const int nbin=32767;
const int sb = 2585;//切り出してくるビンの最初
const int cb = 15725;//探索範囲の中点,ここを境にデータの振幅が変わっているデータがある
const int fb = 28865;//切り出してくるビンの最後
vector<int> sbin = {0,512,-512,-1024};
const vector<int> lsbo = {1,3,5,13,15,17};
const vector<int> lsbi = {2,4,6,14,16,18};
const vector<int> usbi = {7,9,11,19,21,23};
const vector<int> usbo = {8,10,12,20,22,24};
string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
void toge_scan2(double input[nbin],double &sigma,double (&tvalue)[nbin]){
    random_device rnd;
    mt19937 mt(rnd());
    uniform_int_distribution<> rand1000(0, 9999);  
    string hname = "hist"+to_string(rand1000(mt));
    TH1D* hist = new TH1D(hname.c_str(),"hist;ddinput;Count",100,-0.002,0.002);
    TF1* f = new TF1("gaus","gaus",-0.002,0.002);
    double binput = input[sb-1];
    double bbinput = input[sb-2];
    double ddinput;
    vector<double> togevalue(nbin,0);
    for(int bin=sb;bin<fb;bin++){
        ddinput = input[bin] + bbinput -2*binput;
        ddinput /= input[bin];
        togevalue[bin] = ddinput;
        hist -> Fill(ddinput);
        bbinput = binput;
        binput = input[bin];
    }
    //hist -> Draw();
    hist -> Fit(f,"QEN","",-0.005,0.005);
    //st.Hist(hist);
    sigma = f -> GetParameter("Sigma");
    int num = 0;
    for(int bin=sb;bin<fb;bin++){
        tvalue[bin] = togevalue[bin]/sigma;
        num++;
    }
    //normhist -> Fit("gaus","Q");
}
void fitcodetest(){
    FitFunc ff;
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    /*
    要件定義(本関数内でやりたいこと)
    1. 実データからランダムかつ十分な量のサンプルを取り出す(ここの有意性はしっかり検証)
    2. 取り出した区間に対してピークのサイズを変えながらノイズやシグナルを載せ、それぞれ棘検出とピーク検出がどのくらいできるか実験
    */
    random_device rnd;// 非決定的な乱数生成器を生成
    mt19937 mt(rnd());
    uniform_int_distribution<> randband(1,24);//取り出すバンドをランダムに決定
    uniform_int_distribution<> randj(0,7);//取り出す測定番号を決定
    uniform_int_distribution<> randbin(sb,fb);//どの区間を引っ張り出すかを決定
    Setting st;
    st.dot_size = 0.8;
    st.color = kBlue;
    st.lcolor = kBlue;
    //まずは棘検出がどの程度できるか
    TGraphErrors* sigmagraph = new TGraphErrors;
    prep(pv,1,11){
        double Pgiven = pv*pow(10,13);
        double Psum = 0;
        double dPsum = 0;
        int itenum = 1000;
        double Pite[itenum];
        rep(ite,itenum){
            int i = randband(mt);
            int j0 = randj(mt);
            int mbin = randbin(mt);//質量に対応するbin
            int j = j/2;
            int p = 1+j%2;
            TGraph* testgraph = new TGraph;
            filesystem::path path=filesystem::current_path();
            string cdir=dir+"band"+to_string(i);
            filesystem::current_path(cdir);
            //とりあえずcoldだけ取り出せたらラッキー
            double Freq[nbin],cold[nbin];
            double testcold[nbin];
            for(int data=0;data<576;data++){
                string file_name="test"+to_string(data)+".root";
                if(filesystem::exists(file_name)){
                    TFile* file = new TFile(file_name.c_str());
                    TTree* tree = (TTree*)file->Get("tree");
                    int time,xfft,shift,body,num;
                    long double freq,power;
                    tree->SetBranchAddress("freq",&freq);
                    tree->SetBranchAddress("power",&power);
                    tree->SetBranchAddress("num",&num);
                    tree->SetBranchAddress("body",&body);
                    tree->SetBranchAddress("shift",&shift);
                    tree->GetEntry(0);
                    if(shift==sbin.at(j)){
                        //cout << file_name << endl;
                        if(num==(2-p) && body==0){
                            for(int bin=0;bin<nbin;bin++){
                                tree -> GetEntry(bin);
                                Freq[bin]=freq;
                                cold[bin]=power;
                            }
                        }
                    }
                    
                }
            }
            rep(bin,nbin){
                testcold[bin] = cold[bin]+ff.F_sig2(Freq[bin],Freq[mbin],Pgiven,0.5);
                testgraph -> SetPoint(bin,Freq[bin],testcold[bin]);
            }
            double tsigma;
            double sigvalue[nbin];
            toge_scan2(testcold,tsigma,sigvalue);
            if(i%2==1)Pite[ite] = abs(sigvalue[mbin+2]);
            else Pite[ite] = abs(sigvalue[mbin]);
            double fmin = 213.8+2*i;
            double fmax = 216.2+2*i;
            Psum += Pite[ite];
            /*axrange ax = {fmin,fmax,0,pow(10,16),0,1,";Freq[GHz];Power[a.u.]"};
            st.Graph(testgraph,ax);
            testgraph -> Draw("AL");*/
        }
        Psum /= itenum;
        rep(ite,itenum)dPsum += pow(Psum-Pite[ite],2);
        dPsum = sqrt(dPsum/(itenum-1));
        sigmagraph -> SetPoint(pv,Pgiven,Psum);
        sigmagraph -> SetPointError(pv,0,dPsum);
        cout << Psum << " +- " << dPsum << endl;
    }
    axrange axsigma = {0,pow(10,14),0,100,0,1,";P_{given};ddP_{given}"};
    st.GraphErrors(sigmagraph,axsigma);
    sigmagraph -> Draw("AP");
}