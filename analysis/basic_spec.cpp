#include <iostream>

#include "../headers/fitter.h"
using namespace std;
const int sb = 2585;//切り出してくるビンの最初
const int cb = 15725;//探索範囲の中点,ここを境にデータの振幅が変わっているデータがある
const int fb = 28865;//切り出してくるビンの最後
pair<double,double> MeanError(vector<double>data){
    double mean = 0;
    double num = 0;
    for(auto v:data){
        mean += v;
        num++;
    }
    mean /= num;
    double rtn = 0;
    for(auto v:data)rtn += (mean-v)*(mean-v);
    rtn = sqrt(rtn/(num-1));
    return {mean,rtn};
}
void basic_spec(){
    Setting st;
    st.color = kBlue;
    st.lcolor = kBlue;
    //お絵描きの設定など
    Double_t xlo=215;
    Double_t xhi=264;
    Double_t ylo=0;
    Double_t yhi=100;
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    c1 -> SetLogy();
    //TH1F* myHist = new TH1F("my Hist","Prec_distribution",100,-100.,300.);

    //共通の物理量
    double kb=1.38*pow(10,-23);
    double df=88.5*pow(10,3);
    double Tc=76;
    double Th=297;
    int nbin=32767;
    filesystem::path path=filesystem::current_path();
    string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
    string savedir = "/Users/oginokyousuke/data/basic_data/";
    vector<int> sbin={0,512,-512,-1024};
    filesystem::current_path(dir);
    TH1D* good_fit= new TH1D("chi2/NDF","chi2/NDF",100,-2,2);
    TH1D* p_value = new TH1D("p_value","p_value;Prob;Count",100,-0.2,1.2);
    TH1D* white_noise = new TH1D("delta_P","delta_P;delta_P[K];Count",100,-20,20);
    TGraph* fit_result = new TGraph;//帯域ごとにどのくらいフィットがうまくいっているのか評価する
    TGraphErrors* spec_ave = new TGraphErrors;
    int avebin = 0;
    axrange axave = {215,265,0,pow(10,16),0,1,";Freq[GHz];Spec[a.u.]"};
    for(int i=1;i<2;i++){
        string cdir=dir+"band"+to_string(i);
        filesystem::current_path(cdir);
        
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
            TGraph* cgraph = new TGraph;
            TGraph* hgraph = new TGraph;
            TGraph* mgraph = new TGraph;
            TGraph* calgraph = new TGraph;
            rep(bin,nbin){
                
            }
            st.Graph(hmgraph,axave);
            hmgraph -> Draw("AL");
        }
        /*prep(bin,sb,fb){
            spec_ave -> SetPoint(avebin,Freq[bin],MeanError(Cold[bin]).first);
            spec_ave -> SetPointError(avebin,0,MeanError(Cold[bin]).second);
            avebin++;
        }*/
    }
    /*st.GraphErrors(spec_ave,axave);
    spec_ave -> Draw("AP");*/
}