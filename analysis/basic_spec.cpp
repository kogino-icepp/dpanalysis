#include <iostream>
#include "../headers/setting.h"
#include "../headers/fitter.h"
using namespace std;

void basic_spec(){
    //お絵描きの設定など
    Double_t xlo=215;
    Double_t xhi=264;
    Double_t ylo=0;
    Double_t yhi=100;
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    TH1F *frame=gPad->DrawFrame(xlo,ylo,xhi,yhi);
    //TH1F* myHist = new TH1F("my Hist","Prec_distribution",100,-100.,300.);

    //共通の物理量
    double kb=1.38*pow(10,-23);
    double df=88.5*pow(10,3);
    double Tc=76;
    double Th=297;
    int nbin=32767;
    filesystem::path path=filesystem::current_path();
    string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
    vector<int> sbin={0,512,-512,-1024};
    filesystem::current_path(dir);
    TH1D* good_fit= new TH1D("chi2/NDF","chi2/NDF",100,-2,2);
    TH1D* p_value = new TH1D("p_value","p_value;Prob;Count",100,-0.2,1.2);
    TH1D* white_noise = new TH1D("delta_P","delta_P;delta_P[K];Count",100,-20,20);
    TGraph* fit_result = new TGraph;//帯域ごとにどのくらいフィットがうまくいっているのか評価する
    for(int i=2;i<3;i++){
        TGraph* graph1 = new TGraph;
        string cdir=dir+"band"+to_string(i);
        filesystem::current_path(cdir);
        double f[nbin][4],prec[nbin][4],dprec[nbin][4];
        double fmin=213.5+i*2;
        double fmax=fmin+3;
        for(int j=0;j<4;j++){
            TGraph* graph = new TGraph;
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
            long double cmae1=cold1[1999];
            long double cmae2=cold2[1999];
            long double hmae1=hot1[1999];
            long double hmae2=hot2[1999];
            vector<bool> ch1(nbin,true);
            vector<bool> ch2(nbin,true);
            vector<bool> hh1(nbin,true);
            vector<bool> hh2(nbin,true);
            const double csigma = 0.000487;
            const double hsigma = 0.000456;
            int alpha = 20;
            double clim=csigma*alpha;
            double hlim=hsigma*alpha;
            int dist=5;
            double dcc1,dcc2,dhh1,dhh2;
            //ノーマライズして刺抜き
            for(int bin=2000;bin<nbin-2000;bin++){
                long double dcold1=cold1[bin]-cmae1;
                long double dcold2=cold2[bin]-cmae2;
                long double dhot1=hot1[bin]-hmae1;
                long double dhot2=hot2[bin]-hmae2;
                dcc1=dcold1/cmae1;
                dcc2=dcold2/cmae2;
                dhh1=dhot1/hmae1;
                dhh2=dhot2/hmae2;
                cmae1=cold1[bin];
                cmae2=cold2[bin];
                hmae1=hot1[bin];
                hmae2=hot2[bin];
                if(abs(dcc1)>clim){
                    for(int j=-5;j<5;j++) ch1[bin+j]=false;
                }
                if(abs(dcc2)>clim){
                    for(int j=-5;j<5;j++) ch2[bin+j]=false;
                }
                if(abs(dhh1)>hlim){
                    for(int j=-5;j<5;j++) hh1[bin+j]=false;
                }
                if(abs(dhh2)>hlim){
                    for(int j=-5;j<5;j++) hh2[bin+j]=false;
                }
            }
            TGraph* pgraph1 = new TGraph;
            TGraph* pgraph2 = new TGraph;
            TGraphErrors* pegraph1 = new TGraphErrors;
            TGraphErrors* pegraph2 = new TGraphErrors;
            int bin_true1=0;
            int bin_true2=0;
            long double gain1,gain2;
            long double psys1,psys2;
            long double prec1,prec2;
            for(int bin=2000;bin<nbin-2000;bin++){
                if(!ch1[bin] || !hh1[bin])continue;
                gain1=(hot1[bin]-cold1[bin])/(2*kb*(Th-Tc)*df);
                double psys1=(cold1[bin]/gain1)-2*kb*Tc*df;
                double prec1=((mirror1[bin]/gain1)-psys1)/(2*kb*df);
                pgraph1->SetPoint(bin_true1,bin,prec1);
                pegraph1 -> SetPoint(bin_true1,bin,prec1);
                pegraph1 -> SetPointError(bin_true1,0.0,0.1);
                bin_true1++;
            }
            for(int bin=2000;bin<nbin-2000;bin++){
                if(!ch2[bin] || !hh2[bin])continue;
                gain2=(hot2[bin]-cold2[bin])/(2*kb*(Th-Tc)*df);
                psys2=(cold2[bin]/gain2)-2*kb*Tc*df;
                prec2=((mirror2[bin]/gain2)-psys2)/(2*kb*df);
                pgraph2->SetPoint(bin_true2,bin,prec2);
                pegraph2->SetPoint(bin_true2,bin,prec2);
                pegraph2->SetPointError(bin_true2,0.0,0.1);
                bin_true2++;
            }
            
            vector<vector<double>> fit_param1(3);
            vector<vector<double>> fit_param2(3);
            double p01,p11,p21,p31,p02,p12,p22,p32;
            pegraph1->SetTitle("base_fit_test;bin;Prec[K]");
            pegraph1->SetMaximum(100);
            pegraph1->SetMinimum(0);
            pegraph1->GetXaxis()->SetLimits(0,nbin);
            pegraph1->SetLineColor(kGreen);
            pegraph1->Draw("AC");
            gStyle->SetOptFit();
            /*pgraph2->SetTitle("base_fit_test;bin;Prec[K]");
            pgraph2->SetMaximum(100);
            pgraph2->SetMinimum(0);
            pgraph2->GetXaxis()->SetLimits(0,nbin);
            pgraph2->SetLineColor(kGreen);
            pgraph2->Draw("AC");*/
            double dbin=Freq1[1]-Freq1[0];
            cout << dbin << endl;
            //+-二つでフィットして良い方使うとか？
            double chi11,chi12,chi21,chi22;
            TF1 *f1 = new TF1("base1","[0]*(x-[1])*(x-[1])+[2]");
            TF1 *f2 = new TF1("base2","[0]*(x-[1])*(x-[1])+[2]");
            int sx,gx;
            double sy,gy,ty;
            int len1=0;
            int len2=0;
            /*for(int bin=0;bin<bin_true1;bin+=60){
                if(bin+60>bin_true1){
                    len1=bin;
                    break;
                }
                sx=pgraph1->GetPointX(bin);
                gx=pgraph1->GetPointX(bin+59);
                sy=pgraph1->GetPointY(bin);
                gy=pgraph1->GetPointY(bin+59);
                ty=pgraph1->GetPointY(bin+30);
                if((sy+gy)/2>ty){
                    f1->SetParameter(0,0.01);
                    f1->SetParameter(1,sx+30);
                    f1->SetParameter(2,(sy+gy)/2);
                }
                else{
                    f1->SetParameter(0,-0.01);
                    f1->SetParameter(1,sx+30);
                    f1->SetParameter(2,(sy+gy)/2);
                }
                for(int i=0;i<5;i++){
                    pgraph1->Fit("base1","M","CP",sx,gx);
                    p01=f1->GetParameter(0);
                    p11=f1->GetParameter(1);
                    p21=f1->GetParameter(2);
                    f1->SetParameter(0,p01);
                    f1->SetParameter(1,p11);
                    f1->SetParameter(2,p21);
                }
                pgraph1->Fit("base1","M","CP",sx,gx);
                p01=f1->GetParameter(0);
                p11=f1->GetParameter(1);
                p21=f1->GetParameter(2);
                fit_param1[0].push_back(p01);
                fit_param1[1].push_back(p11);
                fit_param1[2].push_back(p21);
            }
            for(int bin=0;bin<bin_true2;bin+=60){
                if(bin+60>bin_true2){
                    len2=bin;
                    break;
                }
                sx=pgraph2->GetPointX(bin);
                gx=pgraph2->GetPointX(bin+59);
                sy=pgraph2->GetPointY(bin);
                gy=pgraph2->GetPointY(bin+59);
                ty=pgraph2->GetPointY(bin+30);
                if((sy+gy)/2>ty){
                    f2->SetParameter(0,0.01);
                    f2->SetParameter(1,sx+30);
                    f2->SetParameter(2,(sy+gy)/2);
                }
                else{
                    f2->SetParameter(0,-0.01);
                    f2->SetParameter(1,sx+30);
                    f2->SetParameter(2,(sy+gy)/2);
                }
                for(int i=0;i<5;i++){
                    pgraph2->Fit("base2","M","CP",sx,gx);
                    p02=f2->GetParameter(0);
                    p12=f2->GetParameter(1);
                    p22=f2->GetParameter(2);
                    f2->SetParameter(0,p02);
                    f2->SetParameter(1,p12);
                    f2->SetParameter(2,p22);
                }
                pgraph2->Fit("base2","M","CP",sx,gx);
                p02=f2->GetParameter(0);
                p12=f2->GetParameter(1);
                p22=f2->GetParameter(2);
                fit_param2[0].push_back(p02);
                fit_param2[1].push_back(p12);
                fit_param2[2].push_back(p22); 
            }*/
            for(int bin=0;bin<bin_true1;bin+=60){
                if(bin+60>bin_true1){
                    len1=bin;
                    break;
                }
                sx=pegraph1->GetPointX(bin);
                gx=pegraph1->GetPointX(bin+59);
                sy=pegraph1->GetPointY(bin);
                gy=pegraph1->GetPointY(bin+59);
                ty=pegraph1->GetPointY(bin+30);
                if((sy+gy)/2>ty){
                    f1->SetParameter(0,0.01);
                    f1->SetParameter(1,sx+30);
                    f1->SetParameter(2,(sy+gy)/2);
                }
                else{
                    f1->SetParameter(0,-0.01);
                    f1->SetParameter(1,sx+30);
                    f1->SetParameter(2,(sy+gy)/2);
                }
                for(int i=0;i<5;i++){
                    pegraph1->Fit("base1","+","CP",sx,gx);
                    p01=f1->GetParameter(0);
                    p11=f1->GetParameter(1);
                    p21=f1->GetParameter(2);
                    f1->SetParameter(0,p01);
                    f1->SetParameter(1,p11);
                    f1->SetParameter(2,p21);
                }
                pegraph1->Fit("base1","","CP",sx,gx);
                p01=f1->GetParameter(0);
                p11=f1->GetParameter(1);
                p21=f1->GetParameter(2);
                fit_param1[0].push_back(p01);
                fit_param1[1].push_back(p11);
                fit_param1[2].push_back(p21);
            }
            for(int bin=0;bin<bin_true2;bin+=60){
                if(bin+60>bin_true2){
                    len2=bin;
                    break;
                }
                sx=pegraph2->GetPointX(bin);
                gx=pegraph2->GetPointX(bin+59);
                sy=pegraph2->GetPointY(bin);
                gy=pegraph2->GetPointY(bin+59);
                ty=pegraph2->GetPointY(bin+30);
                if((sy+gy)/2>ty){
                    f2->SetParameter(0,0.01);
                    f2->SetParameter(1,sx+30);
                    f2->SetParameter(2,(sy+gy)/2);
                }
                else{
                    f2->SetParameter(0,-0.01);
                    f2->SetParameter(1,sx+30);
                    f2->SetParameter(2,(sy+gy)/2);
                }
                for(int i=0;i<5;i++){
                    pegraph2->Fit("base2","","CP",sx,gx);
                    p02=f2->GetParameter(0);
                    p12=f2->GetParameter(1);
                    p22=f2->GetParameter(2);
                    f2->SetParameter(0,p02);
                    f2->SetParameter(1,p12);
                    f2->SetParameter(2,p22);
                }
                pegraph2->Fit("base2","","CP",sx,gx);
                p02=f2->GetParameter(0);
                p12=f2->GetParameter(1);
                p22=f2->GetParameter(2);
                fit_param2[0].push_back(p02);
                fit_param2[1].push_back(p12);
                fit_param2[2].push_back(p22); 
            }
            
            int paranum1=-1;
            int paranum2=-1;
            for(int bin=0;bin<len1;bin++){
                if(bin%60==0)paranum1++;
                int x=pegraph1->GetPointX(bin);
                double prec1=pegraph1->GetPointY(bin);
                double cprec1=fit_param1[0][paranum1]*(x-fit_param1[1][paranum1])*(x-fit_param1[1][paranum1])+fit_param1[2][paranum1];
                double dprec1=prec1-cprec1;
                white_noise->Fill(dprec1);
                
            }
            for(int bin=0;bin<len2;bin++){
                if(bin%60==0)paranum2++;
                double prec2=pegraph2->GetPointY(bin);
                int x=pegraph2->GetPointX(bin);
                double cprec2=fit_param2[0][paranum2]*(x-fit_param2[1][paranum2])*(x-fit_param2[1][paranum2])+fit_param2[2][paranum2];
                double dprec2=prec2-cprec2;
                white_noise->Fill(dprec2);
            }
        }
        //fit_result->SetPoint(i-1,(fmax+fmin)/2,)
    }
    //white_noise->Fit("gaus");
    c1->SetLogy();
    white_noise->Draw();
    //p_value->Draw();
}