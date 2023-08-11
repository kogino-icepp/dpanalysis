#include <iostream>
#include <set>
using namespace std;

void cold_toge_mk2(){
    //お絵描きの設定など
    Double_t xlo=215;
    Double_t xhi=264;
    Double_t ylo=0;
    Double_t yhi=100;
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    TH1F *frame=gPad->DrawFrame(xlo,ylo,xhi,yhi);
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(0.7);
    //TH1F* togeHist = new TH1F("my Hist","Prec_distribution",100,0.,500.);
    vector<int> lsbo={1,3,5,13,15,17};
    vector<int> lsbi={2,4,6,14,16,18};
    vector<int> usbi={7,9,11,19,21,23};
    vector<int> usbo={8,10,12,20,22,24};
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
    TH1D* nchist = new TH1D("nchist","dcold/cold;dcold/cold[a.u.];Count",100,-1000,1000);
    TH1D* gnchist = new TH1D("gnchist","dcold/gain;dcold/gain[a.u.];Count",100,-100,100);
    TH1D* nhhist = new TH1D("nhhist","dhot/hot;dhot/hot;Count",100,-0.01,0.01);
    TH1D* gnhhist = new TH1D("gnhhist","dhot/gain;dhot/gain;Count",100,-100,100);
    TH2D* chist2 = new TH2D("chist2","dcold;power[a.u.x10^12];dpower[a.u.x10^12]",100,0,8000,100,-10,10);
    const double csigma = 0.000487;
    const double hsigma = 0.000456;
    const int alpha = 20;
    for(int i=16;i<17;i++){
        string cdir=dir+"band"+to_string(i);
        filesystem::current_path(cdir);
        deque<double> all_dcold;
        double fmin=213.5+i*2;
        double fmax=fmin+3;
        for(int j=2;j<3;j++){
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
            TGraph* cgraph1 = new TGraph;
            TGraph* cgraph2 = new TGraph;
            TGraph* dcgraph1 = new TGraph;
            TGraph* dcgraph2 = new TGraph;
            TGraph* hgraph1 = new TGraph;
            TGraph* hgraph2 = new TGraph;
            TGraph* mgraph1 = new TGraph;
            TGraph* mgraph2 = new TGraph;
            TGraph* omgraph1 = new TGraph;
            TGraph* omgraph2 = new TGraph;

            double fin=min(Freq1[0],Freq1[nbin-1]);
            //cout << fin << endl;
            long double cmae1=cold1[1999];
            long double cmae2=cold2[1999];
            long double hmae1=hot1[1999];
            long double hmae2=hot2[1999];
            long double mmae1=mirror1[1999];
            long double mmae2=mirror2[1999];
            long double dcold1,dcold2,dhot1,dhot2;
            long double dcc1,dcc2,dcg1,dcg2;
            long double dhh1,dhh2,dhg1,dhg2;
            int bin_true1=0;
            int bin_true2=0;
            int bin_false1=0;
            int bin_false2=0;
            long double clim = csigma*alpha;
            long double hlim = hsigma*alpha;

            for(int bin=2000;bin<nbin-2000;bin++){
                dcold1=cold1[bin]-cmae1;
                dcold2=cold2[bin]-cmae2;
                dhot1=hot1[bin]-hmae1;
                dhot2=hot2[bin]-hmae2;

                if(bin==19506)cout << Freq1[bin] << endl;
                dcc1=dcold1/cmae1;
                dcc2=dcold2/cmae2;
                dcg1=dcold1/(hmae1-cmae1);
                dcg2=dcold2/(hmae2-cmae2);
                dhh1=dhot1/hmae1;
                dhh2=dhot2/hmae2;
                dhg1=dhot1/(hmae1-cmae1);
                dhg2=dhot2/(hmae2-cmae2);

                if(dcc1 < clim && dhh1 < hlim){
                    cgraph1 -> SetPoint(bin_true1,Freq1[bin],cold1[bin]);
                    dcgraph1 -> SetPoint(bin_true1,Freq1[bin],dcc1);
                    mgraph1 -> SetPoint(bin_true1,Freq1[bin],mirror1[bin]);
                    bin_true1++;
                }
                else{
                    omgraph1 -> SetPoint(bin_false1,Freq1[bin],cold1[bin]);
                    bin_false1++;
                }
                if(dcc2 < clim && dhh2 < hlim){
                    cgraph2 -> SetPoint(bin_true2,Freq2[bin],cold2[bin]);
                    dcgraph2 -> SetPoint(bin_true2,Freq2[bin],dcc2);
                    mgraph2 -> SetPoint(bin_true2,Freq2[bin],mirror2[bin]);
                    bin_true2++;
                }
                else{
                    omgraph2 -> SetPoint(bin_false2,Freq2[bin],cold2[bin]);
                    bin_false2++;
                }
                cmae1=cold1[bin];
                cmae2=cold2[bin];
                hmae1=hot1[bin];
                hmae2=hot2[bin];
                

                nchist->Fill(dcc1);
                nchist->Fill(dcc2);
                gnchist->Fill(dcg1);
                gnchist->Fill(dcg2);
                nhhist->Fill(dhh1);
                nhhist->Fill(dhh2);
                gnhhist->Fill(dhg1);
                gnhhist->Fill(dhg2);
            }
            cgraph1 -> SetMarkerColor(kBlue);
            cgraph1 -> SetMarkerSize(0.3);
            cgraph2 -> SetMarkerColor(kBlue);
            cgraph2 -> SetMarkerSize(0.3);
            cgraph1 -> SetMaximum(3*pow(10,15));
            cgraph1 -> SetMinimum(5*pow(10,14));
            cgraph1 -> GetXaxis()->SetLimits(246.67,246.68);
            omgraph1 -> SetMarkerColor(kRed);
            omgraph1 -> SetMarkerSize(0.3);
            omgraph2 -> SetMarkerColor(kRed);
            omgraph2 -> SetMarkerSize(0.3);

            if(i==16 && j==2){
                cgraph1 -> SetTitle("cold_power;Freq[GHz];cold_power[a.u.]");
                cgraph1 -> Draw("AP");
            }
            else cgraph1 -> Draw("P");
            cgraph2 -> Draw("P");
            omgraph1 -> Draw("P");
            omgraph2 -> Draw("P");
        }
    }
    //c1->SetLogy();
    //chist2->Draw("colz");
    //nchist->Fit("gaus");
    //nchist->Draw();
    //nhhist->Fit("gaus");
    //nhhist->Draw();
}
//band18,23,24 周波数は正しいけど答えがおかしい