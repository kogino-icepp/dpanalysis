#include <iostream>
#include <string>
#include "../headers/fitter.h"
const int nbin=32767;
const int sb = 2585;//切り出してくるビンの最初
const int cb = 15725;//探索範囲の中点,ここを境にデータの振幅が変わっているデータがある
const int fb = 28865;//切り出してくるビンの最後

vector<int> sbin = {0,512,-512,-1024};

void calib_fit(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
    Setting st;
    st.dot_size = 0.8;
    st.markerstyle = 20;
    st.color = kBlue;
    Fitter ft;
    prep(i,1,2){
        string cdir = dir + "band" + to_string(i);
        filesystem::current_path(cdir);
        prep(j,0,1){
            double Freq1[nbin],cold1[nbin],hot1[nbin],mirror1[nbin],Freq2[nbin],cold2[nbin],hot2[nbin],mirror2[nbin];
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
            //試しにベースラインフィットを敢行してみる？
            TGraph* gcold = new TGraph;
            TGraph* ghot = new TGraph;
            TGraph* gmirror = new TGraph;
            rep(bin,nbin){
                gcold -> SetPoint(bin,Freq1[bin],cold1[bin]);
                ghot -> SetPoint(bin,Freq1[bin],hot1[bin]);
                gmirror -> SetPoint(bin,Freq1[bin],mirror1[bin]);
            }
            axrange axscale = {0,1,0,1,0,1,"cold(scale);xscale;yscale"};
            TH1D* chihist = new TH1D("hist","hist;chi/ndf;Count",100,0,2);
            TH1D* sigmahist = new TH1D("sigmahist","sigma;Sigma;Count",100,0,0.5);
            prep(bin,sb,fb){
                TF1* coldfunc = new TF1("coldfunc","[0]*(x-[1])*(x-[1])+[2]",0,1);
                TF1* hotfunc = new TF1("hotfunc","[0]*(x-[1])*(x-[1])+[2]",0,1);
                double yMin,ycscale,yhscale;
                TGraphErrors* scgraph = new TGraphErrors;
                ft.make_scale(scgraph,gcold,bin,ycscale);
                ft.make_scale(shgraph,ghot,bin,yhscale);
                st.GraphErrors(scgraph,axscale);
                sgraph -> Draw("AP");
                double res1,res2;
                ft.calfit(scgraph,coldfunc,5,res1);
                ft.calfit(shscale,hotfunc,5,res2)
                sigmahist -> Fill(res1);
                quadfunc -> Draw("same");
                double a,b,c;
                //フィット結果を用いてキャリブレーションしたらどうなるのか確認する
                chihist -> Fill(res1);
            }
            c1 -> SetLogy();
            st.Hist(sigmahist);
            sigmahist -> Draw();
        }
    }
}