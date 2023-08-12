#include <iostream>
#include <queue>
#include "../headers/setting.h"
#include "../headers/fitter.h"
using namespace std;
typedef long long ll;
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
const double kb=1.38*pow(10,-23);
const double df=88.5*pow(10,3);
const double Tc=76;
const double Th=297;

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
vector<vector<long double>> calibrate2(double hot[nbin],double cold[nbin], double mirror[nbin]){
    vector<vector<long double>> calresult(3,vector<long double>(nbin));//0:gain, 1:psys, 2:prec
    rep(bin,nbin){
        long double gain,psys,prec;
        gain=(hot[bin]-cold[bin])/(2*kb*(Th-Tc)*df);
        psys=(cold[bin]/gain)-2*kb*Tc*df;
        prec=((mirror[bin]/gain)-psys)/(2*kb*df);
        calresult[0][bin] = gain;
        calresult[1][bin] = psys;
        calresult[2][bin] = prec;
    }
    return calresult;
}
void basic_spec2(){
    //お絵描きの設定など
    Double_t xlo=215;
    Double_t xhi=264;
    Double_t ylo=0;
    Double_t yhi=100;
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    TH1F *frame=gPad->DrawFrame(xlo,ylo,xhi,yhi);
    //TH1F* myHist = new TH1F("my Hist","Prec_distribution",100,-100.,300.);
    TGraph* after_fit = new TGraph;//ホワイトノイズフィットした後のノイズを載せる
    ll pnum=0;//トータルのポイントNo.
    //共通の物理量
    double kb=1.38*pow(10,-23);
    double df=88.5*pow(10,3);
    double Tc=76;
    double Th=297;
    int nbin=32767;
    filesystem::path path=filesystem::current_path();
    string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
    string savedird = "/Users/oginokyousuke/data/basic_data/";
    vector<int> sbin={0,512,-512,-1024};
    filesystem::current_path(dir);
    Setting st;
    Fitter ft;
    st.dot_size = 0.6;
    st.markerstyle = 20;
    st.color = kGreen;
    queue<int> que;
    for(int i=1;i<25;i++){
        string cdir=dir+"band"+to_string(i);
        filesystem::current_path(cdir);
        double f[nbin][4],prec[nbin][4],dprec[nbin][4];
        double ifmin = 213.9+2*i;
        double ifmax = 216.1+2*i;
        vector<vector<int>> Ite(34400);
        axrange ax = {ifmin,ifmax,-10,100,0,1,""};
        axrange axg = {ifmin,ifmax,0,pow(10,32),0,1,""};
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
            
            vector<vector<long double>> calres1 = calibrate2(hot1,cold1,mirror1);
            vector<vector<long double>> calres2 = calibrate2(hot2,cold2,mirror2);
            TGraphErrors* pgraph1 = new TGraphErrors;
            TGraphErrors* pgraph2 = new TGraphErrors;
            TGraph* ggraph1 = new TGraph;
            TGraph* ggraph2 = new TGraph;
            prep(bin,sb,fb){
                //cout << bin << " " << Freq1[bin] << " " << prec1[bin] << endl;
                pgraph1 -> SetPoint(bin-sb,Freq1[bin],calres1[2][bin]);
                //cout << calres1[0][bin] << endl;
                pgraph1 -> SetPointError(bin-sb,0,0.1);
                pgraph2 -> SetPoint(bin-sb,Freq2[bin],calres2[2][bin]);
                pgraph2 -> SetPointError(bin-sb,0,0.1);
                ggraph1 -> SetPoint(bin-sb,Freq1[bin],calres1[0][bin]);
                ggraph2 -> SetPoint(bin-sb,Freq2[bin],calres2[0][bin]);
            }
            axg.title = "Gain1;Freq[GHz];Gain[a.u]";
            //c1 -> SetLogy();
            st.Graph(ggraph1,axg);
            ggraph1 -> Draw("AL");
            filesystem::current_path(savedird);
            string gname = "gain"+to_string(i)+"_"+to_string(j)+"_1.ps";
            c1 -> SaveAs(gname.c_str());

            axg.title = "Gain2;Freq[GHz];Prec[K]";
            st.Graph(ggraph2,axg);
            ggraph2 -> Draw("AL");
            //c1 -> SetLogy();
            gname = "gain"+to_string(i)+"_"+to_string(j)+"_2.ps";
            c1 -> SaveAs(gname.c_str());
            filesystem::current_path(cdir);
        }
    }
}