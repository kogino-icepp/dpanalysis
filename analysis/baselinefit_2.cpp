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
Double_t chif(double x,double p0,double p1){
    return p0*TMath::Exp(-x/2)*pow(x,(p1/2)-1)/(pow(2,p1/2)*TMath::Gamma(p1/2));
}
Double_t chif_ndf(double x,double p0,double p1){
    return p1*chif(p1*x,p0,p1);
}
Double_t chiF(double x,double p0,double p1){
    return p0*TMath::Gamma(p1/2,x/2);
}
Double_t chiF_ndf(double x,double p0,double p1){
    return chiF(p1*x,p0,p1);
}
Double_t chiF_ndf2(double x,double p0,double p1,double bin){
    return (chiF_ndf((x+bin/2),p0,27)-chiF_ndf((x-bin/2),p0,27));
}
Double_t chif_free(double x,double p0,double k,double p1){
    return p0*TMath::Exp(-x/(2*p1))*pow(x,(k/2)-1)/(pow(2*p1,k/2)*TMath::Gamma(k/2));
}
Double_t chiF_free(double x,double p0,double k,double p1){
    return p0*TMath::Gamma(k/2,x/(2*p1));
}
Double_t chiF_freefit(double x,double p0,double p1,double k,double bin){
    return (chiF_free((x+bin/2),p0,k,p1)-chiF_free((x-bin/2),p0,k,p1));
}
bool toge_hantei(vector<bool>ctoge, vector<bool>htoge, int bin, queue<int>& que){
    prep(i,bin,bin+dbin){
        if(ctoge[i] || htoge[i]){
            que.push(i);
            return true;
        }
    }
    return false;
}
Double_t gausf(double x,double p0,double mean,double sigma){
    return p0*TMath::Exp(-(x-mean)*(x-mean)/(2*sigma*sigma));
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
double v_conv(double f,double f0){
    double rtn=c*sqrt(1-((f0/f)*(f0/f)));
    return rtn;
}
double F_nu(double f,double f0){
    double rtn;
    double v=v_conv(f,f0);
    double p_kata=(v+vE)/v0;
    double m_kata=(v-vE)/v0;
    rtn = (vc/(2*sqrt(M_PI)*vE))*(exp(-(p_kata*p_kata))-exp(-(m_kata*m_kata)));
    rtn += 0.5*(erf(p_kata)+erf(m_kata));
    //rtn *= (c*f0*f0)/(f*f*f*sqrt(1-(f0/f)*(f0/f)));
    return rtn;
}
double F_sig1(double f,double f0,double P,double r){
    double rtn = P*(F_nu(f+r*dNu,f0)-F_nu(f-r*dNu,f0));
    return rtn;
}
double F_sig2(double f,double f0,double P,double r){
    if(f+dnu*r<=f0)return 0;
    else if(f+r*dNu>f0 && f-(1-r)*dNu<=f0){
        return P*(F_nu(f+r*dNu,f0)-F_nu(f0,f0));
    }
    else if(f-dnu*(1-r)>f0){
        return P*(F_nu(f+r*dNu,f0)-F_nu(f-(1-r)*dNu,f0));
    }
}
//パラメータ数を可変にしたときのフィッティング関数 できればフィット後のステータスも知りたい
/*
要件定義
1. coldとhotでアウトなものは無条件に抜く
2. mirrorに関しては4回分データを走査し、
*/
void baselinefit_2(){
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
    string savedir = "/Users/oginokyousuke/data/baseline_change/";
    string savedir2 = "/Users/oginokyousuke/data/rand_fit/";
    string savedird = "/Users/oginokyousuke/data/basic_data/";
    string savedirp = "/Users/oginokyousuke/data/peak_data/";

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
    //sigmaとchi2/ndfのエラーを調べる
    TGraphErrors* sigma_scale = new TGraphErrors;
    int ssbin = 0;
    axrange axss = {0,0.5,0,5,0,1};
    for(int i=1;i<25;i++){
        double dym = 0;
        int outnum = 0;
        TH1D* plus_ratio = new TH1D("plus_ratio","log10(plus_raito;d(chi2/ndf)/(chi2/ndf));Count",100,-20,0);
        TH1D* dist = new TH1D("dist","after-before;dist[d(chi2/ndf)];count",100,-0.000001,0.000001);
        
        TH1D* white_hzen = new TH1D("white_hzen","white_noise;dT[K];Count",100,-1,1);
        TH1D* white_hkou = new TH1D("white_hkou","white_noise;dT[K];Count",100,-1,1);
        TH1D* yscale_hist = new TH1D("yscale_hist","yscale;yscale[K];Count",100,0,15);
        vector<tuple<int,int,double>> revresult;
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
        TH1D* chi_hist = new TH1D("chi_hist","chi_hist;chi2/Ndf;Count",100,0,10);
        //要件定義: カイ二乗分布でフィットするのがなぜかうまく行かない理由を探る
        for(int j=0;j<4;j++){
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
            //setting.GraphErrors(pgraph1,ifmin,ifmax,0,100);
            //pgraph1 -> SetTitle("Prec;Freq[GHz];Prec[K]");
            //pgraph1 -> Draw("AP");
            
            TH1D* chihist = new TH1D("chihist","chihist;chi2/NDF;Coutn",100,0,10);
            double fm,fM;
            double chi2,chi2_d,chi22;
            int ndf,ndf_d,ndf2;
            double ysmax = -10;
            TH1D* ketah = new TH1D("hetah","ketah;PrecisionKeta;Count",21,0,20);
            //一回目のデータの処理
            TGraphErrors* wgraph1 = new TGraphErrors;
            for(int bin=sb;bin<fb;bin+=dbin){
                TH1D* whist = new TH1D("white","white_noise;dT[K];Count",100,-1,1);
                bool hantei = false;
                prep(k,bin,bin+dbin){
                    if(c_toge1[j][k] || h_toge1[j][k]){
                        hantei = true;
                        break;
                    }
                }
                if(hantei){
                    continue;
                }
                double yscale = 1000000;
                TGraphErrors* spgraph = new TGraphErrors;
                ft.make_scale(spgraph,pgraph1,bin-sb,yscale);
                
                TF1* f2 = new TF1("f2","[0]*(x-[1])*(x-[1])+[2]",0,1);
                TF1* gausfit = new TF1("gausfit","gaus",-1,1);
                
                //ft.rand_conv2(spgraph,f2,100,ketah,bin);
                double res2;
                ft.rand_fit(spgraph,f2,100,5,0,1,res2);
                chihist -> Fill(res2);
                rep(k,dbin){
                    double xValue = spgraph -> GetPointX(k);
                    double yValue = f2 -> Eval(xValue);
                    double yTrue = spgraph -> GetPointY(k);
                    whist -> Fill((yValue-yTrue)*yscale);
                    wgraph1 -> SetPoint(bin+k,Freq1[bin+k],(yValue-yTrue)*yscale);
                }
                st.Hist(whist);
                whist -> Draw();
                whist -> Fit(gausfit,"MQ","",-1,1);
                double sigma = gausfit -> GetParameter("Sigma");
                double esigma = gausfit -> GetParError(2);
                cout << sigma << " " << esigma << endl;
                for(int rb=bin;rb<bin+dbin;rb++)wgraph1 -> SetPointError(rb,0,sigma);
                
            }
            //pgraph1についてピークサーチしてその感度も出す
            /*
            どうフィットするか
            1. 注目する周波数を一つフィックスして前後数ビン(-10~10とか？)でフィット
            2. 丸々データがない時や統計量が少ない時(自由度でもカウントして)はスルー？
            3. フィットした時のピークの値、chi2/NDF、各データのシグマの値、なんでもいいので片っ端から保存？
            */
            
            int dameten = 0;
            TH1D* phist1 = new TH1D("phist","phist;P[];Count",100,-1,1);
            TH1D* chihist1 = new TH1D("chihist1","chihist1;Chi2/NDF;Count",100,0,10);
            for(int bin=sb+10;bin<fb-10;bin++){
                TF1* peakf = new TF1("peakf","F_sig2(x,[0],[1],0.5)",Freq1[bin-10],Freq1[bin+10]);
                peakf -> FixParameter(0,Freq1[bin]);
                peakf -> SetParameter(1,0.1);
                rep(k,5)wgraph1 -> Fit(peakf,"EMQ","",Freq1[bin-10],Freq1[bin+10]);
                int ndf = peakf -> GetNDF();
                double chi2 = peakf -> GetChisquare();
                if(ndf<10){
                    cout << Freq1[bin] << "GHz Can't Fit" << endl;
                    dameten++;
                    continue;
                }
                else{
                    double p = peakf -> GetParameter(1);
                    double perr = peakf -> GetParError(1);
                    cout << "p: " << p << " <=> perr: " << perr << endl;
                    if(abs(p)>1)phist1 -> Fill(1);
                    else phist1 -> Fill(p);
                    chihist1 -> Fill(chi2/ndf);
                }
                
            }
            st.Hist(phist1);
            st.Hist(chihist1);
            c1 -> SetLogy();
            filesystem::current_path(savedirp);
            phist1 -> Draw();
            phist1 -> Fit("gaus");
            string gname = "phist"+to_string(i)+"_"+to_string(j)+"_1.ps";
            c1 -> SaveAs(gname.c_str());
            chihist1 -> Draw();
            gname = "chihist"+to_string(i)+"_"+to_string(j)+"_1.ps";
            c1 -> SaveAs(gname.c_str());
            cout << "dameten : " << dameten << endl;



            //二回目のデータの処理　分離するのはいいとしてその基準と
            TGraphErrors* wgraph2 = new TGraphErrors;
            for(int bin=sb;bin<fb;bin+=dbin){
                TH1D* whist = new TH1D("white","white_noise;dT[K];Count",100,-1,1);
                //cout << Freq2[bin] << endl;
                bool hantei = false;
                prep(k,bin,bin+dbin){
                    if(c_toge2[j][k] || h_toge2[j][k]){
                        hantei = true;
                        break;
                    }
                }
                if(hantei)continue;
                double yscale = 100000;
                TGraphErrors* spgraph = new TGraphErrors;
                ft.make_scale(spgraph,pgraph2,bin-sb,yscale);
                TF1* f2 = new TF1("f2","[0]*(x-[1])*(x-[1])+[2]",0,1);
                TF1* gausfit = new TF1("gausfit","gaus",-1,1);
                double res2;
                ft.rand_fit(spgraph,f2,100,5,0,1,res2);
                chihist -> Fill(res2);
                rep(k,dbin){
                    double xValue = spgraph -> GetPointX(k);
                    double yValue = f2 -> Eval(xValue);
                    double yTrue = spgraph -> GetPointY(k);
                    whist -> Fill((yValue-yTrue)*yscale);
                    wgraph2 -> SetPoint(bin+k,Freq2[bin+k],(yValue-yTrue)*yscale);
                }
                st.Hist(whist);
                whist -> Draw();
                whist -> Fit(gausfit,"MQ","",-1,1);
                double sigma = gausfit -> GetParameter("Sigma");
                double esigma = gausfit -> GetParError(2);
                cout << sigma << " " << esigma << endl;
                for(int rb=bin;rb<bin+dbin;rb++)wgraph2 -> SetPointError(rb,0,sigma);
            }
            TH1D* phist2 = new TH1D("phist2","phist;P[];Count",100,-1,1);
            TH1D* chihist2 = new TH1D("chihist2","chihist1;Chi2/NDF;Count",100,0,10);
            for(int bin=sb+10;bin<fb-10;bin++){
                TF1* peakf = new TF1("peakf","F_sig2(x,[0],[1],0.5)",Freq2[bin-10],Freq2[bin+10]);
                peakf -> FixParameter(0,Freq1[bin]);
                peakf -> SetParameter(1,0.1);
                rep(k,5)wgraph2 -> Fit(peakf,"EMQ","",Freq2[bin-10],Freq2[bin+10]);
                int ndf = peakf -> GetNDF();
                double chi2 = peakf -> GetChisquare();
                if(ndf<10){
                    cout << Freq2[bin] << "GHz Can't Fit" << endl;
                    dameten++;
                    continue;
                }
                else{
                    double p = peakf -> GetParameter(1);
                    double perr = peakf -> GetParError(1);
                    cout << "p: " << p << " <=> perr: " << perr << endl;
                    if(abs(p)>1)phist2 -> Fill(1);
                    else phist2 -> Fill(p);
                    chihist2 -> Fill(chi2/ndf);
                }
                
            }
            st.Hist(phist2);
            st.Hist(chihist2);
            c1 -> SetLogy();
            filesystem::current_path(savedirp);
            phist2 -> Draw();
            phist2 -> Fit("gaus");
            gname = "phist"+to_string(i)+"_"+to_string(j)+"_2.ps";
            c1 -> SaveAs(gname.c_str());
            chihist2 -> Draw();
            gname = "chihist"+to_string(i)+"_"+to_string(j)+"_2.ps";
            c1 -> SaveAs(gname.c_str());
            cout << "dameten : " << dameten << endl;
        }
        
    }
    
}