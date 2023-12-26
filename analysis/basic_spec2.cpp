#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <format>
#include <tuple>
#include "../headers/fitter.h"
#include "../headers/mask_map.h"
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
string xfftname[4] = {"lsbo","lsbi","usbi","usbo"};
vector<double> errorv = {0.063,0.066,0.067,0.078,0.066,0.087,0.067,0.062,0.073,0.061,0.089,0.061,0.091,0.142,0.082,0.097,0.087,0.107,0.131,0.087,0.093,0.077,0.103,0.096};
//共通の物理量
vector<Color_t> gColor = {kBlue,kRed,kGreen,kMagenta};
const int falchan = 24061;
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
const double dnu=0.000076296;//ビン幅[GHz]
const double dNu=0.0000885;//周波数分解能[GHz]
const double kb=1.38*pow(10,-23);
const double df=88.5*pow(10,3);
const double Tc=76;
const double Th=297;
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
double F_sig2(double f,double f0,double P,double r){
    if(f+dnu*r<=f0)return 0;
    else if(f+r*dNu>f0 && f-(1-r)*dNu<=f0){
        return P*(F_nu(f+r*dNu,f0)-F_nu(f0,f0));
    }
    else if(f-dnu*(1-r)>f0){
        return P*(F_nu(f+r*dNu,f0)-F_nu(f-(1-r)*dNu,f0));
    }
    else return 0;
}
string basicdir = "/Users/oginokyousuke/data/basic_data/";
bool toge_hantei(vector<bool>ctoge, vector<bool>htoge, int bin, queue<int>& que){
    prep(i,bin,bin+dbin){
        if(ctoge[i] || htoge[i]){
            que.push(i);
            return true;
        }
    }
    return false;
}
int XFFT(int i){
    int xfft;
    auto result = find(lsbo.begin(),lsbo.end(),i);
    if(result==lsbo.end()){
        result = find(lsbi.begin(),lsbi.end(),i);
        if(result==lsbi.end()){
            result = find(usbo.begin(),usbo.end(),i);
            if(result==usbo.end())xfft = 2;
            else xfft = 3;
        }
        else xfft = 1;
    }
    else xfft = 0;
    return xfft;
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
void toge_scan2(bool (&hantei)[nbin],double input[nbin],double &sigma,double (&tvalue)[nbin]){
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
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
        if(abs(togevalue[bin])>5*sigma)hantei[bin] = true;
        tvalue[bin] = togevalue[bin]/sigma;
        num++;
    }
    
    //normhist -> Fit("gaus","Q");
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

//ミラーのみにしても点単位ではなく探索できなくなるmass単位でダメな点をカウントする
void basic_spec2(){
    //coldとmirrorの二階微分の相関？グラフを作成して有意に差があるか見てみるkernel空間
    Mask ms;
    vector<vector<int>> vec = ms.maskmap;
    bool hantei[4][nbin];
    rep(xfft,4)rep(bin,nbin)hantei[xfft][bin] = false;
    rep(xfft,4){
        for(auto v:vec[xfft])hantei[xfft][v] = true;
    }
    //今回の目玉となる、取り除きたいデータセット
    int togenum = 0;
    int hnomi = 0;
    int mnomi = 0;
    //種々のヘッダー関数を用意
    Setting st;
    st.dot_size = 0.8;
    st.markerstyle = 20;
    st.color = kBlue;
    st.lcolor = kBlue;
    Fitter ft;
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
    int cvec[nbin],hvec[nbin],mvec[nbin];
    TGraph* ctogecount = new TGraph;
    TGraph* htogecount = new TGraph;
    TGraph* mtogecount = new TGraph;
    rep(bin,nbin)cvec[bin]=0,hvec[bin]=0,mvec[bin]=0;
    //二次配列に、分光器ごとに一つでも悪いデータがあったらマークする(0:LSBO,1:LSBI,2:USBI,3:USBO)
    int hchannel[nbin][4],cchannel[nbin][4];
    rep(i,nbin)rep(j,4){
        cchannel[i][j] = 0;
        hchannel[i][j] = 0;
    }
    //ミリ波上で棄却されるエリアを記録する、統計がそもそも何点か+その中で捨てたデータが何点あるか、で見る
    //全周波数統一のマップ作る？(測定を超えた被り領域があるので)->めんどいから一旦そのままな
    bool maskchannel[4][nbin];
    rep(xfft,4)rep(bin,nbin)maskchannel[xfft][bin] = false;
    //最初に24セットでmaskchannel作る、それでもう一回for(i)を見たい帯域で回す
    TH1D* mexcess = new TH1D("mexcess","mexcess;Sigma;Count",100,5,20);
    TGraph* cmgsigma = new TGraph;//ミラーでposだった点のcoldでの値がどうなのか見る
    int cmnum = 0;
    //要件定義：map作り直し、きちんと切れているか確認すること
    bool cmask[4][nbin];
    rep(i,4)rep(j,nbin)cmask[i][j] = false;
    
    for(int i=9;i<10;i++){
        TLegend *legend = new TLegend(0.65, 0.65, 0.85, 0.85); 
        double fmin = 213.5+i*2;
        double fmax = 216.5+i*2;
        cout << i << " : " << fmin+0.5 << " ~ " << fmax-0.5 << endl;
        bool freqmap[nbin];
        int x = 0;
        rep(bin,nbin)freqmap[bin] = false;//一回でも探索したらtrueに変更する、最終的にfalseがあるかを計上？
        //band番号が奇数なら昇順、偶数なら降順のはず、昇順ならずらしたビンだけ引けばいい(どっちも同じじゃなくて？)
        int xfft = XFFT(i);
        string cdir=dir+"band"+to_string(i);
        filesystem::current_path(cdir);
        bool Ctoge[4][nbin],Htoge[4][nbin],Mtoge[4][nbin];
        rep(j,4){
            rep(k,nbin){
                // c_toge1[j][k] = false;
                // c_toge2[j][k] = false;
                // h_toge1[j][k] = false;
                // h_toge2[j][k] = false;
                // m_toge1[j][k] = false;
                // m_toge2[j][k] = false;
                Ctoge[j][k] = false;
                Htoge[j][k] = false;
                Mtoge[j][k] = false;
            }
        }
        for(int j=0;j<1;j++){
            //cout << j << endl;
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
            //if(xfft%2==1)cout << Freq1[2000+sbin[j]] << " : ";
            //else cout << Freq1[2000-sbin[j]] << " : "; 
            //積分値(平均値)
            double Cold[nbin],Hot[nbin],Mirror[nbin];
            rep(bin,nbin){
                Cold[bin] = (cold1[bin] + cold2[bin])/2;
                Hot[bin] = (hot1[bin] + hot2[bin])/2;
                Mirror[bin] = (mirror1[bin] + mirror2[bin])/2;
            }
            int outbin = sb;
            if(xfft%2==1)outbin += sbin[j];
            else outbin -= sbin[j];
            //cout << Freq1[outbin] << " ";
            //テスト関数(デルタ関数とシグナル)を用意してヒストグラム化,FFTで変換してその物性を確かめる
            axrange axtemp = {fmin,fmax,0,100,0,1,";Freq[GHz];Gain[a.u]"};
            axrange axgain = {fmin,fmax,0,pow(10,31),0,1,";Freq[GHz];Prec[K]"};
            axrange axdd = {fmin,fmax,-10,10,0,1,"ddmirror;bin;sigma"};
            axrange axraw = {fmin,fmax,0,pow(10,16),0,1,"Spectrum;Freq[GHz];Spectrum[a.u]"};
            axrange axd = {12439,12469,0,pow(10,16),0,1,"dcold;Bin;dcold[a.u]"};
            double Csigma,Hsigma,Msigma;
            double ddcold[nbin],ddhot[nbin],ddmirror[nbin];
            toge_scan2(Ctoge[j],Cold,Csigma,ddcold);
            toge_scan2(Htoge[j],Hot,Hsigma,ddhot);
            toge_scan2(Mtoge[j],Mirror,Msigma,ddmirror);
            rep(bin,nbin){
                if(Ctoge[j][bin])cmask[xfft][bin] = true;
            }
            
            TGraph* gcold = new TGraph;
            TGraph* ghot = new TGraph;
            TGraph* gmirror = new TGraph;
            TGraph* ggain = new TGraph;
            TGraph* gsys = new TGraph;
            TGraph* ddgcold = new TGraph;
            TGraph* ddghot = new TGraph;
            TGraph* ddgmirror = new TGraph;
            TGraph* dgcold = new TGraph;
            TGraph* dgmirror = new TGraph;
            TGraph* pgraph = new TGraph;
            double gain,psys;
            double ptemp;
            prep(bin,sb,fb){
                gcold -> SetPoint(bin-sb,Freq1[bin],Cold[bin]);
                dgcold -> SetPoint(bin,bin,(Cold[bin]-Cold[bin-1])/Cold[bin]);
                dgmirror -> SetPoint(bin,bin,(Mirror[bin]-Mirror[bin-1])/Mirror[bin]);
                ghot -> SetPoint(bin-sb,Freq1[bin],Hot[bin]);
                gmirror -> SetPoint(bin-sb,Freq1[bin],Mirror[bin]);
                ddgcold -> SetPoint(bin,Freq1[bin],ddcold[bin]);
                ddghot -> SetPoint(bin,bin,ddhot[bin]);
                ddgmirror -> SetPoint(bin,Freq1[bin],ddmirror[bin]);
                gain = (Hot[bin]-Cold[bin])/(2*kb*(Th-Tc)*df);
                psys = (Cold[bin]/gain)-2*kb*Tc*df;
                ptemp = (Mirror[bin]/gain-psys)/(2*kb*df);
                ggain -> SetPoint(bin-sb,Freq1[bin],gain);
                gsys -> SetPoint(bin-sb,Freq1[bin],psys);
                pgraph -> SetPoint(bin-sb,Freq1[bin],ptemp);
            }
            axrange axsys = {fmin,fmax,0,pow(10,-14),0,1,";Freq[GHz];Psys[W]"};
            st.Graph(pgraph,axtemp);
            st.Graph(gcold,axraw);
            st.Graph(ggain,axgain);
            st.Graph(gsys,axsys);
            c1 -> SetLogy();
            gcold -> SetLineColor(kBlue);
            ghot -> SetLineColor(kRed);
            gmirror -> SetLineColor(kGreen);
            gcold -> Draw("AL");
            ghot -> Draw("L");
            gmirror -> Draw("L");
            //ggain -> Draw("AL");
            //string lname = to_string(sbin[j]*76.2939453125*0.001)+"MHz";
            //legend->AddEntry(gcold, lname.c_str(), "l");
            //pgraph -> SetLineColor(gColor[j]);
            //ggain -> SetMarkerColor(kBlack);
            //ggain -> Draw("AP");
            /*if(j==0){
                pgraph -> Draw("AL");
            }
            else pgraph -> Draw("L");*/
            //cout << Freq1[18596]  << " : " << Freq1[22258] <<endl; 
            //怪しいチャンネルがどこで反応しているのか、何点反応しているのかなどを確かめる
            //cout << ddcold[falchan-1024] << " " << ddhot[falchan-1024] << endl;
            //st.Graph(gcold,axraw);
            
            //c1 -> SetLogy();
            
            TLegend *legend = new TLegend(0.7, 0.3, 0.85, 0.5);
            legend->AddEntry(gcold, "Cold", "l");
            legend->AddEntry(ghot, "Hot", "l");
            legend->AddEntry(gmirror, "Mirror", "l"); 
            legend->SetBorderSize(0); // 凡例のボーダーを非表示に
            legend->Draw();
            //axrange axg = {fmin,fmax,0,pow(10,31),0,1,"Gain;Freq[GHz];Gain[a.u]"};
            //axrange axsys = {fmin,fmax,0,pow(10,-15),0,1,"Psys;Freq[GHz];Psys[W]"};
            //st.Graph(ggain,axg);
            //st.Graph(gsys,axsys);
            //ggain -> Draw("AC");
            //ddgmirror -> Draw("AL");
            //ghot -> Draw("L");*/
            
        }
        //legend -> Draw();
    }
    
    
    
    
    /*rep(xfft,4){
        int xoutnum = 0;
        cout << xfftname[xfft] << ": " << endl;
        prep(bin,sb,fb-30){
            //if(cmask[xfft][bin])cout << bin << " " ;
            //4種類のLOでずらしていくとスキャンできない点が何点あるかを見てみる
            bool hantei = false;
            rep(dd,30){
                bool phantei = true;
                if(xfft%2==1){
                    rep(ite,4){
                        if(!cmask[xfft][bin+dd+sbin[ite]])phantei = false;
                    }
                }
                else{
                    rep(ite,4){
                        if(!cmask[xfft][bin+dd-sbin[ite]])phantei = false;
                    }
                }
                if(phantei){
                    hantei = true;
                    break;
                }
            }
            if(hantei)xoutnum++;
        }
        cout << xoutnum << endl;
    }*/
}