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
const double dnu=0.000076296;//ビン幅
const double dNu=0.0000885;//周波数分解能
const double kb=1.38*pow(10,-23);
const double df=88.5*pow(10,3);
const double Tc=76;
const double Th=297;
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
    TH1D* hist = new TH1D(hname.c_str(),"hist;ddinput;Count",100,-0.005,0.005);
    TF1* f = new TF1("gaus","gaus",-0.005,0.005);
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
    hist -> Fit(f,"QEN","",-0.005,0.005);
    //st.Hist(hist);
    sigma = f -> GetParameter("Sigma");
    int num = 0;
    for(int bin=sb;bin<fb;bin++){
        if(abs(togevalue[bin])>5*sigma)hantei[bin] = true;
        tvalue[bin] = togevalue[bin]/sigma;
        num++;
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

//ミラーのみにしても点単位ではなく探索できなくなるmass単位でダメな点をカウントする
void basic_spec2(){
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
    bool maskchannel[nbin][4];
    rep(bin,nbin)rep(xfft,4)maskchannel[bin][xfft] = false;
    //最初に24セットでmaskchannel作る、それでもう一回for(i)を見たい帯域で回す
    TH1D* mexcess = new TH1D("mexcess","mexcess;Sigma;Count",100,5,20);
    //for(auto)
    //要件定義(周波数のスキャンが潰れたところを精査したい)チャンネルを全パターンboolで持っておくしかないンゴ！
    for(int i=1;i<25;i++){
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
            
            //積分値(平均値)
            double Cold[nbin],Hot[nbin],Mirror[nbin];
            rep(bin,nbin){
                Cold[bin] = (cold1[bin] + cold2[bin])/2;
                Hot[bin] = (hot1[bin] + hot2[bin])/2;
                Mirror[bin] = (mirror1[bin] + mirror2[bin])/2;
            }
            
            axrange axdd = {17390,17430,-10,10,0,1,"ddcold;bin;sigma"};
            axrange axraw = {17390,17430,0,pow(10,16),0,1,"Cold;Bin;Cold[a.u]"};
            double Csigma,Hsigma,Msigma;
            double ddcold[nbin],ddhot[nbin],ddmirror[nbin];
            toge_scan2(Ctoge[j],Cold,Csigma,ddcold);
            toge_scan2(Htoge[j],Hot,Hsigma,ddhot);
            toge_scan2(Mtoge[j],Mirror,Msigma,ddmirror);
            //ここで作ったboolの配列を30binごとに区切って区間の中に棘を入れないようにする→何データ死ぬかを確認
            prep(bin,sb,fb){
                bool chantei = false;
                bool hhantei = false;
                bool mhantei = false;
                bool maskhantei = false;
                rep(ad,dbin){
                    if(hantei[xfft][bin+ad]){
                        maskhantei = true;
                        break;
                    }
                    if(Ctoge[j][bin+ad])chantei = true;
                    if(Htoge[j][bin+ad])hhantei = true;
                    if(Mtoge[j][bin+ad])mhantei = true;
                }
                if(!maskhantei){
                    if(i%2==1)freqmap[bin+11-sbin[j]] = true;
                    else freqmap[bin+11+sbin[j]] = true;
                }
                /*if(maskhantei)continue;
                if(chantei || hhantei)togenum++;
                else if(!chantei && hhantei)hnomi++;
                else if(!chantei && !hhantei && mhantei){
                    //cout << j << " " << bin << endl;
                    mnomi++;
                }*/
            }
            
            TGraph* gcold = new TGraph;
            TGraph* ghot = new TGraph;
            TGraph* gmirror = new TGraph;
            TGraph* ddgcold = new TGraph;
            TGraph* ddghot = new TGraph;
            TGraph* ddgmirror = new TGraph;
            rep(bin,nbin){
                gcold -> SetPoint(bin,bin,Cold[bin]);
                ghot -> SetPoint(bin,bin,Hot[bin]);
                gmirror -> SetPoint(bin,bin,Mirror[bin]);
                ddgcold -> SetPoint(bin,bin,ddcold[bin]);
                ddghot -> SetPoint(bin,bin,ddhot[bin]);
                ddgmirror -> SetPoint(bin,bin,ddmirror[bin]);
            }

            //怪しいチャンネルがどこで反応しているのか、何点反応しているのかなどを確かめる
            //cout << ddcold[falchan-1024] << " " << ddhot[falchan-1024] << endl;
            st.Graph(gcold,axraw);
            st.Graph(gmirror,axraw);
            st.Graph(ddgcold,axdd);
            st.Graph(ddgmirror,axdd);
            //ddgcold -> SetLineColor(gColor[j]);
            gmirror -> SetLineColor(kGreen);
            ddgcold -> SetLineColor(kBlue);
            ddgmirror -> SetLineColor(kGreen);
            ghot -> SetLineColor(kRed);
            gcold -> Draw("AL");
            //ddgcold -> Draw("AL");
            //ghot -> Draw("L");
            
        }
        
        prep(bin,sb+11,fb){
            if(!freqmap[bin]){
                //cout << bin << endl;
                x++;
            }
        }
        cout << "x : " << x << endl;
    
    }
    
    /*cout << "before : " << mnomi << endl;
    //xfft4パターン,lsbo,usbiは昇順、lsbi,usboは降順
    ofstream writing_file;
    prep(xfft,0,4){
        //cout << xfftname[xfft] << " : ";
        int outnum = 0;
        prep(bin,sb,fb){
            if(maskchannel[bin][xfft]){
                //cout << bin  << " ," ;
                outnum++;
            }
        }
        //探索範囲をどこからどこまでに設定するべきかは非常に悩ましい問題である、がとりあえず2GHz基準で
        /*prep(bin,sb,fb){
            if(xfft%2==1){
                if(maskchannel[bin][xfft] && maskchannel[bin-512][xfft] && maskchannel[bin+512][xfft] && maskchannel[bin+1024][xfft]){
                    cout << bin << " ";
                }
            }
            else{
                if(maskchannel[bin][xfft] && maskchannel[bin+512][xfft] && maskchannel[bin-512][xfft] && maskchannel[bin-1024][xfft]){
                    cout << bin << " ";
                }
            }
        }
        //cout << outnum << endl;
    }
    /*int mnomi2 = 0;
    //フィルタリングした上で001のパターンがカウントされていないかどうかを確認する
    prep(i,2,3){
        int xfft = XFFT(i);
        string cdir=dir+"band"+to_string(i);
        filesystem::current_path(cdir);
        bool Ctoge[4][nbin],Htoge[4][nbin],Mtoge[4][nbin];
        rep(j,4){
            rep(k,nbin){
                Ctoge[j][k] = false;
                Htoge[j][k] = false;
                Mtoge[j][k] = false;
            }
        }
        prep(j,0,1){
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
            axrange axrawc = {2620,2630,0,pow(10,16),0,1,"Cold;Bin;Cold[a.u]"};
            axrange axrawh = {2620,2630,0,pow(10,16),0,1,"Hot;Bin;Hot[a.u]"};
            axrange axrawm = {2620,2630,0,pow(10,16),0,1,"Mirror;Bin;Mirror[a.u]"};
            axrange axdd = {2620,2630,-10,10,0,1,"ddmirror;Bin;ddmirror"};
            TGraph* gcold = new TGraph;
            TGraph* ghot = new TGraph;
            TGraph* gmirror = new TGraph;
            double Cold[nbin],Hot[nbin],Mirror[nbin];
            rep(bin,nbin){
                Cold[bin] = (cold1[bin] + cold2[bin])/2;
                Hot[bin] = (hot1[bin] + hot2[bin])/2;
                Mirror[bin] = (mirror1[bin] + mirror2[bin])/2;
                gcold -> SetPoint(bin,bin,Cold[bin]);
                ghot -> SetPoint(bin,bin,Hot[bin]);
                gmirror -> SetPoint(bin,bin,Mirror[bin]);
            }
            st.Graph(gcold,axrawc);
            st.Graph(ghot,axrawh);
            st.Graph(gmirror,axrawm);
            
            gcold -> SetMarkerColor(kBlue);
            ghot -> SetMarkerColor(kRed);
            gmirror -> SetMarkerColor(kGreen);
            gmirror -> Draw("AP");

            TGraph* ddcold = new TGraph;
            TGraph* ddhot = new TGraph;
            TGraph* ddmirror = new TGraph;
            double Csigma,Hsigma,Msigma;
            toge_scan2(Ctoge[j],Cold,Csigma,ddcold);
            toge_scan2(Htoge[j],Hot,Hsigma,ddhot);
            toge_scan2(Mtoge[j],Mirror,Msigma,ddmirror);
            st.Graph(ddcold,axdd);
            st.Graph(ddhot,axdd);
            st.Graph(ddmirror,axdd);
            ddcold -> SetLineColor(kBlue);
            ddhot -> SetLineColor(kRed);
            ddmirror -> SetLineColor(kGreen);
            ddmirror -> Draw("AL");
            prep(bin,sb,fb){
                if(maskchannel[bin][xfft])continue;
                bool hantei = false;//いくつか棘の判定を行う部署を
                if(!Ctoge[j][bin] && !Htoge[j][bin] && Mtoge[j][bin]){
                    //例えその点でnegativeでも周辺でpositiveであれば見逃す
                    if(Ctoge[j][bin-1] || Htoge[j][bin-1] || Ctoge[j][bin+1] || Htoge[j][bin+1])continue;
                    cout << xfftname[xfft] <<  " i:" << i << " <=> j: " << j << " <=> bin: " << bin << endl;
                    mnomi2++;
                }
            }
        }
    }*/
    //cout << "after : " <<mnomi2 << endl;
}