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
    
    hist -> Draw();
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
    for(int i=1;i<2;i++){
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
        for(int j=3;j<4;j++){
            cout << j << endl;
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
            int outbin = sb;
            if(xfft%2==1)outbin += sbin[j];
            else outbin -= sbin[j];
            //cout << Freq1[outbin] << " ";
            //テスト関数(デルタ関数とシグナル)を用意してヒストグラム化,FFTで変換してその物性を確かめる
            axrange axdd = {fmin,fmax,-10,10,0,1,"ddmirror;bin;sigma"};
            axrange axraw = {fmin,fmax,0,pow(10,16),0,1,"Cold;Bin;Cold[a.u]"};
            axrange axd = {12439,12469,0,pow(10,16),0,1,"dcold;Bin;dcold[a.u]"};
            double Csigma,Hsigma,Msigma;
            double ddcold[nbin],ddhot[nbin],ddmirror[nbin];
            toge_scan2(Ctoge[j],Cold,Csigma,ddcold);
            toge_scan2(Htoge[j],Hot,Hsigma,ddhot);
            toge_scan2(Mtoge[j],Mirror,Msigma,ddmirror);
            prep(bin,sb,fb){
                if(abs(ddmirror[bin])>5){
                    bool rhantei = false;
                    prep(ad,-15,15){
                        if(hantei[xfft][bin+ad]){
                            rhantei = true;
                            break;
                        }
                    }
                }
            }
            
            double jfmin = min(Freq1[0],Freq1[nbin-1]);
            double jfmax = max(Freq1[0],Freq1[nbin-1]);
            TH1* thist = new TH1D("thist",";Freq[GHz];Spec",nbin,jfmin,jfmax);
            TH1* onhist = new TH1D("onhist",";Freq[GHz];Spec",nbin,0,nbin);
            TH1* delhist = new TH1D("delhist","test;Freq[GHz];Spec",nbin,jfmin,jfmax);
            TH1* fhist = new TH1D("fhist","MAG;1/Freq[ns];",nbin,0,nbin/2.5);//棘抜く前のFFT
            TH1* tfhist = new TH1D("tfhist","PH;1/Freq[ns];",nbin,0,nbin/2.5);//棘抜いた後のFFT
            TH1* sighist = new TH1D("sighist","sig;Freq[GHz];Spec",30,Freq1[0],Freq1[30]);
            TF1* hamf = new TF1("hamf","");
            TH1* sintest = new TH1D("sintest","test;test;test",30,0,30);
            TH1* sinfft = new TH1D("sinfft","test;test;test",15,0,15);
            TH1* sigtest = nullptr;
            
            rep(bin,30)sintest -> SetBinContent(bin,TMath::Cos(2*TMath::Pi()*bin/5));
            st.Hist(sintest);
            sintest -> Draw();
            sinfft = sintest -> FFT(sinfft,"MAG");
            st.Hist(sinfft);
            sinfft -> Draw();
            //TVirtualFFT::SetTransform(nullptr);
            /*TH1* parhist = new TH1D("parhsit","partition;Freq[GHz];Spec",30,Freq1[2940],Freq1[2970]);
            rep(k,30)parhist -> SetBinContent(k,Cold[25940+k]/pow(10,13));
            st.Hist(parhist);
            parhist -> Draw();
            fhist = parhist -> FFT(fhist,"MAG");
            st.Hist(fhist);*/
            //fhist -> Draw();
            //30binの
            /*rep(bin,30){
                thist -> SetBinContent(bin,Cold[bin]/pow(10,13));
                sighist -> SetBinContent(bin,F_sig2(Freq1[bin],Freq1[1],10,0.5));
            }
            st.Hist(sighist);
            sighist -> Draw();
            sigtest = sighist -> FFT(sigtest,"MAG");
            st.Hist(sigtest);
            sigtest -> Draw();

            //filter作り、カットオフ有無で二通り
            TH1* filhist = new TH1D("filhist","wiener filter[not cutoff];1/Freq[ns];",nbin,0,nbin/2.5);
            TH1* filhist2 = new TH1D("filhist2","wiener filter[cutoff];1/Freq[ns];",nbin,0,nbin/2.5);
            TH1* tfilhist = new TH1D("tfilhist","after filtered[not cutoff];",nbin,0,nbin/2.5);
            TH1* tfilhist2 = new TH1D("tfilhist2","after filtered[cutoff];",nbin,0,nbin/2.5);

            Double_t *re_full = new Double_t[nbin];//シグナルのFFTデータ
            Double_t *im_full = new Double_t[nbin];//上の虚部
            Double_t *tre_full = new Double_t[nbin];//元データのFFTデータ
            Double_t *tim_full = new Double_t[nbin];//上の虚部
            //シグナルと生データをそれぞれFFTにかける
            fhist = thist -> FFT(fhist,"MAG");
            TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();//直近でFFTしたデータのあれこれを引っ張り出すポインタ
            fft -> GetPointsComplex(re_full,im_full);
            double cre = re_full[100];
            double cim = im_full[100];
            tfhist = thist -> FFT(tfhist,"MAG");
            TVirtualFFT *tfft = TVirtualFFT::GetCurrentTransform();
            tfft -> GetPointsComplex(tre_full,tim_full);
        
            rep(bin,nbin){
                double deno = (re_full[bin]*re_full[bin]+im_full[bin]*im_full[bin]);
                double deno2 = (re_full[bin]*re_full[bin]+im_full[bin]*im_full[bin])+(cre*cre+cim*cim);
                double re = re_full[bin]/deno;
                double im = -im_full[bin]/deno;
                double re2 = re_full[bin]/deno2;
                double im2 = -im_full[bin]/deno2;
                double tre = tre_full[bin]*re-tim_full[bin]*im;
                double tim = tre_full[bin]*im+tim_full[bin]*re;
                double tre2 = tre_full[bin]*re2-tim_full[bin]*im2;
                double tim2 = tre_full[bin]*im2+tim_full[bin]*re2;
                tfilhist -> SetBinContent(bin,sqrt(tre*tre+tim*tim));
                tfilhist2 -> SetBinContent(bin,sqrt(tre2*tre2+tim2*tim2));
                filhist -> SetBinContent(bin,sqrt(re*re+im*im));
                filhist2 -> SetBinContent(bin,sqrt(re2*re2+im2*im2));
            }
            st.Hist(filhist2);
            filhist2 -> Draw();
            st.Hist(tfilhist);
            tfilhist -> Draw();

            /*;//各点のDCとナイキスト周波数を取得
            st.Hist(fhist);
            fhist -> Draw();
            
            //これを踏まえた時にフィルタがどう機能しているのかを確認する
            fhist2 = sighist -> FFT(fhist2,"MAG");
            TVirtualFFT *fft2 = TVirtualFFT::GetCurrentTransform();//直近でFFTしたデータのあれこれを引っ張り出すポインタ
            fft2 -> GetPointsComplex(mre_full,mim_full);//各点のDCとナイキスト周波数を取得
            double mcr = mre_full[100];
            double mci = mim_full[100];
            TH1* reshist = new TH1D("reshist","res;1/Freq[ns];Spec",nbin,0,nbin/2.5);
            rep(bin,nbin){
                double tr = tre_full[bin];
                double ti = tim_full[bin];
                double mr = mre_full[bin];
                double mi = mim_full[bin];
                double K = (mcr*mcr+mci*mci)+(mr*mr+mi*mi);
                tre_full[bin] = (tr*mr+ti*mi)/K;
                tim_full[bin] = (ti*mr-tr*mi)/K;
                reshist -> SetBinContent(bin,sqrt(tre_full[bin]*tre_full[bin]+tim_full[bin]*tim_full[bin]));
            }
            st.Hist(reshist);
            reshist -> Draw();
            Int_t Nbin = 32767;
            // TH1* backhist = new TH1D("backhist","back;Freq[GHz];Spec",nbin,jfmin,jfmax);
            // TVirtualFFT *fft_back = TVirtualFFT::FFT(1,&Nbin,"C2R M K");
            // fft_back->SetPointsComplex(tre_full,tim_full);
            // fft_back -> Transform();
            // backhist = TH1::TransformHisto(fft_back,backhist,"Re");
            // st.Hist(backhist);
            // backhist -> Draw();
            // 
            // 
            // fhist2 = onhist -> FFT(fhist2,"MAG");
            // TVirtualFFT *fft2 = TVirtualFFT::GetCurrentTransform();
            // fft2 -> GetPointsComplex(mre_full,mim_full);
            // rep(bin,nbin){
            //     double tr = tre_full[bin];
            //     double ti = tim_full[bin];
            //     double mr = mre_full[bin];
            //     double mi = mim_full[bin];
            //     double deno = mr*mr+mi*mi;
            //     tre_full[bin] = (tr*mr+ti*mi)/deno;
            //     tim_full[bin] = (ti*mr-tr*mi)/deno;
            // }
            // 
            // TVirtualFFT *fft_back = TVirtualFFT::FFT(1,&Nbin,"C2R M K");
            // fft_back->SetPointsComplex(tre_full,tim_full);
            // fft_back->Transform();
            // TH1 *hb = new TH1D("hb","back;Freq[GHz];Spec",nbin,jfmin,jfmax);
            // //Let's look at the output
            // hb = TH1::TransformHisto(fft_back,hb,"Re");
            // hb->SetTitle("The backward transform result");
            // st.Hist(hb);
            // hb->Draw();
            
            /*for(int bin=nbin-1;bin>=nbin;bin++){
                re_full[bin] = re_full[nbin-1-bin];
                im_full[bin] = im_full[nbin-1-bin];
            }
            Int_t Nbin = 32767;
            TVirtualFFT *fft_back = TVirtualFFT::FFT(1,&Nbin,"C2R M K");
            fft_back->SetPointsComplex(re_full,im_full);
            fft_back->Transform();
            TH1 *hb = nullptr;
            //Let's look at the output
            hb = TH1::TransformHisto(fft_back,hb,"Re");
            hb->SetTitle("The backward transform result");
            st.Hist(hb);
            hb->Draw();
            */
            /*TGraph* gcold = new TGraph;
            TGraph* ghot = new TGraph;
            TGraph* gmirror = new TGraph;
            TGraph* ggain = new TGraph;
            TGraph* gsys = new TGraph;
            TGraph* ddgcold = new TGraph;
            TGraph* ddghot = new TGraph;
            TGraph* ddgmirror = new TGraph;
            TGraph* dgcold = new TGraph;
            TGraph* dgmirror = new TGraph;
            double gain,psys;
            rep(bin,nbin){
                gcold -> SetPoint(bin,Freq1[bin],Cold[bin]);
                dgcold -> SetPoint(bin,bin,(Cold[bin]-Cold[bin-1])/Cold[bin]);
                dgmirror -> SetPoint(bin,bin,(Mirror[bin]-Mirror[bin-1])/Mirror[bin]);
                ghot -> SetPoint(bin,Freq1[bin],Hot[bin]);
                gmirror -> SetPoint(bin,Freq1[bin],Mirror[bin]);
                ddgcold -> SetPoint(bin,Freq1[bin],ddcold[bin]);
                ddghot -> SetPoint(bin,bin,ddhot[bin]);
                ddgmirror -> SetPoint(bin,Freq1[bin],ddmirror[bin]);
                gain = (Hot[bin]-Cold[bin])/(2*kb*(Th-Tc)*df);
                psys = (Cold[bin]/gain)-2*kb*Tc*df;
                ggain -> SetPoint(bin,Freq1[bin],gain);
                gsys -> SetPoint(bin,Freq1[bin],psys);
            }
            //cout << Freq1[18596]  << " : " << Freq1[22258] <<endl; 
            //怪しいチャンネルがどこで反応しているのか、何点反応しているのかなどを確かめる
            //cout << ddcold[falchan-1024] << " " << ddhot[falchan-1024] << endl;
            st.Graph(gcold,axraw);
            st.Graph(gmirror,axraw);
            st.Graph(ghot,axraw);
            st.Graph(ddgcold,axdd);
            st.Graph(ddgmirror,axdd);
            st.Graph(dgcold,axd);
            st.Graph(dgmirror,axd);
            string ctitle = "Cold"+to_string(i)+"_"+to_string(j)+";Freq[GHz];Cold[a.u]";
            string htitle = "Hot"+to_string(i)+"_"+to_string(j)+";Freq[GHz];Hot[a.u]";
            string mtitle = "Mirror"+to_string(i)+"_"+to_string(j)+";Freq[GHz];Mirror[a.u]";
            
            //c1 -> SetLogy();
            ghot -> SetLineColor(kRed);
            gcold -> SetLineColor(kBlue);
            gmirror -> SetLineColor(kGreen);
            ddgmirror -> SetLineColor(kGreen);
            //gcold -> SetTitle(ctitle.c_str());
            gcold -> SetTitle("Spectrum;Freq[GHz];Spectrum[a.u]");
            ghot -> SetTitle(htitle.c_str());
            gmirror -> SetTitle(mtitle.c_str());
            ctitle = "Cold"+to_string(i)+"_"+to_string(j)+".ps";
            htitle = "Hot"+to_string(i)+"_"+to_string(j)+".ps"; 
            mtitle = "Mirror"+to_string(i)+"_"+to_string(j)+".ps";
            gcold -> Draw("AL");
            ghot -> Draw("L"); 
            gmirror -> Draw("L");
            TLegend *legend = new TLegend(0.7, 0.5, 0.85, 0.7);
            legend->AddEntry(gcold, "Cold", "l");
            legend->AddEntry(ghot, "Hot", "l");
            legend->AddEntry(gmirror, "Mirror", "l"); 
            legend->SetBorderSize(0); // 凡例のボーダーを非表示に
            //legend->Draw();
            axrange axg = {fmin,fmax,0,pow(10,31),0,1,"Gain;Freq[GHz];Gain[a.u]"};
            axrange axsys = {fmin,fmax,0,pow(10,-15),0,1,"Psys;Freq[GHz];Psys[W]"};
            st.Graph(ggain,axg);
            st.Graph(gsys,axsys);
            //ggain -> Draw("AC");
            //ddgmirror -> Draw("AL");
            //ghot -> Draw("L");*/
            
        }
    }
    /*cout << cmnum << endl;
    axrange axcm = {0,20,0,20,0,1,";Msigma;Csigma"};
    st.Graph(cmgsigma,axcm);
    cmgsigma -> Draw("AP");

    int mnomi2 = 0;
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
    /*TF1* corf = new TF1("corf","x",0.0003,0.0004);
    TF1* corfu = new TF1("corfu","1.03*x",0.0003,0.0004);
    TF1* corfl = new TF1("corfl","0.97*x",0.0003,0.0004);
    axrange axsig = {0.00033,0.00037,0.00033,0.00037,0,1,"cold-hot;Csigma;Hsigma"};
    st.Graph(cmgsigma,axsig);
    cmgsigma -> Draw("AP");
    corf -> Draw("same");
    corfu -> Draw("same");
    corfl -> Draw("same");*/
    //cout << "after : " <<mnomi2 << endl;
}