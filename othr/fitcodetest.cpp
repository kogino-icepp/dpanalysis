//スプリアス検出コード並びにピークサーチのコードが実データ内でしっかり機能しているかどうかを、テストシグナルを用いて検証する
#include <iostream>
#include "../headers/fitter.h"
#include <random>
#include "../headers/mask_map.h"
using namespace std;
const int nbin=32767;
const int sb = 2585;//切り出してくるビンの最初
const int cb = 15725;//探索範囲の中点,ここを境にデータの振幅が変わっているデータがある
const int fb = 28865;//切り出してくるビンの最後
const double c= 3.0*pow(10,8);
const double v0=220000.0;
const double vE=v0;
const double vc=v0;
const double dnu=0.000076296;//ビン幅[GHz]
const double dNu=0.0000885;//周波数分解能
const double kb=1.38*pow(10,-23);
const double df=88.5*pow(10,3);
const double Tc=76;
const double Th=297;
axrange axscale = {0,1,0,1,0,1,"test_data;xscale;yscale"};
vector<int> sbin = {0,512,-512,-1024};
const vector<int> lsbo = {1,3,5,13,15,17};
const vector<int> lsbi = {2,4,6,14,16,18};
const vector<int> usbi = {7,9,11,19,21,23};
const vector<int> usbo = {8,10,12,20,22,24};
string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
Fitter ft;

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
double v_conv(double f,double f0){
    double rtn=c*sqrt(1-((f0/f)*(f0/f)));
    return rtn;
}
double PtoChi(double P){
    return 4.5*pow(10,-14)*sqrt(P*pow(10,23))*sqrt(1/0.385);
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
double F_sigscale(double x,double P,double r,double fmin){
    double bin = x/0.0344827586;
    
    double f = fmin+bin*dnu;//ここをxとfminの関数に変える
    double f0 = fmin+10*dnu;//ここは自明ではあるがx0からf0へ
    return F_sig2(f,f0,P,r);
}
void fitcodetest(){
    Mask ms;
    vector<vector<int>> mskmap = ms.maskmap;
    bool mskhantei[4][nbin];
    for(int j=0;j<4;j++){
        for(int bin=0;bin<nbin;bin++){mskhantei[j][bin] = false;}
    }
    for(int j=0;j<4;j++){//mskmapに記録されているものだけtrueに変更
        for(int bin:mskmap[j])mskhantei[j][bin] = true;
    }
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
    st.markerstyle = 20;
    st.color = kBlue;
    st.lcolor = kBlue;
    //まずは棘検出がどの程度できるか
    TGraph* sigmagraph = new TGraph;
    double ratio[10] = {0.692,0.97,0.995,0.997,0.999,1,1,1,1,1};
    int itenum = 300;
    /*
    考えるべき事項
    1. 載せるシグナルの積分値が大体どのくらいであるのが良いのか(感度に直すとどのくらいといえば良いのだろうか)
    2. 今回はキャリブレーションデータに載せるのが対象なのでまずはそこの取り出しからやるべき
    */
    //TGraphErrors* resgraph = new TGraphErrors;
    TGraph* resgraph = new TGraph;
    for(int pg=1;pg<15;pg+=2){
        double pgiven = pg*0.1;
        int rnum = 0;
        vector<double> reslist;
        rep(ite,itenum){
            int i = randband(mt);
            int j0 = randj(mt);
            int mbin = randbin(mt);//質量に対応するbin
            int j = j0/2;
            int p = 1+j0%2;
            
            int xfft = XFFT((i));
            bool hantei = false;
            prep(bin,-30,31){
                if(mskhantei[xfft][mbin+bin]){
                    hantei = true;
                    break;
                }
            }
            if(hantei)continue;
            filesystem::path path=filesystem::current_path();
            string cdir=dir+"band"+to_string(i);
            filesystem::current_path(cdir);
            //とりあえずcoldだけ取り出せたらラッキー
            double Freq[nbin],cold[nbin],hot[nbin],mirror[nbin];
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
                        else if(num==(2-p) && body==1){
                            for(int bin=0;bin<nbin;bin++){
                                tree->GetEntry(bin);
                                hot[bin]=power;
                            }
                        }
                        else if(num==(2-p) && body==2){
                            for(int bin=0;bin<nbin;bin++){
                                tree->GetEntry(bin);
                                mirror[bin]=power;
                            }
                        }
                    }
                    
                }
            }
            double caldata[nbin];
            double gain,psys,ptemp;
            TGraph* pgraph = new TGraph;
            TGraphErrors* spgraph = new TGraphErrors;
            
            int sbin,fbin;
            if(i%2==1){
                sbin = mbin-10;
                fbin = mbin+19;
            }
            else{
                sbin = mbin-19;
                fbin = mbin+10;
            }
            
            prep(bin,sbin,fbin+1){
                gain = (hot[bin]-cold[bin])/(2*kb*(Th-Tc)*df);
                psys = (cold[bin]/gain)-2*kb*Tc*df;
                ptemp = (mirror[bin]/gain-psys)/(2*kb*df);
                caldata[bin] = ptemp+F_sig2(Freq[bin],Freq[mbin],pgiven,0.5);
                pgraph -> SetPoint(bin-sbin,Freq[bin],caldata[bin]);
            }
            double yscale;
            ft.make_scale(spgraph,pgraph,0,yscale);
            TGraphErrors* spgraphk = new TGraphErrors;
            rep(spbin,dbin){
                if(spbin<10 || spbin>=20){
                    double x = spgraph -> GetPointX(spbin);
                    double y = spgraph -> GetPointY(spbin);
                    double ye = spgraph -> GetErrorY(spbin);
                    spgraphk -> SetPoint(spbin,x,y);
                    spgraphk -> SetPointError(spbin,0,ye);
                }
            }
            double res1;
            TF1* f1 = new TF1("f1","[0]*(x-[1])*(x-[1])+[2]",0,1);
            TF1* scalepeak = new TF1("scalepeak","[0]*(x-[1])*(x-[1])+[2]+F_sigscale(x,[3],0.5,[4])",0,1);
            ft.allfit2(spgraphk,f1,5,res1);
            double a = f1 -> GetParameter(0);
            double b = f1 -> GetParameter(1);
            double c = f1 -> GetParameter(2);
            double terr = 0;
            rep(spbin,dbin){
                if(spbin<10 || spbin>=20){
                    double x = spgraphk -> GetPointX(spbin);
                    double y = spgraphk -> GetPointY(spbin);
                    terr += pow(y-a*(x-b)*(x-b)-c,2);
                }
            }
            terr = sqrt(terr/17);
            rep(spbin,30){
                spgraph -> SetPointError(spbin,0,terr);
                spgraphk -> SetPointError(spbin,0,terr);
            }
            double fmin = min(Freq[sbin],Freq[fbin]);
            double fmax = max(Freq[sbin],Freq[fbin]);
            axrange axite = {fmin,fmax,0,100,0,1,";Freq[GHz];Prec[K]"};
            
            scalepeak -> SetParameter(0,a);
            scalepeak -> SetParameter(1,b);
            scalepeak -> SetParameter(2,c);
            scalepeak -> FixParameter(4,fmin);
            scalepeak -> SetParameter(3,0.1);
            rep(ite,20)spgraph -> Fit(scalepeak,"IQ0","",0,1);
            rep(ite,20)spgraph -> Fit(scalepeak,"IMQ0","",0,1);
            rep(ite,1)spgraph -> Fit(scalepeak,"IEQ","",0,1);
            double res = scalepeak -> GetParameter(3);
            st.GraphErrors(spgraphk,axscale);
            spgraphk -> SetMarkerColor(kRed);
            spgraphk -> SetLineColor(kRed);
            //cout << res*yscale << endl;
            st.GraphErrors(spgraph,axscale);
            spgraph -> Draw("AP");
            spgraphk -> Draw("P");
            scalepeak -> Draw("same");
            //st.SetStatInfo(c1,spgraph,scalepeak);
            // TPaveText* paveText = new TPaveText(0.15, 0.75, 0.35, 0.95, "brNDC");
            // paveText -> SetName("fitResults");
            // paveText -> SetBorderSize(1);
            // paveText -> SetFillColor(0);
            // paveText -> SetTextAlign(12);
            // paveText -> SetTextFont(42);
            // c1 -> Update();  // Make sure the canvas is up-to-date for gPad to work
            // c1 -> Modified();
            // paveText -> AddText("Fit Results:");
            // paveText -> AddText(Form(" Chi: %.4f", scalepeak->GetChisquare()));
            // paveText -> AddText(Form(" NDF: %d", scalepeak->GetNDF()));
            // paveText -> AddText(Form(" Error: %.4f",spgraph -> GetErrorY(0)));
            // paveText->Draw();
        
            reslist.push_back(res*yscale);
            //cout << res*yscale << endl;
            rnum++;

            //cout << Freq[mbin] << endl;
            //これもしかして面倒臭いけどまたスケール→フィット→リスケール？？
        }
        pair<double,double> fin = MeanError(reslist);
        resgraph -> SetPoint(pg-1,pgiven,fin.first/fin.second);
        //resgraph -> SetPointError(pg-1,0,fin.second);
        //cout<< fin.first << " +- " << fin.second << endl;
        /*rep(i,10)sigmagraph -> SetPoint(i,i+1,ratio[i]);
        axrange axtest = {0,10,0,10,0,1,";;"};
        axrange axsigma = {0,11,0,1,0,1,";P_{given}[#times 10^{13}a.u];correct_ratio"};
        st.Graph(sigmagraph,axsigma);
        
        //st.GraphErrors(sigmagraph,axsigma);
        sigmagraph -> Draw("AP");*/
    }
    axrange axres = {0,11,0,100,0,1,";P_{given}[K*kHz];P_{fit}/#DeltaP_{fit}[K*kHz]"};
    st.Graph(resgraph,axres);
    resgraph -> Draw("AP");
    //TF1* f = new TF1("f","x",0,11);
    //f -> Draw("same");
}