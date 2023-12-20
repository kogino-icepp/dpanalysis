#include <iostream>
#include <queue>
#include "../headers/fitter.h"
#include "../headers/mask_map.h"
using namespace std;
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
const double scadnu = 0.033333;
const double dNu=0.0000885;//周波数分解能
const double kb=1.38*pow(10,-23);
const double df=88.5*pow(10,3);
const double Tc=76;
const double Th=297;
const double DINF=1e9;
const double DeltaP = 0;
const double Aeff = 0.385;
//特定のデータだけ引っ張ってくる用
const vector<int> lsbo = {1,3,5,13,15,17};
const vector<int> lsbi = {2,4,6,14,16,18};
const vector<int> usbi = {7,9,11,19,21,23};
const vector<int> usbo = {8,10,12,20,22,24};
string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
string savedirp = "/Users/oginokyousuke/data/peak_data/";
string saveexe = "/Users/oginokyousuke/data/search_exe/";
Fitter ft;
axrange axscale = {0,1,0,1,0,1,"test_data;xscale;yscale"};
vector<int> sbin = {0,512,-512,-1024};
Color_t vcolor[5] = {kBlue,kRed,kGreen,kCyan,kMagenta};//色変えて複数パターンお絵描きしたい時の配列
double psigma[4] = {0.18068,0.183154,0.194241,0.186464};
double psigma2[4] = {0.162231,0.161409,0.172565,0.166457};//ベースラインとの合わせ技でフィット　した場合の温度幅　　ここで作ったσでnormalizeすると原義1に本当になるのか
double wsigma[4] = {0.066,0.066,0.066,0.066};

//前もって隠しておく

//フィットに使う諸々の関数群
Double_t chiF_free(double x,double p0,double k,double p1){
    return p0*TMath::Gamma(k/2,x/(2*p1));
}
Double_t chiF_freefit(double x,double p0,double p1,double k,double bin){
    return (chiF_free((x+bin/2),p0,k,p1)-chiF_free((x-bin/2),p0,k,p1));
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
    else return 0;
}
//スケールするとそのままだとdnuの影響でピーク0になるので修正
/*
要件定義
1. 中身にあまりにも色々な定数が使われていて全部作り直しは不可能
2. 周波数の情報(両端の値など)を渡して中で適切に変換を行いフィットしてくれるような関数に作り替えたい
*/
double F_sigscale(double x,double P,double r,double fmin){
    double bin = x/0.0344827586;
    
    double f = fmin+bin*dnu;//ここをxとfminの関数に変える
    double f0 = fmin+10*dnu;//ここは自明ではあるがx0からf0へ
    return F_sig2(f,f0,P,r);
}
double CoupConst2(double p,double dp){
    if(p<0)return 2*1.3458*4.5*pow(10,-14)*sqrt(dp*1.96*pow(10,23))*sqrt(1/Aeff);
    else return 2*1.3458*4.5*pow(10,-14)*sqrt((dp*1.96+p)*pow(10,23))*sqrt(1/Aeff);
}
//ノーマライズされていることは前提としてそのカイ二乗値のp値をお手軽に返してくれる関数
double PValue(int k,double x){
    return 1-TMath::Gamma(k/2,(k/2)*x);
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
void ChiCheck2(int i,int j,int p,TH1D*hist){
    int xfft = XFFT(i);
    Mask ms;
    vector<vector<int>> mskmap = ms.maskmap;
    bool mskhantei[4][nbin];
    for(int j=0;j<4;j++){
        for(int bin=0;bin<nbin;bin++){mskhantei[j][bin] = false;}
    }
    for(int j=0;j<4;j++){//mskmapに記録されているものだけtrueに変更
        for(int bin:mskmap[j])mskhantei[j][bin] = true;
    }
    filesystem::path path=filesystem::current_path();
    filesystem::current_path(saveexe);
    string fname = "allbinbase"+to_string(i)+"_"+to_string(j)+".root";
    string tname = "test_tree"+to_string(p);
    const char* filename = fname.c_str();
    const char* treeName = tname.c_str();
    TFile*file = TFile::Open(filename);
    if (!file) {
        cerr << "Error opening file " << filename << std::endl;
        return;
    }
    // Get the TTree
    TTree*tree = dynamic_cast<TTree*>(file->Get(treeName));
    if (!tree) {
        cerr << "Error getting tree " << treeName << " from file " << filename << std::endl;
        file->Close();
        return;
    }
    Int_t numEntries = tree->GetEntries();
    vector<vector<double>> vparas(3,vector<double>(nbin));
    vector<double> vpfreq(nbin,DINF);
    vector<int> vbin(nbin);
    rep(i,numEntries){
        tree -> GetEntry(i);
        Double_t a, b, c, chi,freq;
        int bin;
        tree->SetBranchAddress("a", &a);
        tree->SetBranchAddress("b", &b);
        tree->SetBranchAddress("c", &c);
        tree->SetBranchAddress("chi", &chi);
        tree->SetBranchAddress("bin", &bin);
        tree->SetBranchAddress("freq",&freq);
        vparas[0][bin]=a;
        vparas[1][bin]=b;
        vparas[2][bin]=c;
        vpfreq[bin]=freq;
        if(chi==DINF)continue;
        hist -> Fill(chi);
        
    }
}

void GetBasicData(int i,int j,int p,TGraph*prec){
    filesystem::path path=filesystem::current_path();
    string cdir=dir+"band"+to_string(i);
    filesystem::current_path(cdir);
    Setting st;
    st.dot_size=0.8;
    st.markerstyle=20;
    st.color = kGreen;
    st.lcolor = kBlue;
    double Freq[nbin],cold[nbin],hot[nbin],mirror[nbin];
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
    rep(bin,nbin){
        long double gain=(hot[bin]-cold[bin])/(2*kb*(Th-Tc)*df);
        double psys=(cold[bin]/gain)-2*kb*Tc*df;
        
        prec -> SetPoint(bin-sb,Freq[bin],((mirror[bin]/gain)-psys)/(2*kb*df));
        //cout << Freq[bin] << " " <<((mirror[bin]/gain)-psys)/(2*kb*df) << endl;
        
    }
    return;
}
//グローバルにエラーをつけてフィッティング
void gloGetDPfit(int i,int j,int p,double (&dlist)[nbin],double (&deltaP)[nbin]){
    Mask ms;
    vector<vector<int>> mskmap = ms.maskmap;
    bool mskhantei[4][nbin];
    for(int j=0;j<4;j++){
        for(int bin=0;bin<nbin;bin++){mskhantei[j][bin] = false;}
    }
    for(int j=0;j<4;j++){//mskmapに記録されているものだけtrueに変更
        for(int bin:mskmap[j])mskhantei[j][bin] = true;
    }
    filesystem::path path=filesystem::current_path();
    filesystem::current_path(saveexe);
    Setting st;
    string fname = "allbinbase"+to_string(i)+"_"+to_string(j)+".root";
    string tname = "test_tree"+to_string(p);
    const char* filename = fname.c_str();
    const char* treeName = tname.c_str();
    TFile*file = TFile::Open(filename);
    if (!file) {
        cerr << "Error opening file " << filename << std::endl;
        return;
    }
    // Get the TTree
    TTree*tree = dynamic_cast<TTree*>(file->Get(treeName));
    if (!tree) {
        cerr << "Error getting tree " << treeName << " from file " << filename << std::endl;
        file->Close();
        return;
    }
    Int_t numEntries = tree->GetEntries();
    vector<vector<double>> vparas(3,vector<double>(nbin));
    vector<double> vpfreq(nbin,DINF);
    vector<double> vpchi(nbin,DINF);
    vector<int> vbin(nbin);
    rep(i,numEntries){
        tree -> GetEntry(i);
        Double_t a, b, c, chi,freq;
        int bin;
        tree->SetBranchAddress("a", &a);
        tree->SetBranchAddress("b", &b);
        tree->SetBranchAddress("c", &c);
        tree->SetBranchAddress("chi", &chi);
        tree->SetBranchAddress("bin", &bin);
        tree->SetBranchAddress("freq",&freq);
        vparas[0][bin]=a;
        vparas[1][bin]=b;
        vparas[2][bin]=c;
        vpfreq[bin]=freq;
        vpchi[bin] = chi;
    }
    file -> Close();
    TGraph* pgraph = new TGraph;
    GetBasicData(i,j,p,pgraph);
    int xfft = XFFT(i);
    double rtn = 0;
    double gsigma = 0;
    double samplenum = 0;
    prep(bin,sb,fb){
        if(vparas[0][bin]==DINF &&vparas[1][bin]==DINF  &&vparas[2][bin]==DINF)continue;
        bool togemask = false;
        rep(ad,29){
            if(mskhantei[xfft][bin+ad]){
                togemask = true;
                break;
            }
        }
        if(togemask)continue;
        double sfreq,ffreq;
        if(i%2==1){
            sfreq = pgraph -> GetPointX(bin-sb);
            ffreq = pgraph -> GetPointX(bin-sb+29);
        }
        else{
            ffreq = pgraph -> GetPointX(bin-sb);
            sfreq = pgraph -> GetPointX(bin-sb+29);
        }
        TF1* peakquad = new TF1("peakquad","[0]*(x-[1])*(x-[1])+[2]+F_sig2(x,[3],[4],0.5)",sfreq,ffreq);
        TGraphErrors* spgraph = new TGraphErrors;
        double yMin,yscale;
        double xmin;
        ft.make_scale2(pgraph,spgraph,bin-sb,yMin,yscale);
        ft.rescale_para(vparas[0][bin],vparas[1][bin],vparas[2][bin],sfreq,yMin,ffreq-sfreq,yscale,peakquad);
        double a = peakquad -> GetParameter(0);
        double b = peakquad -> GetParameter(1);
        double c = peakquad -> GetParameter(2);
        rep(k,dbin){
            double x = pgraph -> GetPointX(bin-sb+k);
            double y = pgraph -> GetPointY(bin-sb+k);
            if(k<10 || k>=20){
                gsigma += (y-(a*(x-b)*(x-b)+c))*(y-(a*(x-b)*(x-b)+c));
                samplenum++;
            }
        }
    }
    gsigma = sqrt(gsigma/(samplenum-1));
    for(int bin=sb;bin<fb;bin++){
        if(vparas[0][bin]==DINF &&vparas[1][bin]==DINF  &&vparas[2][bin]==DINF)continue;
        bool togemask = false;
        rep(ad,29){
            if(mskhantei[xfft][bin+ad]){
                togemask = true;
                break;
            }
        }
        if(togemask)continue;
        double sfreq,ffreq,mfreq;
        int bpos;
        if(i%2==1){
            sfreq = pgraph -> GetPointX(bin-sb);
            mfreq = pgraph -> GetPointX(bin-sb+10);
            ffreq = pgraph -> GetPointX(bin-sb+29);
            bpos = 10;
        }
        else{
            ffreq = pgraph -> GetPointX(bin-sb);
            mfreq = pgraph -> GetPointX(bin-sb+19);
            sfreq = pgraph -> GetPointX(bin-sb+29);
            bpos = 19;
        }
        TF1* peakquad = new TF1("peakquad","[0]*(x-[1])*(x-[1])+[2]+F_sig2(x,[3],[4],0.5)",sfreq,ffreq);
        //F_sigscale(double x,double P,double r,double fmin)
        TF1* scalepeak = new TF1("scalepeak","[0]*(x-[1])*(x-[1])+[2]+F_sigscale(x,[3],0.5,[4])",0,1);
        //TF1* cubfunc = new TF1("cubfunc","[0]*x*x*x+[1]*x*x+[2]*x+[3]",sfreq,ffreq);
        TGraphErrors* spgraph = new TGraphErrors;
        double yMin,yscale;
        double xmin;
        axrange axfit = {sfreq,ffreq,40,80,0,1,"fitgraph;Freq[GHz];Prec[K]"};
        ft.make_scale2(pgraph,spgraph,bin-sb,yMin,yscale);
        ft.rescale_para(vparas[0][bin],vparas[1][bin],vparas[2][bin],sfreq,yMin,ffreq-sfreq,yscale,peakquad);//ここで
        axrange axtest = {sfreq,ffreq,0,100,0,1,"SignalFit;Freq[GHz];Prec[K]"};
        //一からプロットし直してエラーをグローバルにつけ
        //パラメータ自体はscaleした後のデータを格納しているのでそのまま渡せばいいのでは？？
        scalepeak -> SetParameter(0,vparas[0][bin]);
        scalepeak -> SetParameter(1,vparas[1][bin]);
        scalepeak -> SetParameter(2,vparas[2][bin]);
        rep(k,dbin)spgraph -> SetPointError(k,0,gsigma/yscale);

        st.GraphErrors(spgraph,axscale);
        double smfreq = spgraph -> GetPointX(bpos);
        peakquad -> FixParameter(3,mfreq);
        scalepeak -> FixParameter(4,sfreq);
        rep(ite,10)spgraph -> Fit(scalepeak,"Q0","",0,1);
        rep(ite,10)spgraph -> Fit(scalepeak,"MQ0","",0,1);
        
        //フィットがある程度収束するまでこれ続ける
        //fitgraph -> Fit(peakquad,"EQ","",sfreq,ffreq);
        spgraph -> Fit(scalepeak,"EQ","",0,1);
        st.GraphErrors(spgraph,axscale);
        spgraph -> Draw("AP");

        double spout = scalepeak -> GetParameter(3);
        double sdpout = scalepeak -> GetParError(3);
        spout *= yscale;
        sdpout *= yscale;
        dlist[bin] = spout;
        deltaP[bin] = sdpout;
        //cout << "P/dP: " << spout/sdpout << endl;
        //hist -> Fill(spout);*/
    }
    cout << "gsigma = " << gsigma << endl;
    return;
}
//scaleした後のフィット結果を元の値に復元する、できるかな？
void GetDPfit(int i,int j,int p,double (&dlist)[nbin],double (&deltaP)[nbin],TH1D* hist,double &sighosei){
    Mask ms;
    vector<vector<int>> mskmap = ms.maskmap;
    bool mskhantei[4][nbin];
    for(int j=0;j<4;j++){
        for(int bin=0;bin<nbin;bin++){mskhantei[j][bin] = false;}
    }
    for(int j=0;j<4;j++){//mskmapに記録されているものだけtrueに変更
        for(int bin:mskmap[j])mskhantei[j][bin] = true;
    }
    
    filesystem::path path=filesystem::current_path();
    filesystem::current_path(saveexe);
    Setting st;
    st.dot_size=0.8;
    st.markerstyle=20;
    st.color = kGreen;
    st.lcolor = kGreen;
    
    string fname = "allbinbase"+to_string(i)+"_"+to_string(j)+".root";
    string tname = "test_tree"+to_string(p);
    const char* filename = fname.c_str();
    const char* treeName = tname.c_str();
    TFile*file = TFile::Open(filename);
    if (!file) {
        cerr << "Error opening file " << filename << std::endl;
        return;
    }
    // Get the TTree
    TTree*tree = dynamic_cast<TTree*>(file->Get(treeName));
    if (!tree) {
        cerr << "Error getting tree " << treeName << " from file " << filename << std::endl;
        file->Close();
        return;
    }
    Int_t numEntries = tree->GetEntries();
    vector<vector<double>> vparas(3,vector<double>(nbin));
    vector<double> vpfreq(nbin,DINF);
    vector<double> vpchi(nbin,DINF);
    vector<int> vbin(nbin);
    rep(i,numEntries){
        tree -> GetEntry(i);
        Double_t a, b, c, chi,freq;
        int bin;
        tree->SetBranchAddress("a", &a);
        tree->SetBranchAddress("b", &b);
        tree->SetBranchAddress("c", &c);
        tree->SetBranchAddress("chi", &chi);
        tree->SetBranchAddress("bin", &bin);
        tree->SetBranchAddress("freq",&freq);
        vparas[0][bin]=a;
        vparas[1][bin]=b;
        vparas[2][bin]=c;
        vpfreq[bin]=freq;
        vpchi[bin] = chi;
    }
    file -> Close();
    TGraph* pgraph = new TGraph;
    GetBasicData(i,j,p,pgraph);

    //vectorにとりあえずのフィット結果を詰める
    //どこかにグローバルなエラーと局所的なエラーを両方出すプログラムが欲しい
    int xfft = XFFT(i);
    int outbin[10];
    int outindex = 0;
    
    int bin95 = 0;
    int index = 0;
    TH1D* firsthist = new TH1D("firsthist",";P_{fit}/#Delta P_{fit};Count",100,-10,10);
    
    /*
    要件定義
    最終結果ならびに結果の精査に必要な情報を保存するためにrootファイルを作成して保存する機能を追加する
    ・保存したいパラメータ：bin,周波数,Pfit,dPfit,chi/ndf
    ・注意事項：rootファイルを開いた状態で別のアレを開くとアレするらしいので競合しないように適宜closeしたり順番を工夫する
    ・
    */
    for(int bin=sb;bin<fb;bin++){
        if(vparas[0][bin]==DINF &&vparas[1][bin]==DINF  &&vparas[2][bin]==DINF){
            continue;
        }
        bool togemask = false;
        rep(ad,29){
            if(mskhantei[xfft][bin+ad]){
                togemask = true;
                break;
            }
        }
        if(togemask)continue;
        double sfreq,ffreq,mfreq;
        int bpos;
        if(i%2==1){
            sfreq = pgraph -> GetPointX(bin-sb);
            mfreq = pgraph -> GetPointX(bin-sb+10);
            ffreq = pgraph -> GetPointX(bin-sb+29);
            bpos = 10;
        }
        else{
            ffreq = pgraph -> GetPointX(bin-sb);
            mfreq = pgraph -> GetPointX(bin-sb+19);
            sfreq = pgraph -> GetPointX(bin-sb+29);
            bpos = 19;
        }
        TF1* peakquad = new TF1("peakquad","[0]*(x-[1])*(x-[1])+[2]+F_sig2(x,[3],[4],0.5)",sfreq,ffreq);
        //F_sigscale(double x,double P,double r,double fmin)
        TF1* scalepeak = new TF1("scalepeak","[0]*(x-[1])*(x-[1])+[2]+F_sigscale(x,[3],0.5,[4])",0,1);
        //TF1* cubfunc = new TF1("cubfunc","[0]*x*x*x+[1]*x*x+[2]*x+[3]",sfreq,ffreq);
        TGraphErrors* spgraph = new TGraphErrors;
        double yMin,yscale;
        double xmin;
        axrange axfit = {sfreq,ffreq,40,80,0,1,"fitgraph;Freq[GHz];Prec[K]"};
        ft.make_scale2(pgraph,spgraph,bin-sb,yMin,yscale);
        ft.rescale_para(vparas[0][bin],vparas[1][bin],vparas[2][bin],sfreq,yMin,ffreq-sfreq,yscale,peakquad);//ここで
        axrange axtest = {sfreq,ffreq,0,100,0,1,"SignalFit;Freq[GHz];Prec[K]"};
        //一からプロットし直してエラーをグローバルにつける
        TGraphErrors* fitgraph = new TGraphErrors;
        TGraphErrors* fitgraph2 = new TGraphErrors;
        //パラメータ自体はscaleした後のデータを格納しているのでそのまま渡せばいいのでは？？
        scalepeak -> SetParameter(0,vparas[0][bin]);
        scalepeak -> SetParameter(1,vparas[1][bin]);
        scalepeak -> SetParameter(2,vparas[2][bin]);
        
        double a = peakquad -> GetParameter(0);
        double b = peakquad -> GetParameter(1);
        double c = peakquad -> GetParameter(2);
        double sigma = 0;
        rep(k,dbin){
            double x = pgraph -> GetPointX(bin-sb+k);
            double y = pgraph -> GetPointY(bin-sb+k);
            fitgraph -> SetPoint(k,x,y);
            if(k<10 || k>=20){
                sigma += (y-(a*(x-b)*(x-b)+c))*(y-(a*(x-b)*(x-b)+c));
            }
        }
        //scalepeak -> SetParameter(4,1);
        sigma /= 17;
        sigma = sqrt(sigma);
        rep(k,dbin){
            fitgraph -> SetPointError(k,0,sigma);
            spgraph -> SetPointError(k,0,sigma/yscale);
        }
        st.GraphErrors(spgraph,axscale);
        //spgraph -> Draw("AP");
        //scalepeak -> Draw("same");
        st.GraphErrors(fitgraph,axtest);
        st.GraphErrors(fitgraph2,axtest);
        double smfreq = spgraph -> GetPointX(bpos);
        //cout << "smfreq: " << smfreq << endl;
        //if(sigma>0.2)cout << bin << " " << sigma << endl;
        //fitgraph -> Draw("AP");
        peakquad -> FixParameter(3,mfreq);
        scalepeak -> FixParameter(4,sfreq);
        rep(ite,10){
            //fitgraph -> Fit(peakquad,"Q0","",sfreq,ffreq);
            spgraph -> Fit(scalepeak,"Q0","",0,1);
        }
        rep(ite,10){
            //fitgraph -> Fit(peakquad,"MQ0","",sfreq,ffreq);
            spgraph -> Fit(scalepeak,"MQ0","",0,1);
        }
        //フィットがある程度収束するまでこれ続ける
        //fitgraph -> Fit(peakquad,"EQ","",sfreq,ffreq);
        spgraph -> Fit(scalepeak,"EQ","",0,1);
        //spgraph -> Draw("AP");
        //scalepeak -> Draw("same");
        //fitgraph -> Draw("AP");
        //peakquad -> Draw("same");
        
       // double pout = peakquad -> GetParameter(4);
        //double dpout = peakquad -> GetParError(4);
        double spout = scalepeak -> GetParameter(3);
        double sdpout = scalepeak -> GetParError(3);
        spout *= yscale;
        sdpout *= yscale;
        dlist[bin] = spout;
        deltaP[bin] = sdpout;
        index++;
        firsthist -> Fill(spout/(sdpout));
    }
    double fmin = 213.5+i*2;
    double fmax = 216.5+i*2;
    TF1* fgaus = new TF1("fgaus","gaus");
    firsthist -> Fit(fgaus);
    st.Hist(firsthist);
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    c1 -> SetLogy();
    firsthist -> Draw();
    double sigma = fgaus -> GetParameter("Sigma");
    sighosei = sigma;
    rep(k,nbin){
        if(dlist[k]!=DINF)hist -> Fill(dlist[k]/(deltaP[k]*sigma));
    }
}


void MakeLimit(double (&dlist)[8][nbin],double (&deltaP)[8][nbin],int i){
    //獲得したpoutのデータ一覧から平均した値を出力+配列の確保を行う
    //ここで統計的に薄めて過度なexcessがないかどうかを確認する
    st.lcolor = kBlue;
    st.dot_size = 0.8;
    int xfft = XFFT(i);
    int shift[4];
    if(xfft%2==1)shift[0] = 0,shift[1]=512,shift[2]=-512,shift[3]=-1024;
    else shift[0] = 0,shift[1]=-512,shift[2]=512,shift[3]=1024;
    TGraphErrors* sigraph = new TGraphErrors;
    TGraph* glimit = new TGraph;
    int gbin = 0;
    //deltaPの平均化
    //TH1* delhist = nullptr;
    prep(bin,sb,fb){
        int num = 0;
        double dlim = 0;
        double vardeltaP = 0;
        double sigsum = 0;
        double stdsig = 0;
        rep(ite,8){
            int j = ite/2;
            if(dlist[ite][bin+shift[j]] != DINF){
                num++;
                dlim += dlist[ite][bin+shift[j]];
                vardeltaP += deltaP[ite][bin+shift[j]];
                //σの方もこちらで計算しておく
                sigsum += dlist[ite][bin+shift[j]]/deltaP[ite][bin+shift[j]];
            }
        }
        sigsum /= num;
        dlim /= num;
        vardeltaP /= num;
        rep(ite,8){
            int j = ite/2;
            if(dlist[ite][bin+shift[j]] != DINF)stdsig += pow(sigsum-dlist[ite][bin+shift[j]]/deltaP[ite][bin+shift[j]],2);
        }
        stdsig = sqrt(stdsig/num);
        
        
        if(num>0){
            //横軸周波数◯ →　縦軸を[K]からχに変換したい
            double freq;
            //周波数の導出
            if(i%2==1)freq = (213.8+i*2)+0.0000762939*bin;
            else freq = (216.2+i*2)-0.0000762939*bin;
            //limの計算、正負で対応が異なる
            double chilim;
            if(dlim>0){
                chilim = PtoChi((dlim+1.96*vardeltaP)*2*kb*df);
            }
            else{
                chilim = PtoChi(1.96*vardeltaP*2*kb*df);
            }
            dlist[0][bin] = chilim;
            
            //glimit -> SetPoint(gbin,freq,chilim);
            //sigraph -> SetPoint(gbin,freq,sigsum);
            //sigraph -> SetPointError(gbin,0,stdsig);
            //cout << chilim << " " << freq << endl;
            //cout << gbin << " : " << freq << " : " << sigsum << endl;
            gbin++;
        }
    }
    //cout << "fit num: " << 1750 << " ~ "  << 2000<< endl;
    /*TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    double fmin = 213.5+i*2;
    double fmax = 216.5+i*2;
    
    axrange axsig = {fmin,fmax,-10,10,0,1,";Freq;P_{fit}/#Delta P_{fit}"};
    axrange axlim = {fmin,fmax,pow(10,-11),pow(10,-7),0,1,"#chi"+to_string(i)+";Freq[GHz];#chi"};
    st.Graph(glimit,axlim);
    st.GraphErrors(sigraph,axsig);
    glimit -> SetLineColor(kBlue);
    c1 -> SetLogy();
    glimit -> Draw("AL");
    cout << "sucess " << endl;
    //sigraph -> Draw("AL");*/
    //横軸の周波数だけなんとか引っ張れるかどうか
}
//一旦現実逃避で、全データに対してプロットを作るプログラムを作成する。
void alllimit(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    c1 -> SetLogy();
    TGraph* limgraph = new TGraph;//limitのお絵描きをするフィールド
    axrange axall = {215,265,pow(10,-11),pow(10,-8),0,1,";Freq[GHz];"};
    int allbin = 0;
    st.lcolor = kBlue;
    st.color = kBlue;
    int lbin = 0;
    prep(i,1,25){
        cout << i << endl;
        double deltaP[8][nbin];//フィットで得られたエラーを格納するための二次配列
        double pfitlist[8][nbin];
        rep(k,8)rep(bin,nbin)pfitlist[k][bin] = DINF;
            
        rep(j,4){
            prep(p,1,3){
                double testdouble;
                TH1D* testhist = new TH1D("testhist",";P/#Delta P;Count",100,-1,1);
                //ここではlistに値さえ格納できれば良い
                GetDPfit(i,j,p,pfitlist[2*j+(p-1)],deltaP[2*j+(p-1)],testhist,testdouble);
                //testhist -> Draw();
                //testhist -> Fit("gaus");
            }
        }
        MakeLimit(pfitlist,deltaP,i);
        if(i%2==1){
            prep(bin,sb,fb){
                double freq;
                //周波数の導出
                if(i%2==1)freq = (213.8+i*2)+0.0000762939*bin;
                else freq = (216.2+i*2)-0.0000762939*bin;
                if(pfitlist[0][bin]!=DINF){
                    limgraph -> SetPoint(lbin,freq,pfitlist[0][bin]);
                    lbin++;
                }
            }
        }
        else{
            rprep(bin,fb,sb){
                double freq;
                //周波数の導出
                if(i%2==1)freq = (213.8+i*2)+0.0000762939*bin;
                else freq = (216.2+i*2)-0.0000762939*bin;
                if(pfitlist[0][bin]!=DINF){
                    limgraph -> SetPoint(lbin,freq,pfitlist[0][bin]);
                    //cout << "bin: " << lbin << " | freq: " << freq << " <=> " << pfitlist[0][bin] << endl;
                    lbin++;
                }
            }
        }
        
    }

    st.Graph(limgraph,axall);
    limgraph -> Draw("AL");
}
//これがメイン関数
void scale_peakfit(){
    //まずはそれぞれ別々の情報が詰められているか確認する
    double maxbin = 0.5;
    string whitedir = "/Users/oginokyousuke/data/white_noise";
    vector<int> excess;
    //alllimit();
    for(int fn=15;fn<16;fn++){
        TH1D* chihist = new TH1D("chihist","chihist;#chi^{2}/NDF;Count",100,0,5);
        double deltaP[8][nbin];
        double pfitlist[8][nbin];
        double gdeltaP[8][nbin];
        double gpfitlist[8][nbin];
        rep(i,8)rep(j,nbin){
            pfitlist[i][j] = DINF;
            gpfitlist[i][j] = DINF;
        }
        string roofilename = "peakfitdata"+to_string(fn)+".root";
        //TFile * savefile = new TFile(roofilename.c_str(),"recreate");
        double sighosei[8];
        for(int j=0;j<4;j++){
            prep(p,1,3){
                double testlist[nbin];
                double testdeltaP[nbin];
                double gtestlist[nbin];
                double gtestdeltaP[nbin];
                rep(bin,nbin)testlist[bin] = DINF;
                TH1D* whitehist1 = new TH1D("whitehist1",";P_{fit}[kHz*W];Count",100,-1,1);
                TH1D* whitehist2 = new TH1D("whitehist2",";P_{fit}[kHz*W];Count",100,-1,1);
                TH1D* scalehist = new TH1D("scalehist",";P_{fit}/#Delta P_{fit};Count",100,-10,10);
                TF1* fgaus = new TF1("fgaus","gaus",-10,10);
                TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
                c1 -> SetMargin(0.14,0.11,0.2,0.1);
                //fitterを毎回回さなくてもいいように確定版のデータでなくてもいいのでrootファイルを作成して保存しておきたい
                GetDPfit(fn,j,p,testlist,testdeltaP,scalehist,sighosei[j*2+(p-1)]);
                //c1 -> SetLogy();
                //st.Hist(scalehist);
                //scalehist -> Draw();
                /*gloGetDPfit(fn,j,p,gtestlist,gtestdeltaP);
                TGraph* gfitratio = new TGraph;
                TGraph* gefitratio = new TGraph;
                int rbin = 0;
                prep(bin,sb,fb){
                    if(testlist[bin]==DINF)continue;
                    gfitratio -> SetPoint(rbin,bin,testlist[bin]/gtestlist[bin]);
                    cout << bin << " : " << testlist[bin]/gtestlist[bin] << endl;
                    rbin++;
                }
                axrange axratio = {0,nbin,0,2,0,1,";Bin;Ratio"};
                st.Graph(gfitratio,axratio);
                gfitratio -> Draw("AP");
                //このfor文内でコメントアウトするときにはここを使おう！*/
            }
        }
        cout << fn << endl;
        rep(bin,8)cout << "hosei: " << sighosei[bin] << endl;
        
    }
    return;
}