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
const double dfnarrow = 17.7*pow(10,3);
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
double hosei[24][8] = {
    {1.13245,1.142,1.13562,1.1432,1.14706,1.11157,1.1375,1.15458},//1
    {1.13183,1.14521,1.13852,1.14507,1.12591,1.13597,1.12438,1.14009},
    {1.14608,1.12456,1.12808,1.14717,1.13378,1.14064,1.11851,1.12468},
    {1.12888,1.13746,1.14985,1.13999,1.13494,1.12805,1.13739,1.13202},
    {1.11944,1.14097,1.14338,1.14011,1.14622,1.14715,1.12934,1.13564},//5
    {1.14708,1.15226,1.14045,1.13297,1.1338,1.14008,1.12353,1.13136},
    {1.13668,1.12672,1.14737,1.14072,1.12925,1.12689,1.13128,1.12511},
    {1.12074,1.11436,1.14591,1.11498,1.11501,1.14038,1.12464,1.13384},
    {1.1509,1.13773,1.13295,1.14587,1.15243,1.1385,1.14089,1.13442},
    {1.14157,1.13675,1.15343,1.11896,1.12901,1.13417,1.14798,1.13203},//10
    {1.14526,1.13279,1.14707,1.14046,1.14375,1.14713,1.13475,1.1641},
    {1.1603,1.14312,1.1541,1.13668,1.15352,1.13539,1.14394,1.136},
    {1.15067,1.13055,1.14324,1.13867,1.1304,1.14212,1.12163,1.13575},//13(ほとんど空)
    {1.15356,1.13391,1.13923,1.13294,1.14626,1.14372,1.14041,1.14722},
    {1.12344,1.13704,1.15244,1.15078,1.15253,1.15252,1.16004,1.16154},//15(data0だけ破損)
    {1.13971,1.14376,1.15491,1.1535,1.13437,1.17363,1.16383,1.12419},
    {1.15111,1.14801,1.1397,1.1435,1.14912,1.15738,1.13241,1.13482},
    {1.15442,1.15381,1.1429,1.17753,1.15116,1.15322,1.15399,1.15426},
    {1.14752,1.1495,1.1605,1.1463,1.16988,1.15173,1.15258,1.14483},
    {1.13979,1.17492,1.13996,1.14347,1.15231,1.14718,1.1463,1.14697},//20
    {1.14357,1.15448,1.15859,1.16217,1.16527,1.14983,1.14865,1.15189},
    {1.14171,1.14207,1.14164,1.15007,1.14245,1.13256,1.14801,1.14365},
    {1.14133,1.16191,1.14414,1.14913,1.14094,1.16344,1.14205,1.15776},
    {1.13344,1.14777,1.16311,1.14814,1.14276,1.15348,1.16192,1.16808}//24
};
//前もって隠しておく
pair<double,double> MeanError(vector<double>data){
    double mean = 0;
    double num = 0;
    for(auto v:data){
        mean += v;
        num++;
    }
    if(num==0)return {mean,0};
    mean /= num;
    double rtn = 0;
    for(auto v:data)rtn += (mean-v)*(mean-v);
    rtn = sqrt(rtn/(num));
    return {mean,rtn};
}
double FtoMass(double freq){
    return freq*pow(10,9)*4.1357*pow(10,-15)*pow(10,6);
}
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
    return 4.5*pow(10,-14)*sqrt(P*pow(10,23)*2*1.208)*sqrt(1/0.279);
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
//[W]の単位でP_95%を算出する関数,limitではない
double P95(double pfit,double deltap){
    if(pfit>0)return (pfit+1.96*deltap)*2*kb*df;
    else return 1.96*deltap*2*kb*df;
}
double P95narrow(double pfit,double deltap){
    if(pfit>0)return (pfit+1.96*deltap)*2*kb*dfnarrow;
    else return 1.96*deltap*2*kb*dfnarrow;
}
double Freq(int i,int bin){
    if(i%2==1)return (213.8+i*2)+0.0000762939*bin;
    else return (216.2+i*2)-0.0000762939*bin;
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
        //peakquad -> FixParameter(3,mfreq);
        scalepeak -> FixParameter(4,sfreq);
        scalepeak -> SetParameter(3,0.1);
        
        rep(ite,30)spgraph -> Fit(scalepeak,"Q0","",0,1);
        rep(ite,30)spgraph -> Fit(scalepeak,"MQ0","",0,1);
        
        //フィットがある程度収束するまでこれ続ける
        //fitgraph -> Fit(peakquad,"EQ","",sfreq,ffreq);
        rep(ite,100){
            spgraph -> Fit(scalepeak,"EQ","",0,1);
            double sdpout = scalepeak -> GetParError(3);
            if(sdpout*yscale<1)break;
        }
        double pfitmae = scalepeak -> GetParameter(3);
        cout << pfitmae << endl;
        st.GraphErrors(spgraph,axscale);
        spgraph -> Draw("AP");

        double spout = scalepeak -> GetParameter(3);
        double sdpout = scalepeak -> GetParError(3);
        double chi = scalepeak -> GetChisquare();
        double ndf = scalepeak -> GetNDF();
        spout *= yscale;
        sdpout *= yscale;
        dlist[bin] = spout;
        deltaP[bin] = sdpout;
        //cout << "P/dP: " << spout/sdpout << endl;
        //hist -> Fill(spout);*/
    }
    
    return;
}
//scaleした後のフィット結果を元の値に復元する、できるかな？
void GetDPfit(int i,int j,int p,double (&dlist)[nbin],double (&deltaP)[nbin]){
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

    //どこかにグローバルなエラーと局所的なエラーを両方出すプログラムが欲しい
    int xfft = XFFT(i);
    TH1D* firsthist = new TH1D("firsthist",";P_{fit}/#Delta P_{fit};Count",100,-10,10);
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
        //peakquad -> FixParameter(3,mfreq);
        scalepeak -> FixParameter(4,sfreq);
        scalepeak -> SetParameter(3,0.1);
        
        rep(ite,30){
            //fitgraph -> Fit(peakquad,"Q0","",sfreq,ffreq);
            spgraph -> Fit(scalepeak,"Q0","",0,1);
        }
        rep(ite,30){
            //fitgraph -> Fit(peakquad,"MQ0","",sfreq,ffreq);
            spgraph -> Fit(scalepeak,"MQ0","",0,1);
        }
        
        //フィットがある程度収束するまでこれ続ける
        //fitgraph -> Fit(peakquad,"EQ","",sfreq,ffreq);
        rep(ite,100){
            spgraph -> Fit(scalepeak,"EQ","",0,1);
            double sdpout = scalepeak -> GetParError(3);
            if(sdpout*yscale<1)break;
        }
        double pfitmae = scalepeak -> GetParameter(3);
        //cout << pfitmae << endl;
        double spout = scalepeak -> GetParameter(3);
        double sdpout = scalepeak -> GetParError(3);
        spout *= yscale;
        sdpout *= yscale;
        dlist[bin] = spout;
        deltaP[bin] = sdpout;
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
            
            gbin++;
        }
    }
}
//一旦saveが成功した扱いでファイルを読み出すコードを作っておく
void ReadFile(int i,double (&plist)[8][nbin],double (&deltaP)[8][nbin]){
    filesystem::current_path("/Users/oginokyousuke/data/test/");
    string fname = "pfitres"+to_string(i)+".root";
    TFile* file = new TFile(fname.c_str());
    
    rep(j,8){
        string tname = "rtree"+to_string(j);
        TTree* tree = (TTree*)file->Get(tname.c_str());
        int binF;
        double pfitF,delpfitF;
        tree -> SetBranchAddress("bin",&binF);
        tree -> SetBranchAddress("pfit",&pfitF);
        tree -> SetBranchAddress("delpfit",&delpfitF);
        int entnum = tree -> GetEntries();
        rep(ite,entnum){
            tree -> GetEntry(ite);
            plist[j][binF] = pfitF;
            deltaP[j][binF] = delpfitF;
            //cout << binF << " : " << pfitF << " <=> " << delpfitF << endl;
        }
    }
    file -> Close();
    //本当は上で二次元配列に詰めて描画までやっちゃいたい(一旦別で関数を作ってみる)
}
void ReadGfile(int i,double (&pfitlist)[8][nbin],double (&deltaP)[8][nbin]){
    filesystem::current_path("/Users/oginokyousuke/data/test/");
    string fname = "pfitres"+to_string(i)+"global.root";
    TFile* file = new TFile(fname.c_str());
    
    rep(j,8){
        string tname = "rtree"+to_string(j);
        TTree* tree = (TTree*)file->Get(tname.c_str());
        int binF;
        double pfitF,delpfitF;
        tree -> SetBranchAddress("bin",&binF);
        tree -> SetBranchAddress("pfit",&pfitF);
        tree -> SetBranchAddress("delpfit",&delpfitF);
        int entnum = tree -> GetEntries();
        rep(ite,entnum){
            tree -> GetEntry(ite);
            pfitlist[j][binF] = pfitF;
            deltaP[j][binF] = delpfitF;
            //cout << binF << " : " << pfitF << " <=> " << delpfitF << endl;
        }
    }
    file -> Close();
}
//平均化したP/dPのグラフ作成の図 なんかエラーが大きい気がするのは気のせい？
void P_dP(){
    //ここを任意の区間まで拡張して書きたい(多分大丈夫、知らんけど何処かで誤魔化す)
    TGraphErrors* graph = new TGraphErrors;//描画グラフ
    TGraphErrors* sgraph = new TGraphErrors;//描画グラフ
    int pbin = 0;
    axrange axall = {215.8,264.2,-10,10,0,1,";Freq[GHz];P_{fit}/#DeltaP_{fit}"};
    axrange axsig = {215.8,264.2,0,2,0,1,";Freq[GHz];Sigma"};
    int si = 1;
    int fi = 25;
    prep(i,si,fi){
        TH1D* nhist = new TH1D("nhist",";P_{fit}/#DeltaP_{fit};Count",100,-10,10);
        cout << "band" << i << endl;
        int xfft = XFFT(i);
        double pfitlist[8][nbin],deltaP[8][nbin];
        rep(ite,8)rep(bin,nbin)pfitlist[ite][bin] = DINF;
        ReadFile(i,pfitlist,deltaP);
        int shift[4];
        if(xfft%2==1)shift[0] = 0,shift[1]=512,shift[2]=-512,shift[3]=-1024;
        else shift[0] = 0,shift[1]=-512,shift[2]=512,shift[3]=1024;
        prep(bin,sb,fb){
            int num = 0;
            vector<double> list;
            rep(ite,8){
                int j = ite/2;
                if(pfitlist[ite][bin+shift[j]]==DINF)continue;
                num++;
                list.push_back(pfitlist[ite][bin+shift[j]]/(deltaP[ite][bin+shift[j]]*hosei[i-1][ite]));
                //cout << pfitlist[ite][bin+shift[j]] << endl;
            }
            //if(num==0)continue;
            double prave,prerr;
            prave = MeanError(list).first;
            prerr = MeanError(list).second;
            double freq;
            //周波数の導出
            if(i%2==1)freq = (213.8+i*2)+0.0000762939*bin;
            else freq = (216.2+i*2)-0.0000762939*bin;
            graph -> SetPoint(pbin,freq,prave);
            nhist -> Fill(prave*sqrt(list.size()));
            graph -> SetPointError(pbin,0,prerr);
            /*if(prave+prerr>5 || prave-prerr<-5){
                cout << "num ->" << num << endl;
                cout << "out ->" << prave-prerr <<endl;
            }*/
            pbin++;
        }
        st.dot_size = 0.6;
        st.lcolor = kBlue;
        st.color = kBlue;
        st.markerstyle = 20;
        st.GraphErrors(graph,axall);
        graph -> Draw("AP");
        st.Hist(nhist);
        nhist -> Draw();
        TF1* fgaus = new TF1("fgaus","[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",-10,10);
        fgaus -> SetParameter(0,20000);
        fgaus -> SetParameter(1,0);
        fgaus -> SetParameter(2,1);

        nhist -> Fit(fgaus,"Q","",-10,10);
        double sigma = fgaus -> GetParameter(2);
        double dsigma = fgaus -> GetParError(2);
        sgraph -> SetPoint(i,215+2*i,sigma);
        sgraph -> SetPointError(i,0,dsigma);
        if(i==24)sgraph -> SetPointError(i,0,0.005);
    }
    st.dot_size = 0.7;
    st.color = kBlue;
    st.markerstyle = 20;
    st.GraphErrors(sgraph,axsig);
    sgraph -> Draw("AP");

}
void DrawPower(TH1D* hist){
    //ヒストグラムに詰めてどのくらいのパワーレベルになったか検証する
    double pfitlist[8][nbin],deltaP[8][nbin];
    rep(bin,nbin)pfitlist[0][bin] = DINF;
    ReadFile(1,pfitlist,deltaP);
    prep(bin,sb,fb){
        if(pfitlist[0][bin]==DINF)continue;
        double pw = pfitlist[0][bin]*2*kb*df;
        hist -> Fill(pw);
    }
}
//一旦bandごとにdrawできるプログラムを書く
void DrawLimit(){
    TGraph* chigraph = new TGraph;//limitのグラフ
    TGraph* chinarrow = new TGraph;
    TH1D* whitehist = new TH1D("whitehist",";Freq[GHz];#chi",627378,216,264);
    //axrange axall = {216,264,0,1,0,1,";Freq[GHz];coupling constant #chi"};
    axrange axall = {0,pow(10,5),0,pow(10,-8),0,1,";Freq[GHz];#DeltaP_{fit}[K]"};
    int chibin = 0;
    int chinbin = 0;
    
    TH1D* ehist = new TH1D("ehist",";#DeltaP_{fit}[K];Count",100,0.05,0.2);
    double xmin = 1000;
    double ymin = 1000;
    double xmax = -1;
    double ymax = 0;
    prep(i,1,25){
        double pfitlist[8][nbin],deltaP[8][nbin];
        double gpfitlist[8][nbin],gdeltaP[8][nbin];
        rep(ite,8)rep(bin,nbin){
            pfitlist[ite][bin] = DINF;
            gpfitlist[ite][bin] = DINF;
        }
        //ReadFile(i,pfitlist,deltaP);
        ReadFile(i,pfitlist,deltaP);
        int xfft = XFFT(i);
        int shift[4];
        if(xfft%2==1)shift[0] = 0,shift[1]=512,shift[2]=-512,shift[3]=-1024;
        else shift[0] = 0,shift[1]=-512,shift[2]=512,shift[3]=1024;
        vector<pair<double,double>> res;
        vector<pair<double,double>> resn;
        vector<pair<double,double>> resg;
        prep(bin,sb+40,fb-40){
            int num = 0;
            vector<double> pvec;
            vector<double> evec;
            vector<double> gpvec;
            vector<double> gevec;
            rep(ite,8){
                int j = ite/2;
                if(pfitlist[ite][bin+shift[j]]==DINF)continue;
                if(deltaP[ite][bin+shift[j]]>1)continue;
                num++;
                
                pvec.push_back(pfitlist[ite][bin+shift[j]]);
                evec.push_back(deltaP[ite][bin+shift[j]]);
                //gpvec.push_back(gpfitlist[ite][bin+shift[j]]);
                //gevec.push_back(gdeltaP[ite][bin+shift[j]]);
                //cout << pfitlist[ite][bin+shift[j]] << endl;
            }
            if(num==0)continue;
            //今回はpfitの平均値とdeltaPfitの平均値をそれぞれ値として使っているが本当にそれでいいのか？
            double pave,errave;
            double gpave,gerrave;
            pave = MeanError(pvec).first;
            errave = MeanError(evec).first;
            //gpave = MeanError(gpvec).first;
            //gerrave = MeanError(gevec).first;
            double freq;
            //周波数の導出
            if(i%2==1)freq = (213.8+i*2)+0.0000762939*bin;
            else freq = (216.2+i*2)-0.0000762939*bin;
            /*
            //もし横軸を周波数ではなくてmassにしたい場合はこのコメントアウトを外す
            freq = FtoMass(freq);*/
            res.push_back({(freq),PtoChi(P95(pave,errave))});
            //if(PtoChi(P95(pave,errave))>7*pow(10,-11))cout << pave << " " << errave << " " << freq << " <=> " << PtoChi(P95(pave,errave)) << " "  << bin << endl;
            //res.push_back({(freq),errave});
        }
        sort(res.begin(),res.end());
        //sort(resg.begin(),resg.end());
        for(auto v:res){
            //if(v.second>pow(10,-10))continue;
            if(v.second>0.16)continue;
            chigraph -> SetPoint(chibin,v.first,v.second);
            if(xmin>v.first)xmin = v.first;
            if(ymin>v.second)ymin = v.second;
            if(xmax<v.first)xmax = v.first;
            if(ymax<v.second)ymax = v.second;
            chibin++;
        }
        //cout << chibin << endl;
        /*for(auto v:resg){
            if(v.second>pow(10,-10))continue;
            grochigraph -> SetPoint(chinbin,v.first,v.second);
            chinbin++;
        }*/
    }
    
    st.Graph(chigraph,axall);
    chigraph -> SetLineColor(kBlack);
    chigraph -> Draw("AL");
    cout << "xmin: " << xmin << " , xmax: " << xmax << endl;
    cout << "ymin: " << ymin << " , ymax: " << ymax << endl;
    
}
void DrawRatio(int i,int j){
    double pfitlist[8][nbin],deltaP[8][nbin];
    double gpfitlist[8][nbin],gdeltaP[8][nbin];
    rep(ite,8)rep(bin,nbin){
        pfitlist[ite][bin] = DINF;
        gpfitlist[ite][bin] = DINF;
    }
    ReadFile(i,pfitlist,deltaP);
    ReadGfile(i,gpfitlist,gdeltaP);
    TGraph* rgraph = new TGraph;
    int rbin = 0;
    prep(bin,sb,fb){
        if(pfitlist[j][bin]==DINF)continue;
        
        /*if(pfitlist[j][bin]/gpfitlist[j][bin]>1.01 || pfitlist[j][bin]/gpfitlist[j][bin]<0.99){
            cout << bin << endl;
            continue;
        }*/
        rgraph -> SetPoint(rbin,Freq(i,bin),pfitlist[j][bin]/gpfitlist[j][bin]);
        rbin++;
    }
    axrange axratio = {213.8+2*i,216.2+2*i,0,5,0,1,";Freq[GHz];Ratio"};

    st.Graph(rgraph,axratio);
    rgraph -> Draw("AP");
}
void DrawLogGraphExample() {
    //TCanvas *canvas = new TCanvas("canvas", "Canvas Title", 800, 600);
    //canvas->SetLogx();

    TGraph *graph = new TGraph();
    graph->SetPoint(0, 1, 1);
    graph->SetPoint(1, 10, 10);
    graph->SetPoint(2, 100, 100);

    graph->Draw("APL");
    //canvas->Draw();
}

//これがメイン関数
void savefitres(){
    
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.14,0.2,0.1);
    c1 -> SetLogx();//同じく横軸もログスケールにしたい場合はここを外す
    c1 -> SetLogy();//ログスケールのON/OFFはここだけで！
    
    //まずはそれぞれ別々の情報が詰められているか確認する
    //P_dP();
    // DrawPower(phist);
    // st.Hist(phist);
    // phist -> Draw();
    // phist -> Fit("gaus");
    DrawLimit();
    //c1 -> Draw();
    //DrawRatio(1,0);
    //DrawLogGraphExample();
    //フィット結果そのものの比較,エラーとのコンシステンシー
    /*for(int fn=1;fn<25;fn++){
        axtest = {213.8+2*fn,216.2+2*fn,0,2,0,1,";Freq[GHz];ratio"};
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
        for(int j=0;j<4;j++){
            prep(p,1,3){
                double testlist[nbin];
                double testdeltaP[nbin];
                //double gtestlist[nbin];
                //double gtestdeltaP[nbin];
                rep(bin,nbin)testlist[bin] = DINF;
                
                TH1D* scalehist = new TH1D("scalehist",";P_{fit}/#Delta P_{fit};Count",100,-10,10);
                TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
                c1 -> SetMargin(0.14,0.11,0.2,0.1);
                //fitterを毎回回さなくてもいいように確定版のデータでなくてもいいのでrootファイルを作成して保存しておきたい
                //gloGetDPfit(fn,j,p,gtestlist,gtestdeltaP);
                GetDPfit(fn,j,p,testlist,testdeltaP);
                prep(bin,sb,fb){
                    pfitlist[2*j+(p-1)][bin] = testlist[bin];
                    deltaP[2*j+(p-1)][bin] = testdeltaP[bin];
                }
                //このfor文内でコメントアウトするときにはここを使おう！
            }
        }
        filesystem::current_path("/Users/oginokyousuke/data/test/");
        string frname = "pfitres"+to_string(fn)+".root";
        string fgrname = "pfitres"+to_string(fn)+"global.root";
        TFile* frfile = new TFile(frname.c_str(),"recreate");
        rep(j,8){
            string ftname = "rtree"+to_string(j);
            TTree* frtree = new TTree(ftname.c_str(),ftname.c_str());
            int fbinF;
            double pfitF,delpfitF;
            frtree -> Branch("bin",&fbinF,"bin/I");
            frtree -> Branch("pfit",&pfitF,"pfit/D");
            frtree -> Branch("delpfit",&delpfitF,"delpfit/D");
            prep(bin,sb,fb){
                if(pfitlist[j][bin]==DINF)continue;
                fbinF = bin;
                pfitF = pfitlist[j][bin];
                delpfitF = deltaP[j][bin];
                frtree -> Fill();
            }
            frtree -> Write();
        }
        
        frfile -> Write();
        frfile -> Close();
    }*/
    return;
}