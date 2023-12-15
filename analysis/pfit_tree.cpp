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
void ReadBaseInfo(int i,int j,int p,double (&vparas)[3][nbin]){
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
    }
    file -> Close();
}
void MakeTree(){
    TFile* wfile = new TFile;

}
void GetDPfit(int i,int j,int p,string wfileName){
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

    int xfft = XFFT(i);
    /*
    要件定義
    最終結果ならびに結果の精査に必要な情報を保存するためにrootファイルを作成して保存する機能を追加する
    ・保存したいパラメータ：bin,周波数,Pfit,dPfit,chi/ndf
    ・注意事項：rootファイルを開いた状態で別のアレを開くとアレするらしいので競合しないように適宜closeしたり順番を工夫する
    */

    
    TFile* wfile = new TFile(wfileName.c_str(),"update");
    //string wtreeName = "wtree"+to_string(2*j+(p-1));
    //TTree* wtree = new TTree(wtreeName.c_str(),wtreeName.c_str());
    //double testpara = 1;
    //wtree -> Branch("testpara",&testpara,"testpara/D");
    //wtree -> Fill();
    //wtree -> Branch("ffreq",&freqF,"ffreq/D");
    //wtree -> Branch("fchi",&chiF,"fchi/D");
    //wtree -> Branch("Pfit",&PfitF,"Pfit/D");
    //wtree -> Branch("dPfit",&dPfitF,"dPfit/D");
    //wtree -> Branch("fbin",&binF,"fbin/I");

    /*for(int bin=sb;bin<fb;bin++){
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
        sigma /= 17;
        sigma = sqrt(sigma);
        rep(k,dbin){
            fitgraph -> SetPointError(k,0,sigma);
            spgraph -> SetPointError(k,0,sigma/yscale);
        }
        st.GraphErrors(spgraph,axscale);
        st.GraphErrors(fitgraph,axtest);
        st.GraphErrors(fitgraph2,axtest);
        double smfreq = spgraph -> GetPointX(bpos);
        peakquad -> FixParameter(3,mfreq);
        scalepeak -> FixParameter(4,sfreq);
        rep(ite,10)spgraph -> Fit(scalepeak,"Q0","",0,1);
        rep(ite,10)spgraph -> Fit(scalepeak,"MQ0","",0,1);
        
        //フィットがある程度収束するまでこれ続ける
        //fitgraph -> Fit(peakquad,"EQ","",sfreq,ffreq);
        spgraph -> Fit(scalepeak,"EQ","",0,1);
        
        double spout = scalepeak -> GetParameter(3);
        double sdpout = scalepeak -> GetParError(3);
        spout *= yscale;
        sdpout *= yscale;
        firsthist -> Fill(spout/(sdpout));
    }
    double fmin = 213.5+i*2;
    double fmax = 216.5+i*2;
    TF1* fgaus = new TF1("fgaus","gaus");
    firsthist -> Fit(fgaus);
    double sigma = fgaus -> GetParameter("Sigma");*/
    
}

//メイン関数
void pfit_tree(){
    string testname = "test1.root";
    double vparas[3][nbin];
    ReadBaseInfo(1,2,1,vparas);
    MakeTree();
}