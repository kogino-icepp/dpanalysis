#include <iostream>
#include <queue>
#include "../headers/fitter.h"
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
const double DINF=1e9;
const double DeltaP = 0;

string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
string savedirp = "/Users/oginokyousuke/data/peak_data/";
string saveexe = "/Users/oginokyousuke/data/search_exe/";
Fitter ft;
axrange axscale = {0,1,0,1,0,1,"test_data;xscale;yscale"};
vector<int> sbin = {0,512,-512,-1024};
double psigma[4] = {0.18068,0.183154,0.194241,0.186464};
double psigma2[4] = {0.160424,0.161409,0.172565,0.166457};//ベースラインとの合わせ技でフィットした場合の温度幅
double wsigma[4] = {0.066,0.066,0.066,0.066};

//フィットに使う諸々の関数群
Double_t chiF_free(double x,double p0,double k,double p1){
    return p0*TMath::Gamma(k/2,x/(2*p1));
}
/// @brief 
/// @param x 
/// @param p0 
/// @param p1 
/// @param k 
/// @param bin 
/// @return 
Double_t chiF_freefit(double x,double p0,double p1,double k,double bin){
    return (chiF_free((x+bin/2),p0,k,p1)-chiF_free((x-bin/2),p0,k,p1));
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
    else return 0;
}
void PrintEntryInfo(const char* filename, const char* treeName, Int_t numEntriesToShow) {
    filesystem::path path=filesystem::current_path();
    filesystem::current_path(saveexe);
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
    // Get the number of entries in the tree
    Int_t numEntries = tree->GetEntries();
    // Determine the number of entries to display
    Int_t numToShow = (numEntriesToShow <= 0) ? numEntries : std::min(numEntriesToShow, numEntries);
    // Loop over the entries and display information
    for (Int_t entry = 0; entry < numToShow; ++entry) {
        tree->GetEntry(entry);
        // Print entry information here
        // For example, print the values of specific branches
        Double_t a, b, c, chi,freq;
        Int_t bin;
        tree->SetBranchAddress("a", &a);
        tree->SetBranchAddress("b", &b);
        tree->SetBranchAddress("c", &c);
        tree->SetBranchAddress("chi", &chi);
        tree->SetBranchAddress("freq",&freq);
        tree->SetBranchAddress("bin",&bin);
        //tree->SetBranchAddress("ndf", &ndf);
        // Print the values for this entry
        cout << "==================" << endl;
        cout << "Entry " << entry << ": a = " << a << ", b = " << b << ", c = " << c<< endl;
        cout << "Freq : " << freq << " bin :" << bin << endl;
        cout << "#chi : " << chi << endl; 
    }
    // Close the file
    file->Close();
    return;
}
/// @brief 
/// @param i 
/// @param j 
/// @param p 
/// @param prec 
void GetBasicData(int i,int j,int p,TGraph*prec){
    filesystem::path path=filesystem::current_path();
    string cdir=dir+"band"+to_string(i);
    filesystem::current_path(cdir);
    Setting st;
    st.dot_size=0.8;
    st.markerstyle=20;
    st.color = kGreen;
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
                if(num==p && body==0){
                    for(int bin=0;bin<nbin;bin++){
                        tree -> GetEntry(bin);
                        Freq[bin]=freq;
                        cold[bin]=power;
                    }
                }
                else if(num==p && body==1){
                    for(int bin=0;bin<nbin;bin++){
                        tree->GetEntry(bin);
                        hot[bin]=power;
                    }
                }
                else if(num==p && body==2){
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
/// @brief 
/// @param offset 
/// @param i 
/// @param fn 
/// @param hist 
void ChiCheck(int offset,int fn,TH1D*hist){
    filesystem::path path=filesystem::current_path();
    filesystem::current_path(saveexe);
    Setting st;
    string fnhon = "all_baseline"+to_string(fn)+".root";
    string fnkai = "all_basekai"+to_string(fn)+".root";
    string fnura = "all_baseura"+to_string(fn)+".root";
    string tname = "tree"+to_string(offset);
    const char* filename = fnkai.c_str();
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
    TH1D* chihist = new TH1D("chihist","chihist;Chi2/NDF;Count",100,0,5);
    rep(i,numEntries){
        tree -> GetEntry(i);
        Double_t chi;
        tree->SetBranchAddress("chi", &chi);
        hist -> Fill(chi);
        cout << chi << endl;
    }
    return;
}

void WhiteCheck(int offset,int i,int j,int p,TH1D*hist){
    filesystem::path path=filesystem::current_path();
    filesystem::current_path(saveexe);
    /*TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);*/
    Setting st;
    st.dot_size=0.8;
    st.markerstyle=20;
    st.color = kGreen;
    string fntest = "all_baseline_test.root";
    string fnhon = "all_baseline"+to_string(j)+".root";
    string fnkai = "all_basekai"+to_string(j)+".root";
    string tname = "tree"+to_string(offset);
    const char* filename = fnkai.c_str();
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
        
    }
    TGraph* pgraph = new TGraph;
    GetBasicData(i,j,p,pgraph);
    //TH1D* whist = new TH1D("whist","WhiteNoise;white_noise[K];Count",100,-1,1);
    for(int bin=sb;bin<fb;bin+=dbin){
        bin+=offset;
        if(vpfreq[bin]==DINF)continue;
        TF1* quadf = new TF1("quadf","[0]*(x-[1])*(x-[1])+[2]",0,1);
        TGraphErrors* wgraph = new TGraphErrors;
        quadf -> SetParameter(0,vparas[0][bin]);
        quadf -> SetParameter(1,vparas[1][bin]);
        quadf -> SetParameter(2,vparas[2][bin]);
        TGraphErrors* spgraph = new TGraphErrors;
        double xmin,ymin,xscale,yscale;
        //ft.make_scale2(spgraph,pgraph,bin-sb,xmin,ymin,xscale,yscale);
        //ft.rescale_para(vparas[0][bin],vparas[1][bin],vparas[2][bin],xmin,ymin,xscale,yscale);
        st.GraphErrors(spgraph,axscale);
        spgraph -> Draw("AP");
        quadf -> Draw("same");
        // /spgraph -> Fit(quadf,"M","",0,1);
        rep(k,dbin){
            double xValue = spgraph -> GetPointX(k);
            double yValue = quadf -> Eval(xValue);
            double yTrue = spgraph -> GetPointY(k);
            hist -> Fill((yValue-yTrue)*yscale);
            if(abs((yValue-yTrue)*yscale)>1)hist -> Fill(0.99);
        }
    }
    /*st.Hist(whist);
    c1 -> SetLogy();
    whist -> Draw();
    whist -> Fit("gaus","MQ","",-1,1);*/
    return;
}
/// @brief 
/// @param offset 
/// @param i 
/// @param j 
/// @param p 
/// @param hist 
/// @param que 
void PeakFit(int offset,int i,int j,int p,TH1D* hist,queue<int>&que){
    filesystem::path path=filesystem::current_path();
    filesystem::current_path(saveexe);
    /*TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);*/
    Setting st;
    st.dot_size=0.8;
    st.markerstyle=20;
    st.color = kGreen;
    st.lcolor = kGreen;
    string fntest = "all_baseline_test.root";
    string fnhon = "all_baseline"+to_string(j)+".root";
    string fnkai = "all_basekai"+to_string(j)+".root";
    string fura = "all_baseura"+to_string(j)+".root";
    string tname = "tree"+to_string(offset);
    const char* filename = fnkai.c_str();
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
        
    }
    TGraph* pgraph = new TGraph;
    double ifmin = 213.9+2*i;
    double ifmax = 216.1+2*i;
    GetBasicData(i,j,p,pgraph);
    //axrange axraw = {ifmin,ifmax,0,100,0,1};
    //st.Graph(pgraph,axraw);
    //pgraph -> Draw("AP");
    //サーチするビンはベースライン基点からさらに15binあと
    for(int bin=8645;bin<916;bin+=dbin){
        TString histName = Form("whist_%d", bin);
        //cout << histName << endl;
        TH1D* whist = new TH1D(histName.Data(),"whist;white_noise[K];Count",100,-1,1);
        bin+=offset;
        double s1 = pgraph -> GetPointX(bin-sb);
        double s2 = pgraph -> GetPointX(bin-sb+30);
        double mfreq = pgraph -> GetPointX(bin-sb+11);
        double sfreq = min(s1,s2);
        double ffreq = max(s1,s2);
        TF1* quadf = new TF1("quadf","[0]*(x-[1])*(x-[1])+[2]",0,1);
        TF1* peakf = new TF1("peakf","F_sig2(x,[0],[1],0.5)",sfreq,ffreq);
        TGraphErrors* wgraph = new TGraphErrors;
        
        if(vpfreq[bin]==DINF){
            //cout << "Pass " << endl;
            continue;
        }
        //cout << vparas[0][bin] << " " << vparas[1][bin] << " " << vparas[2][bin] << endl;
        quadf -> SetParameter(0,vparas[0][bin]);
        quadf -> SetParameter(1,vparas[1][bin]);
        quadf -> SetParameter(2,vparas[2][bin]);
        TGraphErrors* spgraph = new TGraphErrors;
        double yscale;
        ft.make_scale(spgraph,pgraph,bin-sb,yscale);
        double ave = 0;
        rep(k,dbin){
            double xValue = spgraph -> GetPointX(k);
            double yValue = quadf -> Eval(xValue);
            double yTrue = spgraph -> GetPointY(k);
            if(k<10 ||  k>=20){
                double ys = (yValue-yTrue)*yscale;
                //cout << ys << endl;
                ave+=ys*ys;
            }
            double freqbin = pgraph -> GetPointX(bin-sb+k);
            wgraph -> SetPoint(k,freqbin,(yValue-yTrue)*yscale);
        }
        ave/=20;
        double sigma = sqrt(ave);
        //cout << sigma << endl;
        //hist -> Fill(sigma);
        rep(k,dbin)wgraph -> SetPointError(k,0,sigma);
        //一点のみのフィットを行う
        peakf -> FixParameter(0,mfreq);
        peakf -> SetParameter(1,0.5);
        rep(ite,5)wgraph -> Fit(peakf,"QE","",sfreq,ffreq);
        double pout = peakf -> GetParameter(1);
        double chi = peakf -> GetChisquare();
        int ndf = peakf -> GetNDF();
        //cout << bin << " " << chi/ndf << " " << pout/psigma[j] << endl;
        //if(abs(pout)>1)hist -> Fill(0.99);
        hist -> Fill(pout);
        if(abs(pout)>4*psigma2[j])cout << bin << " " << pout/psigma[j] << endl;
        //cout << bin << " -> " << sigma << endl;
        //cout << "pout : " << pout << " && chi/ndf : " << chi/ndf << endl;
        axrange axw = {sfreq,ffreq,-0.5,0.5,0,1,"fit_example;Freq[GHz];WhiteNoise[K]"};
        st.GraphErrors(wgraph,axw);
        wgraph -> Draw("AP");
        peakf -> Draw("same");
        //ここでデータを格納
        //que.push(pout);
        /*if(abs(pout)>4*psigma[j]){
            cout << bin << " " << pout/psigma[j] << endl;
            que.push(bin);
        }
        if(abs(pout*2*kb*df)<4*pow(10,-18))hist -> Fill(pout*2*kb*df);
        else hist -> Fill(3.99*pow(10,-18));*/
        bin-=offset;
        //st.GraphErrors(wgraph,axw);
        //wgraph -> Draw("AP");
        //peakf -> Draw("same");
        //whist -> Delete();
    }
    
    file -> Close();
    return;
}
//ホワイトノイズ的にきちんと詰められていそう→これに対してピークサーチをかける？
//壊れたら嫌なので信号+ベースラインでのフィットはこっちでやる
void PeakFit2(int offset,int i,int j,int p,TH1D* hist,queue<int> &que){
    filesystem::path path=filesystem::current_path();
    filesystem::current_path(saveexe);
    /*TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);*/
    Setting st;
    st.dot_size=0.8;
    st.markerstyle=20;
    st.color = kGreen;
    st.lcolor = kGreen;
    string fntest = "all_baseline_test.root";
    string fnhon = "all_baseline"+to_string(j)+".root";
    string fnkai = "all_basekai"+to_string(j)+".root";
    string fura = "all_baseura"+to_string(j)+".root";
    string tname = "tree"+to_string(offset);
    const char* filename = fnkai.c_str();
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
        
    }
    TGraph* pgraph = new TGraph;
    double ifmin = 213.9+2*i;
    double ifmax = 216.1+2*i;
    GetBasicData(i,j,p,pgraph);
    //axrange axraw = {ifmin,ifmax,0,100,0,1};
    //st.Graph(pgraph,axraw);
    //pgraph -> Draw("AP");
    //サーチするビンはベースライン基点からさらに15binあと
    for(int bin=23645;bin<23646;bin+=dbin){
        bin+=offset;
        if(vpfreq[bin]==DINF){
            //cout << "Pass " << endl;
            continue;
        }
        //cout << vpfreq[bin] << endl;
        double s1 = pgraph -> GetPointX(bin-sb);
        double s2 = pgraph -> GetPointX(bin-sb+30);
        double mfreq = pgraph -> GetPointX(bin-sb+11);
        double sfreq = min(s1,s2);
        double ffreq = max(s1,s2);
        TF1* peakquad = new TF1("peakquad","[0]*(x-[1])*(x-[1])+[2]+F_sig2(x,[3],[4],0.5)");
        TGraphErrors* spgraph = new TGraphErrors;
        double yMin,yscale;
        ft.make_scale2(pgraph,bin-sb,yMin,yscale);
        ft.rescale_para(vparas[0][bin],vparas[1][bin],vparas[2][bin],sfreq,yMin,ffreq-sfreq,yscale,peakquad);
        axrange axtest = {sfreq,ffreq,30,80,0,1,"SignalFit;Freq[GHz];Prec[K]"};
        //一からプロットし直してエラーをグローバルにつける
        TGraphErrors* fitgraph = new TGraphErrors;
        rep(k,dbin){
            double x = pgraph -> GetPointX(bin-sb+k);
            double y = pgraph -> GetPointY(bin-sb+k);
            fitgraph -> SetPoint(k,x,y);
            fitgraph -> SetPointError(k,0,0.066);
        }
        peakquad -> FixParameter(3,mfreq);
        peakquad -> SetParameter(4,0.1);
        rep(ite,5)fitgraph -> Fit(peakquad,"QE","",sfreq,ffreq);
        double pout = peakquad -> GetParameter(4);
        //pout *= 2*kb*df;
        //cout << pout << endl;
        st.GraphErrors(fitgraph,axtest);
        fitgraph -> Draw("AP");
        peakquad -> Draw("same");
        double chi = peakquad -> GetChisquare();
        double ndf = peakquad -> GetNDF();
        cout << pout/psigma2[j] << " : " << chi/ndf << endl;
        
        //if(chi/ndf<1 &&  chi/ndf >0.7 && pout/psigma2[j]>2)que.push(bin);
        if(chi/ndf<6)hist -> Fill(pout/psigma2[j]);
        else que.push(bin);
        bin-=offset;
    }
    file -> Close();
    return;
}
void GetDPfit(){
    Setting st;
    st.dot_size=0.8;
    st.markerstyle=20;
    st.color = kGreen;
    st.lcolor = kBlue;
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    //どの単位でデータを持ってくる？→どうせ4データで平均は取りたい,メモリの問題は起こってから考える
    //各データ、同一セットアップ2回の測定に対して1ビンずつずらしたベースラインのデータを保持、
    double Pfit[nbin][8];
    double DPfit[8];
    double Freq[nbin];
    //とりあえず初期化
    rep(i,nbin)rep(j,8)Pfit[i][j] = 0;
    prep(i,0,3){
        filesystem::path path=filesystem::current_path();
        filesystem::current_path(saveexe);
        string fnname = "allbinbase5_"+to_string(i)+".root";
        const char* fileName = fnname.c_str();
        TFile*file = TFile::Open(fileName);
        if (!file) {
            cout << "Error opening file " << fileName << endl;
            return;
        }
        prep(j,1,3){
            string tname = "test_tree"+to_string(j);
            const char* treeName = tname.c_str();
            TTree*tree = dynamic_cast<TTree*>(file->Get(treeName));
            if (!tree) {
                cerr << "Error getting tree " << treeName << " from file " << fileName << std::endl;
                file->Close();
                return;
            }
            Int_t numEntries = tree->GetEntries();
            vector<vector<double>> vparas(3,vector<double>(nbin));
            vector<double> vpfreq(nbin,DINF);
            vector<int> vbin(nbin);
            rep(k,numEntries){
                tree -> GetEntry(k);
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
                if(i==0 && j==1)Freq[bin]=freq;
            }
            //freqをdata0のものに合わせたい→+sbin[i]したパートに格納すればオッケー？
            TGraph* pgraph = new TGraph;
            GetBasicData(5,i,(2-j)%2,pgraph);
            prep(bin,sb,fb){
                if(vpfreq[bin]==DINF)continue;
                double sfreq = pgraph -> GetPointX(bin-sb);
                double ffreq = pgraph -> GetPointX(bin-sb+30);
                double mfreq = pgraph -> GetPointX(bin-sb+11);
                TF1* peakquad = new TF1("peakquad","[0]*(x-[1])*(x-[1])+[2]+F_sig2(x,[3],[4],0.5)");
                TGraphErrors* spgraph = new TGraphErrors;
                double yMin,yscale;
                ft.make_scale2(pgraph,bin-sb,yMin,yscale);
                ft.rescale_para(vparas[0][bin],vparas[1][bin],vparas[2][bin],sfreq,yMin,ffreq-sfreq,yscale,peakquad);
                TGraphErrors* fitgraph = new TGraphErrors;
                rep(k,dbin){
                    double x = pgraph -> GetPointX(bin-sb+k);
                    double y = pgraph -> GetPointY(bin-sb+k);
                    fitgraph -> SetPoint(k,x,y);
                    fitgraph -> SetPointError(k,0,0.066);
                }
                peakquad -> FixParameter(3,mfreq);
                peakquad -> SetParameter(4,0.1);
                rep(ite,5)fitgraph -> Fit(peakquad,"QE","",sfreq,ffreq);
                double pout = peakquad -> GetParameter(4);
                Pfit[bin+sbin[i]][i*2+j] = pout;
            }
        }
        /*
        以下でやりたいこと
        これのデバッグやりたくなさすぎ、どうすればええんや→関数をもう少し細分化してみるのはアリ
        1.　各点最大8データあるのでそれぞれピークフィットまでやって結果を格納
        2. 8データで統計とって平均値をΔP_fitに変換して格納
        3. 適宜グラフなどに起こしてreturn
        */
    }
    TGraph* gplim = new TGraph;
    axrange axpl = {224,226,0,pow(10,-14),0,1,"Plim;Freq[GHz];#{Delta}P_{95#{%}C.L.}"};
    prep(bin,0,nbin){
        int nstat = 0;
        double plim = 0;
        rep(i,8){
            if(Pfit[bin][i]==0)continue;
            else{
                nstat++;
                plim += Pfit[bin][i];
            }
        }
        if(nstat<2)continue;
        plim /= nstat;
        plim += DeltaP;
        plim *= 2*kb*df;
        gplim -> SetPoint(bin,Freq[bin],plim);
    }
    c1 -> SetLogy();
    st.Graph(gplim,axpl);
    gplim -> Draw("AL");
    return;
}
//これがメイン関数
/*
次の要件定義
周波数ごとにP_fit+1.96ΔPを導出して95%棄却区間を求めたい
1. まずは4回の測定で周波数を合わせる(sbinを各自引いたものを並べる的な)
2. 各々データは持っているので(最悪TTree作って別途引っ張ってくるか)ピークフィットしてP_fitの平均と分散をとる
3. +1.96ΔP(95%に相当)してPowerレベルでのlimit出してTGraphにプロット
これなんか要らなさそうな雰囲気出てたよね？
*/
void peaksearch_mk2(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    const char* filename = "all_baseura0.root";
    const char* treeName = "tree0";
    //まずはそれぞれ別々の情報が詰められているか確認する
    //TH1D* whitehist = new TH1D("whitehist","WhiteNoise;white_noise;Count",100,-1,1);
    //TH1D* peakhist = new TH1D("peakhist","peakhist;P[W];Count",100,-5*pow(10,-18),5*pow(10,-18));
    TH1D* normhist = new TH1D("normhist","NormHist;Sigma;Count",100,-10,10);
    //TH1D* wpeakhist = new TH1D("wpeakhist","P_fit;P_fit[W];Count",100,-4*pow(10,-18),4*pow(10,-18));
    TH1D* chihist = new TH1D("chihist","chihist;Chi2/NDF;Count",100,0,10);
    //TF1* gausfit = new TF1("gausfit","gaus",-1,1);
    //TF1* chifit = new TF1("chifit","chiF_freefit(x,[0],[1],25,0.1)");
    //PrintEntryInfo(filename,treeName,10);
    queue<int> que;
    //GetDPfit();
    for(int fn=3;fn<4;fn++){
        for(int offset=6;offset<7;offset++){
            //PrintEntryInfo(filename,treeName,10);
            //WhiteCheck(offset,5,fn,1,whitehist);
            PeakFit2(offset,5,fn,1,normhist,que);
            //ChiCheck(offset,fn,whitehist);
        }
    }
    /*while(!que.empty()){
        auto v = que.front();que.pop();
        cout << v << endl;
    }
    c1 -> SetLogy();
    st.Hist(normhist);
    normhist -> Draw();
    normhist -> Fit("gaus");
    //詰められたホワイトノイズやχ^2/ndfが正常かどうかを確認するプログラム*/
    return;
}