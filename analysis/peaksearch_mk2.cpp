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

string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
string savedirp = "/Users/oginokyousuke/data/peak_data/";
string saveexe = "/Users/oginokyousuke/data/search_exe/";
Fitter ft;
axrange axscale = {0,1,0,1,0,1,"test_data;xscale;yscale"};
vector<int> sbin = {0,512,-512,-1024};
double psigma[4] = {0.17664,0.179145,0.190872,0.11982};
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
    }
    // Close the file
    file->Close();
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
}
/// @brief 
/// @param offset 
/// @param i 
/// @param fn 
void ChiCheck(int offset,int i,int fn,TH1D*hist){
    filesystem::path path=filesystem::current_path();
    filesystem::current_path(saveexe);
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    Setting st;
    string fnhon = "all_baseline"+to_string(fn)+".root";
    string fnkai = "all_basekai"+to_string(fn)+".root";
    const char* filename = fnkai.c_str();
    const char* treeName = ("tree"+to_string(offset)).c_str();
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
        chihist -> Fill(chi);
    }
    TF1* chifit = new TF1("chifit","chiF_freefit(x,[0],[1],17,0.05)");
    st.Hist(chihist);
    double entnum = chihist -> GetEntries();
    int maxbin = chihist -> GetMaximumBin();
    double a =0.02;
    chifit -> FixParameter(0,entnum);
    chifit -> SetParameter(1,a);
    c1 -> SetLogy();
    chihist -> Draw();
    chihist -> Fit(chifit,"MQE","",0,5);
    chifit -> Draw("same");
    double chires = chifit -> GetChisquare();
    int ndf = chifit -> GetNDF();
    cout << fn << " " << offset << endl;
    cout << "chi2/ndf : " << chires/ndf << endl;
    //delete[] filename;
    //delete[] treeName;
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
    const char* filename = fnkai.c_str();
    const char* treeName = ("tree"+to_string(offset)).c_str();
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
        double yscale;
        ft.make_scale(spgraph,pgraph,bin-sb,yscale);
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
}
/// @brief 
/// @param offset 
/// @param i 
/// @param j 
/// @param p 
/// @param hist 
/// @param que 
void PeakFit(int offset,int i,int j,int p,TH1D* hist){
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
    const char* filename = fnkai.c_str();
    const char* treeName = ("tree"+to_string(offset)).c_str();
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
    for(int bin=sb;bin<fb;bin+=dbin){
        TString histName = Form("whist_%d", bin);
        //cout << histName << endl;
        TH1D* whist = new TH1D(histName.Data(),"whist;white_noise[K];Count",100,-1,1);
        bin+=offset;
        double s1 = pgraph -> GetPointX(bin-sb-10+15);
        double s2 = pgraph -> GetPointX(bin-sb+10+15);
        double mfreq = pgraph -> GetPointX(bin-sb+15);
        double sfreq = min(s1,s2);
        double ffreq = max(s1,s2);
        TF1* quadf = new TF1("quadf","[0]*(x-[1])*(x-[1])+[2]",0,1);
        TF1* peakf = new TF1("peakf","F_sig2(x,[0],[1],0.5)",sfreq,ffreq);
        TGraphErrors* wgraph = new TGraphErrors;
        
        if(vpfreq[bin]==DINF)continue;
        //cout << vparas[0][bin] << " " << vparas[1][bin] << " " << vparas[2][bin] << endl;
        quadf -> SetParameter(0,vparas[0][bin]);
        quadf -> SetParameter(1,vparas[1][bin]);
        quadf -> SetParameter(2,vparas[2][bin]);
        TGraphErrors* spgraph = new TGraphErrors;
        double yscale;
        ft.make_scale(spgraph,pgraph,bin-sb,yscale);
        st.GraphErrors(spgraph,axscale);
        spgraph -> Draw("AP");
        quadf -> Draw("same");
        // /spgraph -> Fit(quadf,"M","",0,1);
        rep(k,dbin){
            double xValue = spgraph -> GetPointX(k);
            double yValue = quadf -> Eval(xValue);
            double yTrue = spgraph -> GetPointY(k);
            if(k<10 ||  k>=20)whist -> Fill((yValue-yTrue)*yscale);
            double freqbin = pgraph -> GetPointX(bin-sb+k);
            wgraph -> SetPoint(k,freqbin,(yValue-yTrue)*yscale);
        }
        //gausfit -> Draw();
        double xMin = whist->GetXaxis()->GetBinLowEdge(whist->GetMinimumBin());
        double xMax = whist->GetXaxis()->GetBinUpEdge(whist->GetMaximumBin());
        TF1* gausfit = new TF1("gausfit", "gaus",-1,1);
        gausfit -> SetParameter(0,20);
        gausfit -> SetParLimits(1,-0.2,0.2);
        gausfit -> SetParameter(2,0.01);
        
        rep(ite,5)whist -> Fit(gausfit,"MQR","",-1,1);
        whist -> Draw();
        double sigma = gausfit -> GetParameter(2);
        rep(k,dbin)wgraph -> SetPointError(k,0,sigma);
        //一点のみのフィットを行う
        peakf -> FixParameter(0,mfreq);
        peakf -> SetParameter(1,0.5);
        rep(ite,5)wgraph -> Fit(peakf,"MQ","",sfreq,ffreq);
        double pout = peakf -> GetParameter(1);
        double chi = peakf -> GetChisquare();
        int ndf = peakf -> GetNDF();
        if(abs(pout)>1)hist -> Fill(0.99);
        hist -> Fill(pout);
        //if(abs(pout)>4*psigma[j])cout << bin << " " << pout/psigma[j] << endl;
        cout << bin << " -> " << sigma << endl;
        //cout << "pout : " << pout << " && chi/ndf : " << chi/ndf << endl;
        axrange axw = {sfreq,ffreq,-2,1,0,1,"fit_example;Freq[GHz];WhiteNoise[K]"};
        st.GraphErrors(wgraph,axw);
        wgraph -> Draw("AP");
        peakf -> Draw("same");
        //ここでデータを格納
        //que.push(pout);
        if(abs(pout)>4*psigma[j])cout << bin << " " << pout/psigma[j] << endl;
        bin-=offset;
        //st.GraphErrors(wgraph,axw);
        //wgraph -> Draw("AP");
        //peakf -> Draw("same");
        //whist -> Delete();
    }
    
    file -> Close();
    
}
//ホワイトノイズ的にきちんと詰められていそう→これに対してピークサーチをかける？

//これがメイン関数
void peaksearch_mk2(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    const char* filename = "all_baseline0.root";
    const char* fntest = "all_baseline_test.root";
    //まずはそれぞれ別々の情報が詰められているか確認する
    TH1D* whitehist = new TH1D("whitehist","WhiteNoise;white_noise;Count",100,-1,1);
    TH1D* peakhist = new TH1D("peakhist","peakhist;P[kHz*K];Count",100,-1,1);
    //TH1D* chihist = new TH1D("chihist","chihist;Chi2/NDF;Count",100,0,10);
    //TF1* gausfit = new TF1("gausfit","gaus",-1,1);
    //TF1* chifit = new TF1("chifit","chiF_freefit(x,[0],[1],18,0.1)");
    queue<double> que;
    queue<pair<int,double>> pque;

    for(int fn=2;fn<3;fn++){
        for(int offset=0;offset<1;offset++){
            //PrintEntryInfo(filename,tname,10);
            //WhiteCheck(offset,5,fn,1,whitehist);
            PeakFit(offset,5,fn,1,peakhist);
            //ChiCheck(offset,5,fn);
        }
    }
    st.Hist(peakhist);
    peakhist -> Draw();
    c1 -> SetLogy();
    peakhist -> Fit("gaus");
    /*st.Hist(chihist);
    int entnum = chihist -> GetEntries();
    int maxbin = chihist -> GetMaximumBin();
    double a = maxbin*0.1/(18-2);
    chifit -> FixParameter(0,entnum);
    chifit -> SetParameter(1,a);
    c1 -> SetLogy();
    chihist -> Draw();
    chihist -> Fit(chifit);
    chifit -> Draw("same");
    //hist -> Fit(gausfit,"MQE","",-1,1);
    //double normsig = gausfit -> GetParameter("Sigma");
    //cout << normsig << endl;
    /*TH1D* normhist = new TH1D("normhist","normhist;Sigma;Count",100,-5,5);
    while(!que.empty()){
        auto v = que.front();que.pop();
        //cout << v << endl;
        v /= normsig;
        normhist -> Fill(v);
        if(abs(v)>5)cout << v << endl;
    }
    c1 -> SetLogy();
    st.Hist(normhist);
    normhist -> Draw();
    normhist -> Fit("gaus");*/
    //詰められたホワイトノイズやχ^2/ndfが正常かどうかを確認するプログラム
    
}