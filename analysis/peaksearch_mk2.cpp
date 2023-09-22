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
const double Aeff = 0.385;

string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
string savedirp = "/Users/oginokyousuke/data/peak_data/";
string saveexe = "/Users/oginokyousuke/data/search_exe/";
Fitter ft;
axrange axscale = {0,1,0,1,0,1,"test_data;xscale;yscale"};
vector<int> sbin = {0,512,-512,-1024};
double psigma[4] = {0.18068,0.183154,0.194241,0.186464};
double psigma2[4] = {0.162231,0.161409,0.172565,0.166457};//ベースラインとの合わせ技でフィット　した場合の温度幅　　ここで作ったσでnormalizeすると原義1に本当になるのか
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
double CoupConst2(double p,double dp){
    if(p<0)return 2*1.3458*4.5*pow(10,-14)*sqrt(dp*1.96*pow(10,23))*sqrt(1/Aeff);
    else return 2*1.3458*4.5*pow(10,-14)*sqrt((dp*1.96+p)*pow(10,23))*sqrt(1/Aeff);
}
//ノーマライズされていることは前提としてそのカイ二乗値のp値をお手軽に返してくれる関数
double PValue(int k,double x){
    return 1-TMath::Gamma(k/2,(k/2)*x);
}
void PrintEntryInfo(const char* filename, const char* treeName, Int_t numEntriesToShow) {
    filesystem::path path=filesystem::current_path();
    filesystem::current_path(saveexe);
    TFile*file = TFile::Open(filename);
    if (!file) {
        cerr << "Error opening file " << filename << endl;
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
        cout << "chi : " << chi << endl; 
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
        //cout << chi << endl;
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
//ホワイトノイズ的にきちんと詰められていそう→これに対してピークサーチをかける？
//壊れたら嫌なので信号+ベースラインでのフィットはこっちでやる
void PeakFit2(int offset,int i,int j,int p,TH1D* hist,vector<int>&vec){
    filesystem::path path=filesystem::current_path();
    filesystem::current_path(saveexe);
    /*TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);*/
    Setting st;
    st.dot_size=0.8;
    st.markerstyle=20;
    st.color = kGreen;
    st.lcolor = kGreen;
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
    

    for(int bin=23015;bin<23016;bin+=dbin){
        bin+=offset;
        if(vpfreq[bin]==DINF){
            //cout << "Pass " << endl;
            bin-=offset;
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
        axrange axtest = {sfreq,ffreq,0,100,0,1,"SignalFit;Freq[GHz];Prec[K]"};
        //一からプロットし直してエラーをグローバルにつける
        TGraphErrors* fitgraph = new TGraphErrors;
        double sigma = 0;
        double a = peakquad -> GetParameter(0);
        double b = peakquad -> GetParameter(1);
        double c = peakquad -> GetParameter(2);
        rep(k,dbin){
            double x = pgraph -> GetPointX(bin-sb+k);
            double y = pgraph -> GetPointY(bin-sb+k);
            rep(k,dbin)fitgraph -> SetPointError(k,0,0.068);
            fitgraph -> SetPoint(k,x,y);
            if(k<10 || k>=20){
                y -= a*(x-b)*(x-b)+c;
                sigma += y*y;
            }
        }
        sigma /= 17;
        sigma = sqrt(sigma);
        
        //cout << sigma << endl;
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
        if(abs(pout/psigma2[j])>5)vec.push_back(bin);
        cout  << pout/psigma2[j] << " : " << chi/ndf << endl;
        hist -> Fill(chi/ndf);
        
        //else hist -> Fill(9.9);
        
        bin-=offset;
    }
    file -> Close();
    return;
}
//offsetとかいうゴミ概念を除いた解析プログラム、同じように動くので次は4つのデータでもってlimitをつけよう
void GetDPfit(int i,int j,int p,TH1D*hist){
    filesystem::path path=filesystem::current_path();
    filesystem::current_path(saveexe);
    Setting st;
    st.dot_size=0.8;
    st.markerstyle=20;
    st.color = kGreen;
    st.lcolor = kBlue;
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
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
    }
    TGraph* pgraph = new TGraph;
    GetBasicData(i,j,p,pgraph);

    //vectorにとりあえずのフィット結果を詰める
    //どこかにグローバルなエラーと局所的なエラーを両方出すプログラムが欲しい
    vector<double> Pfit(nbin,DINF);
    for(int bin=sb;bin<fb;bin++){
        if(vparas[0][bin]==DINF &&vparas[1][bin]==DINF  &&vparas[2][bin]==DINF)continue;
        //cout << vpfreq[bin] << endl;
        double s1 = pgraph -> GetPointX(bin-sb);
        double s2 = pgraph -> GetPointX(bin-sb+30);
        double mfreq = pgraph -> GetPointX(bin-sb+11);
        double sfreq = min(s1,s2);
        double ffreq = max(s1,s2);
        TF1* peakquad = new TF1("peakquad","[0]*(x-[1])*(x-[1])+[2]+F_sig2(x,[3],[4],0.5)");
        TF1* syokiquad = new TF1("syokiquad","[0]*(x-[1])*(x-[1])+[2]",sfreq,ffreq);
        TGraphErrors* spgraph = new TGraphErrors;
        double yMin,yscale;
        ft.make_scale2(pgraph,bin-sb,yMin,yscale);
        ft.rescale_para(vparas[0][bin],vparas[1][bin],vparas[2][bin],sfreq,yMin,ffreq-sfreq,yscale,peakquad);
        axrange axtest = {sfreq,ffreq,0,100,0,1,"SignalFit;Freq[GHz];Prec[K]"};
        //一からプロットし直してエラーをグローバルにつける
        TGraphErrors* fitgraph = new TGraphErrors;
        
        double a = peakquad -> GetParameter(0);
        double b = peakquad -> GetParameter(1);
        double c = peakquad -> GetParameter(2);
        //cout << a << " " << b << " " << c << endl;
        syokiquad -> SetParameter(0,a);
        syokiquad -> SetParameter(1,b);
        syokiquad -> SetParameter(2,c);
        //ここで解析的にエラーを出す
        double sigma = 0;
        rep(k,dbin){
            double x = pgraph -> GetPointX(bin-sb+k);
            double y = pgraph -> GetPointY(bin-sb+k);
            fitgraph -> SetPointError(k,0,0.068);
            fitgraph -> SetPoint(k,x,y);
            if(k<10 || k>=20){
                sigma += pow(y-(a*(x-b)*(x-b)+c),2);
            }
        }
        sigma /= 17;
        sigma = sqrt(sigma);
        if(sigma<0.5)hist -> Fill(sigma);
        else hist -> Fill(0.48);
        if(sigma>3){
            cout << bin << " " << sigma << endl;
        }
        axrange axfit = {sfreq,ffreq,40,80,0,1,"fitgraph;Freq[GHz];Prec[K]"};
        st.GraphErrors(fitgraph,axfit);
        fitgraph -> Draw("AP");
        syokiquad -> Draw("same");
        /*peakquad -> FixParameter(3,mfreq);
        peakquad -> SetParameter(4,0.1);
        rep(ite,5)fitgraph -> Fit(peakquad,"QE","",sfreq,ffreq);
        
        double pout = peakquad -> GetParameter(4);
        //pout *= 2*kb*df;

        st.GraphErrors(fitgraph,axtest);
        fitgraph -> Draw("AP");
        peakquad -> Draw("same");
        double chi = peakquad -> GetChisquare();
        double ndf = peakquad -> GetNDF();
        if(chi/ndf<3){
            hist -> Fill(pout);
            Pfit[bin+11] = pout;
        }*/
    }
    /*TF1* fgaus = new TF1("fgaus","gaus",-1,1);
    hist -> Fit(fgaus);
    double DPsigma = fgaus -> GetParameter("Sigma");
    DPsigma *= 2*kb*df;
    cout << "#DeltaP_{fit} : " << DPsigma << endl;
    TGraph* Plim = new TGraph;
    int pnum = 0;
    prep(bin,sb,fb){
        if(Pfit[bin]==DINF || vpfreq[bin]==DINF)continue;
        else Pfit[bin] = CoupConst2(Pfit[bin]*2*kb*df,DPsigma);
        Plim -> SetPoint(pnum,vpfreq[bin],Pfit[bin]);
        //vec[bin+sbin[j]]=Pfit[bin];
        pnum++;
        cout << vpfreq[bin] << " " << Pfit[bin] << endl;
    }
    //これは実験的に書いてみた結果 上できちんと格納する
    axrange axlim = {224,226,0,pow(10,-9),0,1,"#chi;Freq[GHz];#chi coupling constant"};
    st.Graph(Plim,axlim);
    c1 -> SetLogy();
    Plim -> Draw("AL");*/
}

//これがメイン関数
void peaksearch_mk2(){
    
    //まずはそれぞれ別々の情報が詰められているか確認する
    
    TH1D* whitehist = new TH1D("whitehist","WhiteNoise;white_noise[K];Count",100,0,0.5);
    TH1D* normhist = new TH1D("normhist","P_{fit}/#DeltaP;P_{fit}/#DeltaP;Count",100,-10,10);
    TF1* chifit = new TF1("chifit","chiF_freefit(x,[0],[1],[2],[3])",0,5);
    //PrintEntryInfo(filename,treeName,10);
    TH1D* chihist = new TH1D("chihist","chihist;#chi^{2}/NDF;Count",100,0,5);
    string whitedir = "/Users/oginokyousuke/data/white_noise";
    //vector<vector<double>> vec(8,vector<double>(nbin,DINF));
    vector<int> excess;
    double deltaP[8][nbin];
    rep(i,8)rep(j,nbin)deltaP[i][j]=0;
    for(int fn=2;fn<3;fn++){
        for(int j=0;j<1;j++){
            prep(p,1,2){
                //ファイル読み出し
                TF1* fgaus = new TF1("fgaus","gaus",-1,1);
                string fname = "allbinbase"+to_string(fn)+"_"+to_string(j)+".root";
                string tname = "test_tree"+to_string(p);
                const char* filename = fname.c_str();
                const char* treeName = tname.c_str();
                //PrintEntryInfo(filename,treeName,4300);
                GetDPfit(fn,j,p,whitehist);
                TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
                c1 -> SetMargin(0.14,0.11,0.2,0.1);
                st.Hist(whitehist);
                c1 -> SetLogy();
                whitehist -> Draw();
                //whitehist -> Fit(fgaus);
                /*double sigma = fgaus -> GetParameter("Sigma");
                deltaP[j*2+(2-p)] = sigma;
                filesystem::current_path(whitedir);
                string figname = "Pfit"+to_string(fn)+"_"+to_string(j)+"_"+to_string(p)+".ps";
                c1 -> SaveAs(figname.c_str());*/
            }
        }
    }
    /*rep(i,8){
        cout << i/2 << " : " << i%2 << " : " << deltaP[i] << endl;
    }
    //ここで処理
    /*int pnum = 0;
    rep(bin,nbin){
        int count = 0;
        rep(j,8){
            if(vec[j][bin]==DINF)continue;

        }
    }*/
    
    //GetDPfit();
    /*for(int fn=0;fn<1;fn++){
        for(int offset=9;offset<10;offset++){
            //PrintEntryInfo(filename,treeName,10);
            //WhiteCheck(offset,5,fn,1,whitehist);
            //PeakFit2(offset,5,fn,1,chihist,vec);
            //ChiCheck(offset,fn,normhist);
        }
    }
    gStyle->SetOptFit(1111);
    c1 -> SetLogy();
    st.Hist(normhist);
    normhist -> Draw();
    normhist -> Fit("gaus");*/
    //c1 -> SetLogy();
    /*st.Hist(chihist);
    chihist -> SetMaximum(2000);
    //gStyle->SetOptFit(1111);
    chihist -> Draw();
    int entnum = chihist -> GetEntries();
    //int xmax = chihist -> GetMaximumBin();
    //double a = xmax*0.1*15.0/17.0;
    chifit -> FixParameter(0,25640);
    chifit -> FixParameter(1,0.038461538461538464);
    chifit -> Draw("same");
    double x = chifit -> Integral(0,10);
    cout << x << endl;
    // 
    // chifit -> SetParameter(1,0.04);
    // chihist -> Fit(chifit);
    //詰められたホワイトノイズやχ^2/ndfが正常かどうかを確認するプログラム*/
    return;
}