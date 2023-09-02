#include <TFile.h>
#include <TTree.h>
#include <iostream>

#include "../headers/fitter.h"
#define rep(i,n) for(int i=0;i<n;i++)

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
const int INF = 1e9;
using namespace std;
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

axrange axscale = {0,1,0,1,0,1,"scale_test;xscale;yscale"};
Fitter ft;
vector<int> sbin = {0,512,-512,-1024};
string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
string savedirp = "/Users/oginokyousuke/data/peak_data/";
//string savedirp = "/Users/oginokyousuke/data/peak_data/";


void PrintEntryInfo(const char* filename, const char* treeName, Int_t numEntriesToShow) {

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

void check_chi(const char* filename,const char* treeName){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
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
    int binnum = 100;
    int binmin = 0;
    int binmax = 10;
    double binhaba = (binmax-binmin)/binnum;
    TH1D* chi_hist = new TH1D("chi_hist","chi_hist;Chi2/NDF;Count",100,0,10);
    vector<double> vchi;
    vector<vector<double>> vpara(3);
    for (Int_t entry = 0; entry < numEntries; ++entry) {
        tree->GetEntry(entry);
        // Print entry information here
        // For example, print the values of specific branches
        Double_t a, b, c, chi;
        int bin;
        tree->SetBranchAddress("a", &a);
        tree->SetBranchAddress("b", &b);
        tree->SetBranchAddress("c", &c);
        tree->SetBranchAddress("chi", &chi);
        tree->SetBranchAddress("bin", &bin);
        chi_hist -> Fill(chi);
    }
    c1 -> SetLogy();
    st.Hist(chi_hist);
    chi_hist -> Draw();
    
    TF1* chifit = new TF1("chifit","chiF_freefit(x,[0],[1],27,0.1)");
    TF1* chifit2 = new TF1("chifit2","chiF_freefit(x,[0],[1],27,0.1)+chiF_freefit(x,[2],[3],27,0.1)");
    int maxbin = chi_hist -> GetMaximumBin();
    chifit -> FixParameter(0,numEntries);
    chifit -> SetParameter(1,maxbin*0.1/(25.0));
    cout << maxbin << endl;
    chi_hist -> Fit(chifit,"EQ","",0,10);
    chifit -> Draw("same");
    double p1 = chifit -> GetParameter(1);
    cout << fixed;
    cout << setprecision(30) << "alpha : " << 1-chiF_free(1.5,1,27,p1) << endl;
    double chichi = chifit -> GetChisquare();
    int chindf = chifit -> GetNDF();
    cout << "chi/ndf : " << chichi/chindf << endl;
}
//とりあえずベーシックなデータを取ってこれるようにする

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
            cout << file_name << endl;
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
//フィット結果と元データを元にホワイトノイズを作ってエラー出し
/*
機能要件(随時追加)
1. 保存しておいたフィット結果を用いてホワイトノイズを導出
2. ホワイトノイズ上をピーク関数でフィットしてシグナルを詰める
3. シグナルを正規化して著しくexceedしているものを探す
4. シグナルやbin、sigmaの大きさ、周波数などをttreeに詰めて保存(ここまでやりたい)
*/

void MakeWhite(int i,int j,int p/*,double (&vec)[nbin]*/){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    Setting st;
    st.dot_size=0.8;
    st.markerstyle=20;
    st.color = kGreen;
    string fname = "fit_res"+to_string(i)+"_"+to_string(j)+"_"+to_string(2-p)+".root";
    const char* filename = fname.c_str();
    const char* treeName = "tree1";
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
    //30ビンごとに生データとフィット結果を結合
    /*要件定義
    1. 関数に保存したパラメータ渡す
    2. 念の為一応フィットしてから差分を取る
    3. ヒストグラムに詰めて経過観察
    */
    int testbin = 3035;
    TGraph* pgraph = new TGraph;
    double ifmin = 213.9+2*i;
    double ifmax = 216.1+2*i;
    GetBasicData(i,j,p,pgraph);
    axrange axraw = {ifmin,ifmax,0,100,0,1};
    st.Graph(pgraph,axraw);
    //pgraph->Draw("AL");
    //棘があってパスしている場合もあるのでそれを排除するアルゴ欲しい
    //ホワイトノイズまでは取れるようになったのでここからエラーつけてピークフィットするように改良する
    

    TGraphErrors* wgraph = new TGraphErrors;
    for(int bin=sb;bin<fb;bin+=dbin){
        TH1D* whist = new TH1D("whist","whist;white_noise[K];Count",100,-1,1);
        TF1* quadf = new TF1("quadf","[0]*(x-[1])*(x-[1])+[2]",0,1);
        TF1* gausfit = new TF1("gausfit","gaus",-1,1);
        if(vpfreq[bin]==DINF)continue;
        cout << vparas[0][bin] << " " << vparas[1][bin] << " " << vparas[2][bin] << endl;
        quadf -> SetParameter(0,vparas[0][bin]);
        quadf -> SetParameter(1,vparas[1][bin]);
        quadf -> SetParameter(2,vparas[2][bin]);
        TGraphErrors* spgraph = new TGraphErrors;
        double yscale;
        ft.make_scale(spgraph,pgraph,bin-sb,yscale);
        // /spgraph -> Fit(quadf,"M","",0,1);
        
        //この辺に差分取るアルゴ欲しい(スケール化した後なのか一旦戻すのか→戻そう)
        rep(k,dbin){
            double xValue = spgraph -> GetPointX(k);
            double yValue = quadf -> Eval(xValue);
            double yTrue = spgraph -> GetPointY(k);
            whist -> Fill((yValue-yTrue)*yscale);
            double freqbin = pgraph -> GetPointX(bin-sb+k);
            wgraph -> SetPoint(bin+k,freqbin,(yValue-yTrue)*yscale);
        }
        whist -> Fit(gausfit,"MQ","",-1,1);
        double sigma = gausfit -> GetParameter("Sigma");
        prep(rb,bin,bin+dbin)wgraph -> SetPointError(rb,0,sigma);
        st.GraphErrors(spgraph,axscale);
        /*st.Hist(whist);
        whist -> Fit("gaus");
        whist -> Draw();*/
    }
    //wgraphに対してピークサーチ
    TH1D* phist = new TH1D("phist","phist;P[K*kHz];Count",100,-1,1);
    TH1D* chihist = new TH1D("chihist","chihist1;Chi2/NDF;Count",100,0,10);
    vector<pair<int,double>> fitdata;
    
    for(int bin=22685;bin<22686;bin++){
        double sfreq = pgraph -> GetPointX(bin-sb-10);
        double ffreq = pgraph -> GetPointX(bin-sb+10);
        double mfreq = pgraph -> GetPointX(bin-sb);
        //if(sfreq>ffreq)swap(sfreq,ffreq);
        int binnum = 0;
        for(int wbin=sb;wbin<fb;wbin++){
            double xvalue = wgraph -> GetPointX(wbin);
            if(xvalue>=ffreq && xvalue<=sfreq)binnum++;
        }
        //cout << "binnum = " << binnum << endl;
        TF1* peakf = new TF1("peakf","F_sig2(x,[0],[1],0.5)",sfreq,ffreq);
        peakf -> FixParameter(0,mfreq);
        peakf -> SetParameter(1,0.1);
        rep(k,1)wgraph -> Fit(peakf,"EQM","",ffreq,sfreq);
        int ndf = peakf -> GetNDF();
        double chi2 = peakf -> GetChisquare();
        if(binnum<10){
            cout << mfreq << "GHz Can't Fit" << endl;
            
            continue;
        }
        else{
            double p = peakf -> GetParameter(1);
            double perr = peakf -> GetParError(1);
            //cout << "freq : " << mfreq << endl;
            //cout << "p: " << p << " <=> perr: " << perr << endl;
            fitdata.push_back({bin,p});
            phist -> Fill(p);
            chihist -> Fill(chi2/ndf);
        }
        axrange axexe = {ffreq,sfreq,-1,1,0,1,"fitexample;Freq[GHz];WhiteNoise[K]"};
        st.GraphErrors(wgraph,axexe);
        wgraph -> Draw("AP");
        peakf -> Draw("same");
    }
    /*TF1* peakgaus = new TF1("peakgaus","gaus",-1,1);
    //c1 -> SetLogy();
    phist -> Draw();
    phist -> Fit(peakgaus);
    //シグマの値で規格化してどれくらいのシグマなのか一旦詰める
    double psigma = peakgaus -> GetParameter("Sigma");
    int L = fitdata.size();
    TH1D* normhist = new TH1D("normhist","P/dP;P/dP;Counts",100,-5,5);

    rep(i,L){
        double np = fitdata[i].second/psigma;
        normhist -> Fill(np);
        fitdata[i].second = np;
        //vec[v.first] = np;
    }
    TF1* normgaus = new TF1("normgaus","gaus",-5,5);
    normhist -> Fit(normgaus);
    st.Hist(normhist);
    c1 -> SetLogy();
    normhist -> Draw();
    normgaus -> Draw("same");
    filesystem::current_path(savedirp);
    string gname = "normhist"+to_string(i)+"_"+to_string(j)+"_1.ps";
    //c1 -> SaveAs(gname.c_str());
    double normerr = normgaus -> GetParameter("Sigma");
    cout << "normgaus = " << normerr;
    for(auto v:fitdata){
        if(v.second>3)cout << v.first << " " << v.second << endl;
    }*/
    //見積もったσが1でない場合にそれを補正して再度詰め直す→その上でexceedがある点を探す

}

void chiroot_check() {
    Setting st;
    st.dot_size=0.8;
    st.markerstyle=20;
    st.color = kGreen;
    axrange axsig = {0,nbin,-5,5,0,1};
    filesystem::path path=filesystem::current_path();
    string root_dir = "/Users/oginokyousuke/data/root_file/";
    filesystem::current_path(root_dir);
    
    Int_t numEntriesToShow = 10; // Specify the number of entries to show
    //PrintEntryInfo(filename, treeName, numEntriesToShow);
    //これをi,j(pはまだ1セットしかないから固定で)全体で回せば一通りの結果は出てくる(問題はそれだけでは意味がないこと)
    //データの格納はよさそうなのでこのまま4データをまとめる
    /*やるべきこと
    1. binごとの被りチェック
    2. 周波数ごとの非被りチェック
    */
    //double sigdata[nbin];
    //rep(i,4)rep(j,nbin)sigdata[i][j] = DINF;
    /*rep(j,1){
        MakeWhite(2,j,1);
    }*/
    
    /*TGraph* siggraph = new TGraph;
    for(int i=sb;i<fb;i++)siggraph -> SetPoint(i-sb,i,sigdata[0][i]);
    st.Graph(siggraph,axsig);
    siggraph -> Draw("AP");*/
    for(int i=5;i<6;i++){
        rep(j,1){
            cout << "i : " << i << " <=> j: " << j << endl;
            string fn = "fit_res"+to_string(i)+"_"+to_string(j)+"_1.root";
            const char* filename = fn.c_str(); // Specify the path to your ROOT file
            const char* treeName = "tree1"; // Specify the name of the TTree
            check_chi(filename,treeName);
            PrintEntryInfo(filename, treeName, numEntriesToShow);
            
        }
    }
    
    
    /*for(int i=1;i<25;i++){
        for(int j=0;j<4;j++){
            MakeWhite(i,j,1);
        }
    }*/

}