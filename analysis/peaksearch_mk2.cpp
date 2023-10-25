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
        TF1* cubfunc = new TF1("cubfunc","[0]*x*x*x+[1]*x*x+[2]*x+[3]",sfreq,ffreq);
        TGraphErrors* spgraph = new TGraphErrors;
        double yMin,yscale;
    
        ft.make_scale2(pgraph,bin-sb,yMin,yscale);
        rep(ite,5)spgraph -> Fit(cubfunc,"QE","",sfreq,ffreq);
        st.GraphErrors(spgraph,axscale);
        spgraph -> Draw("AP");
        cubfunc -> Draw("same");
        /*ft.rescale_para(vparas[0][bin],vparas[1][bin],vparas[2][bin],sfreq,yMin,ffreq-sfreq,yscale,peakquad);
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
        
        bin-=offset;*/
    }
    file -> Close();
    return;
}
//offsetとかいうゴミ概念を除いた上でベースラインフィットがどのくらいうまくいっていたのかの確認
void ChiCheck2(int i,int j,int p,TH1D*hist){
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
        if(chi>3)cout << freq << " : " << chi  << endl;
    }

}

//offsetとかいうゴミ概念を除いた解析プログラム、同じように動くので次は4つのデータでもってlimitをつけよう
void GetDPfit(int i,int j,int p,TH1D*hist){
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
    int xfft = XFFT(i);
    vector<double> Pfit(nbin,DINF);

    TH1D* whitehist = new TH1D("whitehist","P_{fit};P_fit[K*Hz];Count",100,-1,1);
    for(int bin=sb;bin<fb;bin++){
        if(vparas[0][bin]==DINF &&vparas[1][bin]==DINF  &&vparas[2][bin]==DINF){
            //cout << "Pass" << endl;
            continue;
        }
        rep(ad,29){
            if(mskhantei[xfft][bin+ad])continue;
        }
        double sfreq,ffreq,mfreq;
        if(i%2==1){
            sfreq = pgraph -> GetPointX(bin-sb);
            mfreq = pgraph -> GetPointX(bin-sb+10);
            ffreq = pgraph -> GetPointX(bin-sb+29);
        }
        else{
            ffreq = pgraph -> GetPointX(bin-sb);
            mfreq = pgraph -> GetPointX(bin-sb+19);
            sfreq = pgraph -> GetPointX(bin-sb+29);
        }
        //ここから周波数を既にずらしにいく
        TF1* peakquad = new TF1("peakquad","[0]*(x-[1])*(x-[1])+[2]+F_sig2(x,[3],[4],0.5)",sfreq,ffreq);
        //TF1* cubfunc = new TF1("cubfunc","[0]*x*x*x+[1]*x*x+[2]*x+[3]",sfreq,ffreq);

        TGraphErrors* spgraph = new TGraphErrors;
        double yMin,yscale;
        double xmin;
        axrange axfit = {sfreq,ffreq,40,80,0,1,"fitgraph;Freq[GHz];Prec[K]"};
        
        ft.make_scale2(pgraph,bin-sb,yMin,yscale);
        ft.rescale_para(vparas[0][bin],vparas[1][bin],vparas[2][bin],sfreq,yMin,ffreq-sfreq,yscale,peakquad);
        axrange axtest = {sfreq,ffreq,0,100,0,1,"SignalFit;Freq[GHz];Prec[K]"};
        //一からプロットし直してエラーをグローバルにつける
        TGraphErrors* fitgraph = new TGraphErrors;
        
        double a = peakquad -> GetParameter(0);
        double b = peakquad -> GetParameter(1);
        double c = peakquad -> GetParameter(2);
        // peakquad -> SetParameter(0,0.1);
        // peakquad -> SetParameter(1,225.57);
        // peakquad -> SetParameter(2,47);
        //cout << a << " " << b << " " << c << endl;
        //ここで解析的にエラーを出す
        double sigma = 0;
        rep(k,dbin){
            double x = pgraph -> GetPointX(bin-sb+k);
            double y = pgraph -> GetPointY(bin-sb+k);
            fitgraph -> SetPoint(k,x,y);
            if(k<10 || k>=20){
                sigma += pow(y-(a*(x-b)*(x-b)+c),2);
            }
        }
        sigma /= 17;
        sigma = sqrt(sigma);
        //cout << sigma << endl;
        //if(sigma<0.5)hist -> Fill(sigma);
        //else hist -> Fill(0.48);
        st.GraphErrors(fitgraph,axtest);
        
        //if(sigma>0.2)cout << bin << " " << sigma << endl;
        fitgraph -> Draw("AP");
        peakquad -> FixParameter(3,mfreq);
        peakquad -> SetParameter(4,0.1);
        rep(ite,5)fitgraph -> Fit(peakquad,"Q","",sfreq,ffreq);
        rep(ite,5)fitgraph -> Fit(peakquad,"MQ","",sfreq,ffreq);
        rep(ite,4)fitgraph -> Fit(peakquad,"EQ","",sfreq,ffreq);

        double pout = peakquad -> GetParameter(4);
        //pout *= 2*kb*df;
        fitgraph -> Draw("AP");
        peakquad -> Draw("same");
        
        whitehist -> Fill(pout);
        
        //if(abs(pout/0.169015)>5)que.push(bin);
        //deltaP[i*2+j+(2-p)][bin+11] = pout;*/
    }
    TF1* fgaus = new TF1("fgaus","gaus",-0.2,0.2);
    whitehist -> Fit(fgaus);
    double sigma = fgaus -> GetParameter("Sigma");
    cout << "sigma : " << sigma << endl;
}

//これがメイン関数
void peaksearch_mk2(){
    
    //まずはそれぞれ別々の情報が詰められているか確認する
    double maxbin = 0.5;
    
    TF1* chifit = new TF1("chifit","chiF_freefit(x,[0],[1],[2],[3])",0,5);
    //PrintEntryInfo(filename,treeName,10);
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    string whitedir = "/Users/oginokyousuke/data/white_noise";
    //vector<vector<double>> vec(8,vector<double>(nbin,DINF));
    vector<int> excess;
    double deltaP[8][nbin];
    rep(i,8)rep(j,nbin)deltaP[i][j]=0;
    //
    for(int fn=7;fn<8;fn++){
        for(int j=0;j<1;j++){
            prep(p,1,3){
                queue<double> que; 
                //ファイル読み出し
                //描画するヒストグラムの名前を毎回作る
                string nhist = "P_{fit}/#DeltaP(" + to_string(fn)+"_"+to_string(j)+"_"+to_string(p)+");P_{fit}/#DeltaP;Count";
                TH1D* normhist = new TH1D("normhist",nhist.c_str(),100,-10,10);
                TH1D* whitehist = new TH1D("whitehist","P_{fit};P_fit[K*Hz];Count",100,-1,1);
                TH1D* chihist = new TH1D("chihist","chihist;#chi^{2}/NDF;Count",100,0,5);
                TF1* fgaus = new TF1("fgaus","gaus",-1,1);
                string fname = "allbinbase"+to_string(fn)+"_"+to_string(j)+".root";
                string tname = "test_tree"+to_string(p);
                const char* filename = fname.c_str();
                const char* treeName = tname.c_str();
                //PrintEntryInfo(filename,treeName,10);
                //ChiCheck2(fn,j,p,chihist);
                //TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
                //c1 -> SetMargin(0.14,0.11,0.2,0.1);
                GetDPfit(fn,j,p,normhist);
                c1 -> SetLogy();
                st.Hist(normhist);
                normhist -> Draw();
                normhist -> Fit("gaus");
                //ここでファイル保存
                filesystem::current_path(whitedir);
                string figname = "normhist"+to_string(fn)+"_"+to_string(j)+"_"+to_string(p)+".ps";
                c1 -> SaveAs(figname.c_str());
                
                
            }
        }
    }
    
    //詰められたホワイトノイズやχ^2/ndfが正常かどうかを確認するプログラム*/
    return;
}