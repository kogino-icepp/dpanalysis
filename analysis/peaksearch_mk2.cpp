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
/*double F_sig_delta(double f,double f0,double P,double r,double delF){
    if(f-r*dNu-delF<=0){
        return ;
    }
    else if(f-r*dNu-delF>0){
        return 0;
    }

}*/
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

//offsetとかいうゴミ概念を除いた上でベースラインフィットがどのくらいうまくいっていたのかの確認
void ChiCheck2(int i,int j,int p,TH1D*hist,TH1D*hist2){
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
void ChiSep(int i,int j,int p){
    Setting st;
    st.markerstyle = 20;
    st.dot_size = 0.5;
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
    double vchi[nbin];
    rep(i,nbin)vchi[i] = DINF;
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
        vchi[bin] = chi;
    }
    cout << vpfreq[15289] << endl;
    TGraph* graph = new TGraph;
    double chimin = 1000000.0;
    int minnum = 0;
    int num = 0;
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    //binをずらしながらヒストグラムに詰めて各フィット→合計値が一番いい結果と実際の場所を眺めてみる(というかプロット作ってみる)
    // prep(bin,14000,18000){
    //     string f1name = "f1"+to_string(bin);
    //     string f2name = "f2"+to_string(bin);

    //     TF1* f1 = new TF1(f1name.c_str(),"chiF_freefit(x,[0],[1],27,0.05)",0,5);
    //     TF1* f2 = new TF1(f2name.c_str(),"chiF_freefit(x,[0],[1],27,0.05)",0,5);
    //     int num1 = 0;
    //     int num2 = 0;
    //     TH1D* hist1 = new TH1D("hist1","",100,0,5);
    //     TH1D* hist2 = new TH1D("hist2","",100,0,5);
    //     prep(bin1,sb,bin){
    //         if(vchi[bin1]==DINF)continue;
    //         hist1 -> Fill(vchi[bin1]);
    //         num1++;
    //     }
    //     prep(bin2,bin,fb){
    //         if(vchi[bin2]==DINF)continue;
    //         hist2 -> Fill(vchi[bin2]);
    //         num2++;
    //     }
    //     f1 -> FixParameter(0,num1);
    //     f2 -> FixParameter(0,num2);
    //     f1 -> SetParameter(1,0.02);
    //     f2 -> SetParameter(1,0.02);
    //     rep(ite,5)hist1 -> Fit(f1,"Q","",0,5);
    //     rep(ite,5)hist1 -> Fit(f1,"MQ","",0,5); 
    //     rep(ite,5)hist1 -> Fit(f1,"EQ","",0,5); 
    //     rep(ite,5)hist2 -> Fit(f2,"Q","",0,5);
    //     rep(ite,5)hist2 -> Fit(f2,"MQ","",0,5); 
    //     rep(ite,5)hist2 -> Fit(f2,"EQ","",0,5); 

    //     double chi1 = f1 -> GetChisquare();
    //     double chi2 = f2 -> GetChisquare();
    //     cout << chi1+chi2 <<  endl;
    //     if(chi1+chi2<chimin){
    //         chimin = chi1+chi2;
    //         minnum=bin;
    //     }
    //     graph -> SetPoint(num,bin,chi1+chi2);
    //     num++;
    //     st.Hist(hist1);
    //     st.Hist(hist2);
    //     hist1 -> SetLineColor(kBlue);
    //     hist2 -> SetLineColor(kRed);
    //     hist1 -> Draw();
    //     hist2 -> Draw("same");
    // }
    // cout << "minnum : " << minnum << endl;
    /*axrange ax = {0,nbin,0,40000,0,1,"chi1+chi2;Bin;Chi1+Chi2"};
    st.Graph(graph,ax);
    graph -> SetMarkerStyle(20);
    graph -> SetMarkerColor(kBlack);
    graph -> Draw("AP");*/
}
//limitを計算できるように配列を入れるようにする、入らなかったら泣く
void GetDPfit(int i,int j,int p,double (&dlist)[nbin],double (&deltaP)[nbin],TH1D* hist,TH1D* hist2){
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
    TGraph* pgraph = new TGraph;
    GetBasicData(i,j,p,pgraph);

    //vectorにとりあえずのフィット結果を詰める
    //どこかにグローバルなエラーと局所的なエラーを両方出すプログラムが欲しい
    int xfft = XFFT(i);
    //vector<double> Pfit(nbin,DINF);
    int outlist[20];
    double outvalue[20];
    int index = 0;
    int arcount[nbin];
    int cstart = 0;
    double poutlist[nbin];
    TH1D* whitehist = new TH1D("whitehist","P_{fit};P_fit[K*Hz];Count",100,-1,1);
    for(int bin=sb;bin<fb;bin++){
        if(vparas[0][bin]==DINF &&vparas[1][bin]==DINF  &&vparas[2][bin]==DINF){
            //cout << "Pass" << endl;
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
        //
        double sigma = 0;
        
        rep(k,dbin){
            double x = pgraph -> GetPointX(bin-sb+k);
            double y = pgraph -> GetPointY(bin-sb+k);
            fitgraph -> SetPoint(k,x,y);
            if(k<10 || k>=20){
                sigma += pow(y-(a*(x-b)*(x-b)+c),2);
            }
            if(k==bpos){
                if(y>a*(x-b)*(x-b)+c)peakquad -> SetParameter(4,0.1);
                else peakquad -> SetParameter(4,-0.1);
            }
        }
        sigma /= 17;
        sigma = sqrt(sigma);
        rep(k,dbin)fitgraph -> SetPointError(k,0,sigma);//エラーを手で変更
        //cout << sigma << endl;
        //if(sigma<0.5)hist -> Fill(sigma);
        //else hist -> Fill(0.48);
        st.GraphErrors(fitgraph,axtest);
        
        //if(sigma>0.2)cout << bin << " " << sigma << endl;
        fitgraph -> Draw("AP");
        peakquad -> FixParameter(3,mfreq);
        peakquad -> SetParameter(4,0.1);//ベースラインより下なら負、上なら正の初期値にしてみるとか？
        rep(ite,5)fitgraph -> Fit(peakquad,"Q","",sfreq,ffreq);
        rep(ite,5)fitgraph -> Fit(peakquad,"EQ","",sfreq,ffreq);
        //フィットがある程度収束するまでこれ続ける
        double fitres = 10000000;
        double fitres2 = 10000001;
        double fiterr = 500000;//不本意だが
        int count = 0;
        //案1: エラーがある程度小さくなるまで忖度し続ける-> そもそもlimit計算の段階で謎のエラーが出力されているのはなぜ?
        rep(ite,100000){
            fitres2 = fitres;
            fitgraph -> Fit(peakquad,"MQ","",sfreq,ffreq);
            /*if(result.Get() != nullptr && result->IsValid()){
                cout << "sucess!" << endl;
            }
            else cout << "Fail" << endl;*/
            fitres = peakquad -> GetParameter(4);
            fiterr = peakquad -> GetParError(4);
            if(abs(fitres-fitres2)<0.0000001){
                if(fiterr<1){
                    //cout << "ite: " << ite << endl;
                    break;}
                else continue;
            }
            //cout << "ite " << ite << endl;
        }
        /*while(abs(fitres-fitres2)>0.0000001){
            if(count==10000)break;
            fitres2 = fitres;
            fitgraph -> Fit(peakquad,"MQ","",sfreq,ffreq);
            fitres = peakquad -> GetParameter(4);
            fiterr = peakquad -> GetParError(4);
            count++;
        }*/
        arcount[cstart] = count;
        cstart++;
        double pout = peakquad -> GetParameter(4);
        double dpout = peakquad -> GetParError(4);
        //pout *= 2*kb*df;
        fitgraph -> Draw("AP");
        peakquad -> Draw("same");
        dlist[bin] = pout;
        deltaP[bin] = dpout;
        double chi = peakquad -> GetChisquare();
        int ndf = peakquad -> GetNDF();
        hist -> Fill(chi/ndf);
        whitehist -> Fill(pout);
        poutlist[index] = pout;
        index++;
        //cout << pout << " " << dpout << endl;
        //hist -> Fill(pout/0.242115);
        //if(vpchi[bin]<3)hist2 -> Fill(pout/0.242115);
        /*if(abs(pout/0.242115)>5){
            outlist[index] = bin;
            outvalue[index] = pout/0.242115;
            index ++;
        }
        //deltaP[i*2+j+(2-p)][bin+11] = pout;*/
    }
    TF1* fgaus = new TF1("fgaus","gaus",-1,1);
    whitehist -> Fit(fgaus);
    double DeltaP = fgaus -> GetParameter("Sigma");
    rep(bin,index)hist2 -> Fill(poutlist[bin]/DeltaP);
    
    //rep(bin,cstart)cout << "count: "<< arcount[bin] << endl;
    
    //cout << "DeltaP : " << DeltaP << endl;
    //将来的にはピークサーチした後の点を統計取ってlimitつけるようにする、きっとできる(vector<double>かdouble ??[]に格納しておけば良いはず)
}
/// @brief 
/// @param dlist //fitで得られたPfit[K]を格納、
/// @param deltaP //各測定に対し、結果次第で適切につけるもの→これ解釈間違ってた、フィットごとに個別に得られる結果なので8×nbinに修正する
/// @param i 
void MakeLimit(double (&dlist)[8][nbin],double (&deltaP)[8][nbin],int i){
    //獲得したpoutのデータ一覧から平均した値を出力+配列の確保を行う
    int xfft = XFFT(i);
    int shift[4];
    if(xfft%2==1)shift[0] = 0,shift[1]=512,shift[2]=-512,shift[3]=-1024;
    else shift[0] = 0,shift[1]=-512,shift[2]=512,shift[3]=1024;
    TGraph* glimit = new TGraph;
    int gbin = 0;
    //deltaPの平均化
    //周波数をどう引っ張ってくるかが問題、割合どうせ分かってるもんなあ
    //TH1* delhist = nullptr;
    prep(bin,sb,fb){
        //探索されていない場合をパスしてそれ以外を足し上げる(何回足しあげたかはきちんとカウント)
        int num = 0;
        double dlim = 0;
        double vardeltaP = 0;
        rep(ite,8){
            int j = ite/2;
            if(dlist[ite][bin+shift[j]] != DINF){
                num++;
                dlim += dlist[ite][bin+shift[j]];
                vardeltaP += deltaP[ite][bin+shift[j]];
            }
            
        }
        dlim /= num;
        vardeltaP /= num;
        
        if(num>0){
            //横軸周波数◯ →　縦軸を[K]からχに変換したい
            double freq;
            //周波数の導出
            if(i%2==1)freq = (213.8+i*2)+0.0000762939*bin;
            else freq = (216.2+i*2)-0.0000762939*bin;
            //limの計算、正負で対応が異なる
            double chilim;
            if(dlim>0)chilim = PtoChi((dlim+1.96*vardeltaP)*2*kb*df);
            else chilim = PtoChi(1.96*vardeltaP*2*kb*df);
            if(chilim>6*pow(10,-11)){
                cout << chilim << " : " << bin << " : " << num <<  endl;
                cout << "dlim: " << dlim << " && vardeltaP: " << vardeltaP << endl;
                //delhist -> Fill(vardeltaP);
            }
            glimit -> SetPoint(gbin,freq,chilim);
            gbin++;
        }
    }
    //cout << "fit num: " << 1750 << " ~ "  << 2000<< endl;
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    double fmin = 213.5+i*2;
    double fmax = 216.5+i*2;
    axrange axlim = {fmin,fmax,pow(10,-11),pow(10,-7),0,1,"#Delta P_{95%}"+to_string(i)+";Freq[GHz];#Delta P_{95%}"};
    st.Graph(glimit,axlim);
    glimit -> SetLineColor(kBlue);
    c1 -> SetLogy();
    glimit -> Draw("AL");
    filesystem::current_path(saveexe);
    string fname = "limit" + to_string(i) +".ps";
    c1 -> SaveAs(fname.c_str());
    //横軸の周波数だけなんとか引っ張れるかどうか
}
//これがメイン関数
void peaksearch_mk2(){
    //まずはそれぞれ別々の情報が詰められているか確認する
    double maxbin = 0.5;
    TF1* chifit = new TF1("chifit","0.5*chiF_freefit(x,[0],[1],26,0.05)",0,5);
    TF1* chifit2 = new TF1("chifit2","(chiF_freefit(x,[0],[1],27,0.05)+chiF_freefit(x,25992-[0],[2],27,0.05))*0.9",0,5);
    //PrintEntryInfo(filename,treeName,10);
    // TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    // c1 -> SetMargin(0.14,0.11,0.2,0.1);
    string whitedir = "/Users/oginokyousuke/data/white_noise";
    //vector<vector<double>> vec(8,vector<double>(nbin,DINF));
    vector<int> excess;
    double deltaP[8][nbin];//フィットで得られたエラーを格納するための二次配列
    double testlist[8][nbin];
    
    rep(i,8)rep(j,nbin)testlist[i][j] = DINF;
    for(int fn=1;fn<14;fn++){
        for(int j=0;j<4;j++){
            prep(p,1,3){
                //ファイル読み出し
                //描画するヒストグラムの名前を毎回作る
                string nhist = "P_{fit}/#DeltaP(" + to_string(fn)+"_"+to_string(j)+"_"+to_string(p)+");P_{fit}/#DeltaP;Count";
                TH1D* normhist = new TH1D("normhist",nhist.c_str(),100,-10,10);
                TH1D* normhist2 = new TH1D("normhist2",nhist.c_str(),100,-10,10);
                TH1D* sighist = new TH1D("sighist","#Delta P_{fit};#Delta P_{fit}[K*Hz];Count",100,0,0.5);
                TH1D* whitehist = new TH1D("whitehist","P_{fit};P_fit[K*Hz];Count",100,-1,1);
                TH1D* chihist = new TH1D("chihist","chihist;#chi^{2}/NDF;Count",100,0,5);
                TH1D* chihist2 = new TH1D("chihist2","chihist;#chi^{2}/NDF;Count",100,0,5);
                TF1* fgaus = new TF1("fgaus","gaus",-1,1);
                string fname = "allbinbase"+to_string(fn)+"_"+to_string(j)+".root";
                string tname = "test_tree"+to_string(p);
                const char* filename = fname.c_str();
                const char* treeName = tname.c_str();
                //PrintEntryInfo(filename,treeName,10);
                TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
                c1 -> SetMargin(0.14,0.11,0.2,0.1);
                //ChiSep(fn,j,p);
                /*ChiCheck2(fn,j,p,chihist,chihist2);
                st.Hist(chihist);
                c1 -> SetLogy();
                
                chifit -> SetParameter(1,0.01);
                chihist -> Draw();
                chihist -> Fit(chifit);
                chifit2 -> SetParameter(1,0.01);
                chifit2 -> SetParameter(2,0.05);
                chihist -> Draw();
                rep(ite,10)chihist -> Fit(chifit2,"E");
                chihist2 -> SetLineColor(kRed);
                chihist2 -> Draw("same");
                //TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
                //c1 -> SetMargin(0.14,0.11,0.2,0.1);
                //c1 -> SetLogy();*/
                GetDPfit(fn,j,p,testlist[j*2+p-1],deltaP[j*2+p-1],chihist,normhist);
                st.Hist(normhist);
                c1 -> SetLogy();
                sighist -> Draw();
                normhist -> Fit("gaus");
                filesystem::current_path(whitedir);
                string fgausname = "testgaus"+to_string(fn)+"_"+to_string(j)+"_"+to_string(p)+".ps";
                string fchiname = "testchi"+to_string(fn)+"_"+to_string(j)+"_"+to_string(p)+".ps";
                c1 -> SaveAs(fgausname.c_str());
                st.Hist(chihist);
                chihist -> Draw();
                int allnum = chihist -> GetEntries();
                cout << "allnum: " << allnum << endl;
                chifit -> FixParameter(0,allnum);
                chifit -> FixParameter(1,1.0/26);
                chifit -> Draw("same");
                c1 -> SaveAs(fchiname.c_str());
                //normhist -> Fit("gaus");
                
                //ここでファイル保存*/
                //filesystem::current_path(whitedir);
                //string figname = "normhist"+to_string(fn)+"_"+to_string(j)+"_"+to_string(p)+".ps";
                //c1 -> SaveAs(figname.c_str());
                
            }
        }
        MakeLimit(testlist,deltaP,fn);
    }
    //今のうちに得られたリストの処理を考える(後段に作らず別で関数作る？あり！)
    /*要件定義
    1. 2つずつは全く同じ帯域を見ているので安直に平均化
    2. 残りの4つは違いにビンがずれているのでずれを補正しながら平均化
    3. 得られた値に対し、負なら一様に等しい値を、正なら95%水準の値をつける
    */
    return;
}