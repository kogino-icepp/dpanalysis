#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
using namespace std;
#define rep(i,n) for(int i=0;i<n;i++)
#define prep(i,m,n) for(int i=m;i<n;i++)
#define fore(i,v) for(auto& i:v)

//共通の物理量
const double dcsigma = 0.000487;
const double dhsigma = 0.000457;
const double ddcsigma = 0.000508;
const double ddhsigma = 0.000496;
const double dclim = 6.2;
const double dhlim = 6.2;
const double ddclim = 5;
const double ddhlim = 5;
const int nbin=32767;
const int sb = 2621;
const int fb = 28835;
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
const int dbin = 60;//フィッティングするビンの幅
const double kb=1.38*pow(10,-23);
const double df=88.5*pow(10,3);
const double Tc=76;
const double Th=297;

//速度と周波数を変換する関数

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
bool toge_hantei(vector<bool>ctoge, vector<bool>htoge, int bin, queue<int>& que){
    prep(i,bin,bin+dbin){
        if(ctoge[i] || htoge[i]){
            que.push(i);
            return true;
        }
    }
    return false;
}
vector<double> syoki_para(double x0,double x1,double x2,double y0,double y1,double y2){
    vector<double> ans(3);
    //y=a(x-x0)(x-x2)+mx+nが(x1,y1)を通る
    //m=(y2-y0)/(x2-x0), n=
    double m = (y2-y0)/(x2-x0);
    double n = (y0*x2-y2*x0)/(x2-x0);
    double a = (y1-m*x1-n)/((x1-x0)*(x1-x2));
    double b = m-a*(x0+x2);
    double c = a*x0*x2+n;
    ans[0] = a;
    ans[1] = b;
    ans[2] = c;
    return ans;
}
vector<vector<double>> inv_test(int dim,vector<vector<double>> mat){
    vector<vector<double>> inv(dim,vector<double>(dim));
    double sweep[dim][dim*2];
    rep(i,dim){
        rep(j,dim){
            sweep[i][j] = mat[i][j];
            sweep[i][dim+j] = (i==j) ? 1 : 0;
        }
    }
    double a;
    rep(k,dim){
        a = 1 / sweep[k][k];
        rep(j,dim*2){
            sweep[k][j] *= a;
        }
        rep(i,dim) {
            if (i == k)continue;
            a = -sweep[i][k];
            rep(j,dim*2){
                
                sweep[i][j] += sweep[k][j] * a; 
            }
        }
    }
    rep(i,dim) {
        rep(j,dim) {
            inv[i][j] = sweep[i][dim+j];
        }
    }
    return inv;
}
vector<double> syoki_paraN(int dim,vector<vector<double>> point){
    int N = dim+1;
    vector<double> ret(N);
    vector<double> y(N);
    rep(i,N)y[i] = point[i][1];
    /* 逆行列を求める行列用の２次元配列 */
    double mat[N][N];
    /* 逆行列用の２次元配列 */
    vector<vector<double>> inv(N,vector<double>(N));
    /* 掃き出し法を行う行列 */
    double sweep[N][N * 2];
    rep(i,N){
        double x = point[i][0];
        double v = 1.0;
        rep(j,N){
            mat[i][j] = v;
            v *= x;
        }
    }
    rep(i,N){
        rep(j,N){
            /* sweepの左側に逆行列を求める行列をセット */
            sweep[i][j] = mat[i][j];
            /* sweepの右側に単位行列をセット */
            sweep[i][N+j] = (i==j) ? 1 : 0;
        }
    }
    double a;
    rep(k,N){
        a = 1 / sweep[k][k];
        rep(j,N*2)sweep[k][j] *= a;
        rep(i,N) {
            if (i == k)continue;
            a = -sweep[i][k];
            rep(j,N*2){
                sweep[i][j] += sweep[k][j] * a; 
            }
        }
    }
    rep(i,N) {
        rep(j,N) {
            inv[i][j] = sweep[i][N + j];
        }
    }
    
    /*逆行列をかけて係数の初期パラメータを計算*/
    rep(i,N){
        double value = 0;
        rep(j,N)value += inv[i][j]*y[j];
        ret[i] = value;
    }
    return ret;
}

void Setting_Graph(TGraph* graph,double xmin,double xmax,double ymin,double ymax){
    graph -> SetMaximum(ymax);
    graph -> SetMinimum(ymin);
    graph -> GetXaxis() -> SetLimits(xmin,xmax);
    graph -> GetXaxis() -> SetLabelSize(0.08);
    graph -> GetXaxis() -> SetTitleOffset(0.7);
    graph -> GetXaxis() -> SetTitleSize(0.08);
    graph -> GetYaxis() -> SetTitleOffset(0.5);
    graph -> GetYaxis() -> SetLabelSize(0.08);
    graph -> GetYaxis() -> SetTitleSize(0.08);
    graph -> GetXaxis() -> SetNdivisions(505);
    graph -> GetYaxis() -> SetNdivisions(505);
}

void Setting_GraphErrors(TGraphErrors* graph,double xmin,double xmax,double ymin,double ymax){
    graph -> SetMaximum(ymax);
    graph -> SetMinimum(ymin);
    graph -> GetXaxis() -> SetLimits(xmin,xmax);
    graph -> GetXaxis() -> SetLabelSize(0.08);
    graph -> GetXaxis() -> SetTitleOffset(0.9);
    graph -> GetXaxis() -> SetTitleSize(0.08);
    graph -> GetYaxis() -> SetTitleOffset(0.8);
    graph -> GetYaxis() -> SetLabelSize(0.08);
    graph -> GetYaxis() -> SetTitleSize(0.08);
    graph -> GetXaxis() -> SetNdivisions(505);
    graph -> GetYaxis() -> SetNdivisions(505);
}
void Setting_Hist(TH1D* hist){
    hist -> GetXaxis() -> SetLabelSize(0.08);
    hist -> GetXaxis() -> SetTitleOffset(0.8);
    hist -> GetXaxis() -> SetTitleSize(0.08);
    hist -> GetYaxis() -> SetTitleOffset(0.7);
    hist -> GetYaxis() -> SetLabelSize(0.08);
    hist -> GetYaxis() -> SetTitleSize(0.08);
    hist -> GetXaxis() -> SetNdivisions(505);
    hist -> GetYaxis() -> SetNdivisions(505);
}

//パラメータ数を可変にしたときのフィッティング関数
void Fitting_N(TGraphErrors* graph, TF1* f, vector<double> &sp,int ite, double fxmin,double fxmax,Option_t* t){
    int l = sp.size();
    //初期パラメータを関数に渡す処理
    rep(i,l)f -> SetParameter(i,sp[i]);
    rep(i,ite){
        graph -> Fit(f, t, "", fxmin, fxmax);
    }
    rep(j,l)sp[j] = f -> GetParameter(j);
}

void peak_search2(){
    //お絵描きの設定など
    Double_t xlo=215;
    Double_t xhi=264;
    Double_t ylo=0;
    Double_t yhi=100;
 
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    TH1F *frame=gPad->DrawFrame(xlo,ylo,xhi,yhi);
    filesystem::path path=filesystem::current_path();
    string dir="/Users/oginokyousuke/code/work_bench3/data_all/";
    vector<int> sbin={0,512,-512,-1024};
    filesystem::current_path(dir);

    //諸々のヒストグラム
    TH1D* wnoise_hist = new TH1D("wnoise_hist","white_noise;white_noise[K];Count",100,-2,2);
    TH1D* out_hist = new TH1D("out_hist","out_hist;num_of_out;Count",9,-0.5,8.5);
    TH1D* chi_hist2 = new TH1D("chi_hist2","chi_squre/NDF;chi_squre/NDF;Count",100,0,10);
    TH1D* chi_hist3 = new TH1D("chi_hist3","chi_squre/NDF;chi_squre/NDF;Count",100,0,10);
    TH1D* chi_hist4 = new TH1D("chi_hist4","chi_squre/NDF;chi_squre/NDF;Count",100,0,10);
    vector<vector<int>> toge_or_not(25,vector<int>((fb-sb)/60,0));
    TH1D* fit_amp = new TH1D("fi_amp","Fit_Amp;amp;Count",100,-10,10);
    for(int i=1;i<25;i++){
        double ifmin = 213.7+2*i;
        double ifmax = 216.2+2*i;
        TGraph* graph1 = new TGraph;
        string cdir=dir+"band"+to_string(i);
        filesystem::current_path(cdir);
        for(int j=0;j<4;j++){
            TGraph* graph = new TGraph;
            double Freq1[nbin],cold1[nbin],hot1[nbin],mirror1[nbin],Freq2[nbin],cold2[nbin],hot2[nbin],mirror2[nbin];
            //ビンのシフトをここでいじる
            for(int data=0;data<576;data++){
                string file_name="test"+to_string(data)+".root";
                if(filesystem::exists(file_name)){
                    TFile* file = new TFile(file_name.c_str());
                    TTree* tree = (TTree*)file->Get("tree");
                    int time,xfft,shift,body,num;
                    long double freq,power;
                    tree->SetBranchAddress("time",&time);
                    tree->SetBranchAddress("XFFT",&xfft);
                    tree->SetBranchAddress("shift",&shift);
                    tree->SetBranchAddress("freq",&freq);
                    tree->SetBranchAddress("power",&power);
                    tree->SetBranchAddress("num",&num);
                    tree->SetBranchAddress("body",&body);
                    tree->GetEntry(0);
                    if(shift==sbin.at(j)){
                        //cout << file_name << endl;
                        if(num==1 && body==0){
                            for(int bin=0;bin<nbin;bin++){
                                tree -> GetEntry(bin);
                                Freq1[bin]=freq;
                                cold1[bin]=power;
                            }
                        }
                        else if(num==1 && body==1){
                            for(int bin=0;bin<nbin;bin++){
                                tree->GetEntry(bin);
                                hot1[bin]=power;
                            }
                        }
                        else if(num==1 && body==2){
                            for(int bin=0;bin<nbin;bin++){
                                tree->GetEntry(bin);
                                mirror1[bin]=power;
                            }
                        }
                        else if(num==0 && body==0){
                            for(int bin=0;bin<nbin;bin++){
                                tree->GetEntry(bin);
                                cold2[bin]=power;
                                Freq2[bin]=freq;
                            }
                        }
                        else if(num==0 && body==1){
                            for(int bin=0;bin<nbin;bin++){
                                tree->GetEntry(bin);
                                hot2[bin]=power;
                            }
                        }
                        else if(num==0 && body==2){
                            for(int bin=0;bin<nbin;bin++){
                                tree->GetEntry(bin);
                                mirror2[bin]=power;
                            }
                        }
                    }
                }
            }
            long double cmae11,cmae12,cmae21,cmae22,hmae11,hmae12,hmae21,hmae22;
            long double ddc1,ddc2,ddh1,ddh2;
            long double dc1,dc2,dh1,dh2;
            cmae11 = cold1[sb-1];
            cmae12 = cold1[sb-2];
            cmae21 = cold2[sb-1];
            cmae22 = cold2[sb-2];
            hmae11 = hot1[sb-1];
            hmae12 = hot1[sb-2];
            hmae21 = hot2[sb-1];
            hmae22 = hot2[sb-2];
            
            //cold,hotを用いて棘を抜く手法~173
            //JK曰く以外と捨てデータ広めにとっても被りさえしなければセーフなのでブロックごと(30ビン)捨てるのはあり
            vector<bool> toge_cold1(nbin,false);
            vector<bool> toge_cold2(nbin,false);
            vector<bool> toge_hot1(nbin,false);
            vector<bool> toge_hot2(nbin,false);
            
            for(int bin=sb;bin<fb;bin++){
                dc1 = cold1[bin]-cmae11;
                dc2 = cold2[bin]-cmae21;
                ddc1 = cold1[bin]+cmae12-2*cmae11;
                ddc2 = cold2[bin]+cmae22-2*cmae21;
                ddh1 = hot1[bin]+hmae12-2*hmae11;
                ddh2 = hot2[bin]+hmae22-2*hmae21;

                dc1 /= cold1[bin]*dcsigma;
                dc2 /= cold2[bin]*dcsigma;
                ddc1 /= cold1[bin]*ddcsigma;
                ddc2 /= cold2[bin]*ddcsigma;
                ddh1 /= hot1[bin]*ddhsigma;
                ddh2 /= hot2[bin]*ddhsigma;

                cmae12 = cmae11;
                cmae11 = cold1[bin];
                cmae22 = cmae21;
                cmae21 = cold2[bin];
                hmae12 = hmae11;
                hmae11 = hot1[bin];
                hmae22 = hmae21;
                hmae21 = hot2[bin];
                
                if(abs(ddc1)>ddclim)toge_cold1[bin]=true;
                if(abs(ddc2)>ddclim)toge_cold2[bin]=true;
                if(abs(ddh1)>ddhlim)toge_hot1[bin]=true;
                if(abs(ddh2)>ddhlim)toge_hot2[bin]=true;
                if(toge_cold1[bin] && toge_hot1[bin])toge_or_not[i][bin]++;
                if(toge_cold2[bin] && toge_hot2[bin])toge_or_not[i][bin]++;
            }
            
            vector<int> sa1(nbin,-1);
            vector<int> sa2(nbin,-1);
            TGraphErrors* pgraph1 = new TGraphErrors;
            TGraphErrors* pgraph2 = new TGraphErrors;
            //TGraph* pgraph1 = new TGraph;
            //TGraph* pgraph2 = new TGraph;
            long double gain1,psys1,prec1,gain2,psys2,prec2;
            
            for(int bin=sb;bin<fb;bin++){
                //使えないところを60ビンのブロック丸々抜く
                //if(toge_cold1[bin] || toge_cold2[bin]) continue;
                
                gain1=(hot1[bin]-cold1[bin])/(2*kb*(Th-Tc)*df);
                psys1=(cold1[bin]/gain1)-2*kb*Tc*df;
                prec1=((mirror1[bin]/gain1)-psys1)/(2*kb*df);
                pgraph1 -> SetPoint(bin-sb,Freq1[bin],prec1);
                pgraph1 -> SetPointError(bin-sb,0,0.1);
            
                gain2=(hot2[bin]-cold2[bin])/(2*kb*(Th-Tc)*df);
                psys2=(cold2[bin]/gain2)-2*kb*Tc*df;
                prec2=((mirror2[bin]/gain2)-psys2)/(2*kb*df);
                pgraph2->SetPoint(bin-sb,Freq2[bin],prec2);
                pgraph2 -> SetPointError(bin-sb,0,0.1);
            }
            pgraph1->SetTitle("base_fit_test;Freq[GHz];Prec[K]");
            //graph_cold1 -> SetMinimum(-0.01);
            Setting_GraphErrors(pgraph1,ifmin,ifmax,0,100);
            //pgraph2->GetXaxis()->SetLimits(0,nbin);
            pgraph1 -> SetLineColor(kGreen);
            pgraph1 -> Draw("AP");
            pgraph2 -> SetLineColor(kBlue);
            //pgraph2 -> Draw("L");
            //cout << "out freq = " << Freq1[23140] << endl;
            
            vector<vector<double>> fit_param1(3);
            vector<vector<double>> fit_param2(3);
            double p01,p11,p21,p31,p02,p12,p22,p32;
            
            TH1D* white_hist = new TH1D("white_hist","white_hist;white_noise[K];Count",100,-100,100);
            
            Double_t Fmin;
            Double_t fedm;
            Double_t errdef;
            Int_t npari;
            Int_t nparx;
            Int_t istat;
            queue<int> que;
            TH1D* fit_result = new TH1D("fit_result","mnstatus;istat;Count",4,-0.5,3.5);
            double chi2_1,chi2_2;
            int Ndof1, Ndof2;
            vector<bool> heta_fit(nbin,false);
            vector<double> heta_fit2(nbin);
            int dim = 2;
            int N = dim + 1;
            //試しに三次や高次の関数でのフィットが上手くいくか検証
            TGraph* white_noise1 = new TGraph;
            TGraph* white_noise2 = new TGraph;
            double wnoise1,wnoise2;
            int bin_true1 = 0;
            int bin_true2 = 0;

            for(int bin=sb;bin<fb;bin+=dbin){
                if(toge_hantei(toge_cold1,toge_hot1,bin,que))continue;
                double fm1 = min(Freq1[bin],Freq1[bin+dbin]);
                double fM1 = max(Freq1[bin],Freq1[bin+dbin]);
                double fm2 = min(Freq2[bin],Freq2[bin+dbin]);
                double fM2 = max(Freq2[bin],Freq2[bin+dbin]);

                //フィットする関数一覧
                TF1* f2_1 = new TF1("base1","[0]+[1]*x+[2]*x*x",fm1,fM1);
                TF1* f2_2 = new TF1("base2","[0]+[1]*x+[2]*x*x",fm2,fM2);
                TF1* f3_1 = new TF1("base1","[0]+[1]*x+[2]*x*x+[3]*x*x*x",fm1,fM1);
                TF1* f3_2 = new TF1("base2","[0]+[1]*x+[2]*x*x+[3]*x*x*x",fm2,fM2);
                TF1* f4_1 = new TF1("base1","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",fm1,fM1);
                TF1* f4_2 = new TF1("base2","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",fm2,fM2);
                vector<vector<double>> sp1(N,vector<double>(2));
                vector<vector<double>> sp2(N,vector<double>(2));

                rep(k,N){
                    sp1[k][0] = pgraph1 -> GetPointX(bin+(dbin/dim)*k-sb);
                    sp1[k][1] = pgraph1 -> GetPointY(bin+(dbin/dim)*k-sb);
                    sp2[k][0] = pgraph2 -> GetPointX(bin+(dbin/dim)*k-sb);
                    sp2[k][1] = pgraph2 -> GetPointY(bin+(dbin/dim)*k-sb);
                }
                
                vector<double> sp_fit_1 = syoki_paraN(dim,sp1);
                vector<double> sp_fit_2 = syoki_paraN(dim,sp2);

                Fitting_N(pgraph1,f2_1,sp_fit_1,5,fm1,fM1,"+");
                Fitting_N(pgraph2,f2_2,sp_fit_2,5,fm1,fM1,"+");

                prep(cbin,bin,bin+dbin){
                    if(cbin>=fb)break;
                    wnoise1 = pgraph1 -> GetPointY(cbin-sb);
                    wnoise1 -= sp_fit_1[0]+sp_fit_1[1]*Freq1[cbin]+sp_fit_1[2]*Freq1[cbin]*Freq1[cbin];
                    wnoise_hist -> Fill(wnoise1);
                    white_noise1 -> SetPoint(bin_true1,Freq1[cbin],wnoise1);
                    bin_true1++;
                }
                prep(cbin,bin,bin+dbin){
                    if(cbin>=fb)break;
                    wnoise2 = pgraph2 -> GetPointY(cbin-sb);
                    wnoise2 -= sp_fit_2[0]+sp_fit_2[1]*Freq2[cbin]+sp_fit_2[2]*Freq2[cbin]*Freq2[cbin];
                    wnoise_hist -> Fill(wnoise2);
                    white_noise2 -> SetPoint(bin_true2,Freq2[cbin],wnoise2);
                    bin_true2++;
                }
            
                
                /*chi2_1 = f2_1 -> GetChisquare();
                Ndof1 = f2_1 -> GetNDF();
                chi_hist2 -> Fill(chi2_1/Ndof1);
                chi2_2 = f2_2 -> GetChisquare();
                Ndof2 = f2_2 -> GetNDF();
                chi_hist2 -> Fill(chi2_2/Ndof2);

                sp_fit_1.push_back(0);
                sp_fit_2.push_back(0);
                Fitting_N(pgraph1,f3_1,sp_fit_1,5,fm1,fM1,"+");
                Fitting_N(pgraph2,f3_2,sp_fit_2,5,fm1,fM1,"+");

                chi2_1 = f3_1 -> GetChisquare();
                Ndof1 = f3_1 -> GetNDF();
                chi_hist3 -> Fill(chi2_1/Ndof1);
                chi2_2 = f3_2 -> GetChisquare();
                Ndof2 = f3_2 -> GetNDF();
                chi_hist3 -> Fill(chi2_2/Ndof2);

                sp_fit_1.push_back(0);
                sp_fit_2.push_back(0);
                Fitting_N(pgraph1,f4_1,sp_fit_1,5,fm1,fM1,"+");
                Fitting_N(pgraph2,f4_2,sp_fit_2,5,fm1,fM1,"+");

                chi2_1 = f4_1 -> GetChisquare();
                Ndof1 = f4_1 -> GetNDF();
                chi_hist4 -> Fill(chi2_1/Ndof1);
                chi2_2 = f4_2 -> GetChisquare();
                Ndof2 = f4_2 -> GetNDF();
                chi_hist4 -> Fill(chi2_2/Ndof2);*/
    
            }
            /*
            Setting_Graph(white_noise1,215.8,218.2,-5,5);
            white_noise1 -> SetMarkerStyle(20);
            white_noise2 -> SetMarkerStyle(20);
            white_noise1 -> SetMarkerSize(0.5);
            white_noise2 -> SetMarkerSize(0.5);
            white_noise1 -> SetMarkerColor(kGreen);
            white_noise2 -> SetMarkerColor(kGreen);
            white_noise1 -> SetTitle("white_noise;Freq[GHz];Prec[K]");
            white_noise1 -> Draw("AP");
            //else white_noise1 -> Draw("L");
            white_noise2 -> Draw("P");
            
            ホワイトノイズの中にピークがないか探す作業
            TF1* fit_func = new TF1("fit_func","F_sig2(x,[0],[1],0.5)+[2]");
            double base = 0.01;
            for(int bin=sb;bin<fb-dbin;bin++){
                double p = white_noise1 -> GetPointY(bin-sb);
                fit_func -> SetParameter(0,Freq1[bin]);
                fit_func -> SetParameter(1,p);
                fit_func -> SetParameter(2,base);
                white_noise1 -> Fit(fit_func,"","C",Freq1[bin],Freq1[bin+dbin]);
                p = fit_func -> GetParameter(1);
                fit_amp -> Fill(p);
            }*/
            //このままフィットしてしまってもいいが毎回回すのはしんどいので一旦保存して改めて別の場所でフィット回すでも良い
        }
    }
    //c1 -> SetLogy();
    Setting_Hist(wnoise_hist);
    wnoise_hist -> Draw();
    wnoise_hist -> Fit("gaus");
    /*int free_num = 57;vi 
    TF1* chi_test = new TF1("chi_test","1000*pow(x*57,([0]/2)-1)*exp(-x*57/2)/(pow(2,[0]/2)*tgamma([0]/2))");
    chi_teerst -> SetParameter(0,57);soc
    //c1 -> SetLogy();
    Setting_Hist(chi_hist2);
    chi_hist2 -> SetLineColor(kBlue);
    chi_hist3 -> SetLineColor(kRed);
    chi_hist4 -> SetLineColor(kGreen);
    chi_hist2 -> Draw();
    chi_hist3 -> Draw("same");
    chi_hist4 -> Draw("same");*/
    //cout << chi_test -> GetMaximum() << endl;
    //chi_test -> Draw("same");
}