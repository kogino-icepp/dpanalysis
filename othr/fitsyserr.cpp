#include <iostream>
#include <random>
#include "../headers/fitter.h"
using namespace std;
#define PI 3.14159265359
const double c= 3.0*pow(10,8);
const double v0=220000.0;
const double vE=v0;
const double vc=v0;
const double dnu=0.000076296;//ビン幅[GHz]
const double dNu=0.0000885;//周波数分解能[GHz]
const double kb=1.38*pow(10,-23);
const double df=88.5*pow(10,3);
const double Tc=76;
const double Th=297;
const double a = 4.307257;//初期のビームウエスト？(2σ分)
const double Rh = 35.6176;//これ何？？そのうち分かるやろ
const double focul_len = 56.87;//受信機内部のレンズの焦点距離
const double zwindow = 57.32+76.819+171.49;
const double lens_cord[4] = {0,57.32,zwindow,zwindow+1500};

double candpow[6] = {0,200,400,600,800,1000};
double candmass[3] = {220,240,260};
Color_t candcolor[3] = {kBlue,kRed,kGreen};
axrange axscale = {0,1,0,1,0,1,";xscale;yscale"};
double v_conv(double f,double f0){
    double rtn=c*sqrt(1-((f0/f)*(f0/f)));
    return rtn;
}
double initial_w0(double a,double Rh,double freq){
    freq *= pow(10,9);
    double w = a*0.644;
    double lambda = c/freq;
    double rtn = w/(1+pow(PI*w*w/(lambda*Rh),2));
    return rtn;
}
double initial_z0(double a,double Rh,double freq){
    freq *= pow(10,9);
    double w = a*0.644;
    double lambda = c/freq;
    double rtn = Rh/(1+pow(lambda*Rh/(PI*w*w),2));
    return rtn;
}
double gaus_beam_width(double w0,double z0,double freq,double z){
    double lambda = c/(freq*pow(10,9));
    double zc = PI*pow(w0,2)/lambda;
    return w0*sqrt(1+(z-z0)*(z-z0)/(zc*zc));
}
vector<double> lens(double din,double w0_in,double freq){//ここでz0とw0outを出す
    freq *= pow(10,9);
    double lambda = c/freq;
    double zc = PI*w0_in*w0_in/lambda;
    double A = din/focul_len-1;
    double dout = focul_len*(1+A/(pow(A,2)+pow(zc,2)/pow(focul_len,2)));
    double w0_out = w0_in/sqrt(pow(A,2)+pow(zc,2)/pow(focul_len,2));
    vector<double> rtn = {dout,w0_out};
    return rtn;
}

void west_freq(){
    TGraph* graph = new TGraph;
    int index = 0;
    for(double freq=216;freq<261;freq+=1){
        double beam_params[3][2];
        beam_params[0][0] = -initial_z0(a,Rh,freq);
        beam_params[0][1] = initial_w0(a,Rh,freq);
        rep(i,2){
            double din = lens_cord[i+1]-(beam_params[i][0]+lens_cord[i]);
            vector<double> ans = lens(din,beam_params[i][1],freq);
            beam_params[i+1][0] = ans[0];
            beam_params[i+1][1] = ans[1];
        }
        graph -> SetPoint(index,freq,beam_params[2][0]);
        index++;
    }
    axrange ax = {216,260,0,300,0,1,";;"};
    //st.Graph(graph,ax);
    graph -> Draw("AC");
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
double F_sig_delta(double f,double f0,double P,double r,double delF){
    double omega0 = f0+delF;
    if((f-f0)-r*dNu-delF<=0 && f>omega0-r*dNu){
        return P*(F_nu(f+r*dNu,omega0)-F_nu(omega0,omega0));
    }
    else if((f-f0)-r*dNu-delF>0){
        return P*(F_nu(f+r*dNu,omega0)-F_nu(f-r*dNu,omega0));
    }
    else return 0;
}
double F_sigscale(double x,double P,double r,double fmin,double delF){
    double bin = x/0.0344827586;
    
    double f = fmin+bin*dnu;//ここをxとfminの関数に変える
    double f0 = fmin+10*dnu;//ここは自明ではあるがx0からf0へ
    return F_sig_delta(f,f0,P,r,delF);
}
//仮説：f0がf0+delFに変わっただけでは？？


void checkHist(){
    TH1D* hist = new TH1D("hist",";Freq[GHz];Power",20,220-dnu,220+19*dnu);
    rep(i,20){
        hist -> SetBinContent(i,F_sig_delta(220+dnu*(i-1),220,1000,0.5,0));
    }
    st.Hist(hist);
    hist -> Draw();
}
void fitsyserr(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    Fitter ft;
    Setting st;
    st.dot_size = 0.8;
    st.markerstyle = 20;
    st.color = kBlue;
    st.lcolor = kBlue;
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_int_distribution<> rand100(0, 99);
    normal_distribution<double> normrand(0,100);
    uniform_real_distribution<double> frand(-0.5,0.5);
    //m_γ=220,240,260GHzにおいて周波数確度が-1~1kHzの不定性を持っている時のフィット精度を検証する
    double f0 = 220;//DM質量に対応する周波数[GHz]
    //TH1* testhist = new TH1D("testhist","test;Freq[GHz];",20,f0-dnu,f0+19*dnu);
    //TH1* testhist2 = new TH1D("testhist2","test;Freq[GHz];",20,f0-dnu,f0+19*dnu);//実験用のヒストグラム、質量をちょっとずらす
    //testhist2 -> SetLineColor(kRed);
    //フィットのエラーを記録した図;
    //仮説
    
    /*TGraph* wgraph1 = new TGraph;
    TGraph* wgraph2 = new TGraph;
    TGraph* wgraph3 = new TGraph;
    rep(f,3){
        rep(j,3){
            rep(bin,1000){
                double z = bin*(lens_cord[j+1]-lens_cord[j])/1000+lens_cord[j];
                double w = gaus_beam_width(beam_params[f][j][0]+lens_cord[j],beam_params[f][j][1],candmass[f],z);
                if(f==0)wgraph1 -> SetPoint(bin,z,w);
                else if(f==1)wgraph2 -> SetPoint(bin,z,w);
                else wgraph3 -> SetPoint(bin,z,w);
            }
        }
    }
    wgraph1 -> Draw();
    //周波数確度に関するエラーの計算、だいぶ結果怪しいけど一旦できたことにして先すすむ。
    rep(f,3){//大きな周波数の刻み(3パターン)
        double Psum0 = 0;
        rep(ite,3000){
            TH1* fithist = new TH1D("fithist","",20,candmass[f]-dnu,candmass[f]+19*dnu);
            TF1* pfunc = new TF1("pfunc","F_sig2(x,[0],0.94*[1],0.5)",candmass[f]-dnu,candmass[f]+19*dnu);
            rep(bin,20){
                double freq = candmass[f]+(bin-1)*dnu;
                double dp = normrand(mt);
                fithist -> SetBinContent(bin,F_sig_delta(freq,candmass[f],1000,0.5,0)+dp);
                fithist -> SetBinError(bin,100);
            }
            pfunc -> FixParameter(0,candmass[f]);
            fithist -> Fit(pfunc,"Q","",candmass[f]-dnu,candmass[f]+19*dnu);
            Psum0 += pfunc -> GetParameter(1);
            delete fithist;
        }
        Pnoshift[f] = Psum0/3000;
        //次に周波数をランダムにずらしてみた時の値
        double Pres[1000];
        double Pmean = 0;
        rep(d,1000){
            double deltaF = frand(mt);
            double Psum = 0;
            rep(ite,3000){
                TH1* fithist = new TH1D("fithist","",20,candmass[f]-dnu,candmass[f]+19*dnu);
                TF1* pfunc = new TF1("pfunc","F_sig2(x,[0],0.94*[1],0.5)",candmass[f]-dnu,candmass[f]+19*dnu);
                rep(bin,20){
                    double freq = candmass[f]+(bin-1)*dnu;
                    double dp = normrand(mt);
                    fithist -> SetBinContent(bin,F_sig_delta(freq,candmass[f],1000,0.5,deltaF*pow(10,-6))+dp);
                    fithist -> SetBinError(bin,100);
                }
                pfunc -> FixParameter(0,candmass[f]);
                fithist -> Fit(pfunc,"Q","",candmass[f]-dnu,candmass[f]+19*dnu);
                Psum += pfunc -> GetParameter(1);
                delete fithist;
            }
            Psum /= 3000;
            Psum /= Pnoshift[f];
            Pres[d] = Psum;
            Pmean += Psum;
        }
        Pmean /= 1000;
        cout << "f: " << candmass[f] << " <=> " << Pmean << endl;
        double Perr = 0;
        rep(d,1000)Perr += (Pmean-Pres[d])*(Pmean-Pres[d]);
        Perr /= 999;
        Perr = sqrt(Perr);
        testgraph -> SetPoint(f,candmass[f],Pmean);
        testgraph -> SetPointError(f,0,Perr);
    }
    axrange axtest = {216,264,0,2,0,1,";Freq[GHz];P_{fit}(#Delta)/P_{fit}(0)"};
    st.GraphErrors(testgraph,axtest);
    testgraph -> Draw("AP");*/

    //周波数分解能由来のエラー、明日中には絶対完成させます
    /*
    要件定義
    1. 220,240,260GHzに対してそれぞれ調査を行う
    2. 周波数のビン幅±35kHzを5kHzずつずらしながら、5000回シードを変えてホワイトノイズを載せる
    3. 真値からのずれをプロットし
    仮説
    これベースラインフィットを先にしないといけないのでは？
    、*/
    TGraphErrors* graph1 = new TGraphErrors;
    TGraphErrors* graph2 = new TGraphErrors;
    TGraphErrors* graph3 = new TGraphErrors;
    double Pres[15];
    double Perr[15];
    axrange axfin = {-35,35,0.7,1.2,0,1,";#Delta_{freq}[MHz];P/#Delta P"};
    
    rep(f,1){
        //もしちょっと上手く行かなかったらここで真値を一旦出すことにする
        //多分やるべきこと：エラーを先につける→あらかじめベースラインフィット→シグナルを重ねてもう一度フィット
        prep(d,-7,8){
            int itenum = 5000;
            double delF = d*5*pow(10,-6);
            double Psum = 0;
            double Pstock[itenum];
            rep(ite,itenum){
                TGraphErrors* spgraph = new TGraphErrors;
                TGraph* pgraph = new TGraph;
                TGraphErrors* spgraphpar = new TGraphErrors;
                rep(bin,30){
                    double freq = 220+(bin-10)*dnu;
                    double dp = normrand(mt);
                    pgraph -> SetPoint(bin,freq,F_sig_delta(freq,220,1000,0.5,0*pow(10,-6))+dp);
                }
                double yscale;
                
                ft.make_scale(spgraph,pgraph,0,yscale);
                rep(bin,30){
                    if(bin<10 || bin>=20){
                        double x = spgraph -> GetPointX(bin);
                        double y = spgraph -> GetPointY(bin);
                        spgraphpar -> SetPoint(bin,x,y);
                        spgraphpar -> SetPointError(bin,0,1);
                    }
                }
                TF1* f1 = new TF1("f1","[0]*(x-[1])*(x-[1])+[2]",0,1);
                TF1* scalepeak = new TF1("scalepeak","[0]*(x-[1])*(x-[1])+[2]+F_sigscale(x,[3],0.5,[4],[5])",0,1);
                double res1;
                st.GraphErrors(spgraph,axscale);
                ft.allfit2(spgraphpar,f1,5,res1);
                double a = f1 -> GetParameter(0);
                double b = f1 -> GetParameter(1);
                double c = f1 -> GetParameter(2);
                double mfreq = spgraph -> GetPointX(10);
                scalepeak -> SetParameter(0,a);
                scalepeak -> SetParameter(1,b);
                scalepeak -> SetParameter(2,c);
                scalepeak -> SetParameter(3,100);
                scalepeak -> FixParameter(4,candmass[f]-10*dnu);
                scalepeak -> FixParameter(5,delF);
                spgraph -> Fit(scalepeak,"EIQ","",0,1);
                //spgraph -> Draw("AP");
                //scalepeak -> Draw("same");
                double spfit = scalepeak -> GetParameter(3);
                Pstock[ite] = spfit*yscale;
                Psum += Pstock[ite];
                //delete fithist;
                //delete pfunc;
            }
            Psum /= itenum;
            double psig = 0;
            rep(ite,itenum)psig += (Pstock[ite]-Psum)*(Pstock[ite]-Psum);
            
            psig = sqrt(psig/(itenum-1));
            Perr[d+7] = psig;
            Pres[d+7] = Psum;
            cout << "Pfit: " << Psum << " +- " << psig << endl;
        }
        if(f==0){
            rep(bin,15){
                graph1 -> SetPoint(bin,(bin-7)*5,Pres[bin]/1000);
                graph1 -> SetPointError(bin,0,Perr[bin]/1000);
                
            }
            st.GraphErrors(graph1,axfin);
            graph1 -> SetLineColor(kBlue);
            graph1 -> SetMarkerColor(kBlue);

        }
        /*else if(f==1){
            rep(bin,15){
                graph2 -> SetPoint(bin,(bin-7)*5,Pres[bin]/1000);
                graph2 -> SetPointError(bin,0,Perr[bin]/1000);
            }
            st.GraphErrors(graph2,axfin);
            graph2 -> SetLineColor(kRed);
            graph2 -> SetMarkerColor(kRed);
            
        }
        else if(f==2){
            rep(bin,15){
                graph3 -> SetPoint(bin,(bin-7)*5,Pres[bin]/1000);
                graph3 -> SetPointError(bin,0,Perr[bin]/1000);
            }
            st.GraphErrors(graph3,axfin);
            graph3 -> SetLineColor(kGreen);
            graph3 -> SetMarkerColor(kGreen);
            
        }*/
    }

    graph1 -> Draw("APC");
    /*graph2 -> Draw("PC");
    graph3 -> Draw("PC");
    TLegend* legend = new TLegend(0.65, 0.65, 0.85, 0.85);
    legend -> AddEntry(graph1,"0kHz","l");
    legend -> AddEntry(graph2,"+30kHz","l");
    legend -> AddEntry(graph3,"-30kHz","l");
    //legend -> Draw();*/
}