#include <iostream>
#include <random>
#include "../headers/fitter.h"
using namespace std;
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
double candpow[6] = {0,200,400,600,800,1000};
double candmass[3] = {220,240,260};
Color_t candcolor[3] = {kBlue,kRed,kGreen};
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
//仮説：f0がf0+delFに変わっただけでは？？
double F_sig_delta(double f,double f0,double P,double r,double delF){
    if(f+dnu*r<=f0+delF)return 0;
    else if(f+r*dNu>f0+delF && f-(1-r)*dNu<=f0+delF){
        return P*(F_nu(f+r*dNu,f0+delF)-F_nu(f0+delF,f0+delF));
    }
    else if(f-dnu*(1-r)>f0+delF){
        return P*(F_nu(f+r*dNu,f0+delF)-F_nu(f-(1-r)*dNu,f0+delF));
    }
    else return 0;
}
void fitsyserr(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);

    Setting st;
    st.dot_size = 0.8;
    st.markerstyle = 20;
    st.color = kBlue;
    st.lcolor = kBlue;
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_int_distribution<> rand100(0, 99);
    normal_distribution<double> normrand(0,100);
    //m_γ=220,240,260GHzにおいて周波数確度が-1~1kHzの不定性を持っている時のフィット精度を検証する
    double f0 = 220;//DM質量に対応する周波数[GHz]
    //TH1* testhist = new TH1D("testhist","test;Freq[GHz];",20,f0-dnu,f0+19*dnu);
    //TH1* testhist2 = new TH1D("testhist2","test;Freq[GHz];",20,f0-dnu,f0+19*dnu);//実験用のヒストグラム、質量をちょっとずらす
    //testhist2 -> SetLineColor(kRed);
    //フィットのエラーを記録した図;
    TGraphErrors* graph1 = new TGraphErrors;
    TGraphErrors* graph2 = new TGraphErrors;
    TGraphErrors* graph3 = new TGraphErrors;
    axrange ax = {-40,40,0.7,1.3,0,1,";#Delta [kHz];mean(P_{fit})/P_{given}"};
    st.GraphErrors(graph1,ax);
    st.GraphErrors(graph2,ax);
    st.GraphErrors(graph3,ax);
    graph1 -> SetMarkerColor(candcolor[0]);
    graph2 -> SetMarkerColor(candcolor[1]);
    graph3 -> SetMarkerColor(candcolor[2]);
    graph1 -> SetLineColor(candcolor[0]);
    graph2 -> SetLineColor(candcolor[1]);
    graph3 -> SetLineColor(candcolor[2]);
    TH1* signal1 = new TH1D("signal1",";Freq[GHz];Spec[a.u.]",20,f0-dnu,f0+19*dnu);
    TH1* signal2 = new TH1D("signal2",";Freq[GHz];Spec[a.u.]",20,f0-dnu,f0+19*dnu);
    TH1* signal3 = new TH1D("signal3",";Freq[GHz];Spec[a.u.]",20,f0-dnu,f0+19*dnu);
    //周波数220,240,260GHzに対してmassをビン内で変えながら(5kHz刻みで±35kHzまでいけるはず)ピークフィットし、真の値からのずれを評価する
    double delF = 30*pow(10,-6);
    
    rep(f,1){//大きな周波数の刻み(3パターン)
        for(int d=-7;d<8;d++){//massを8パターンにずらして実験する(d=0が0シフト)
            double Psum = 0;
            double Perrsum = 0;
            rep(ite,1){//イテレーションの回数(こんなに回数いる？まあ適度に増やしていくか)  
                TH1D* fithist = new TH1D("fithist",";;",20,f0-dnu,f0+19*dnu);
                TF1* pfunc = new TF1("pfunc","F_sig2(x,[0],0.934*[1],0.5)",candmass[f]-dnu,candmass[f]+19*dnu);
                pfunc -> FixParameter(0,candmass[f]);
                //ヒストグラム作成
                rep(bin,20){
                    fithist -> SetBinContent(bin,F_sig_delta(candmass[f]+(bin-1)*dnu,candmass[f],1000,0.5,d*5*pow(10,-6)));
                    fithist -> SetBinError(bin,100);
                }
                fithist -> Fit(pfunc);
                double pout = pfunc -> GetParameter(1);
                double perr = pfunc -> GetParError(1); 
                Psum+=pout;
                Perrsum += perr;
                delete fithist;
            }
            
            if(f==0)graph1 -> SetPoint(d+7,d*5,Psum/1000);
            else if(f==1)graph2 -> SetPoint(d+7,d*5,Psum/1000);
            else graph3 -> SetPoint(d+7,d*5,Psum/1000);

            if(f==0)graph1 -> SetPointError(d+7,0,Perrsum/1000);
            else if(f==1)graph2 -> SetPointError(d+7,0,Perrsum/1000);
            else graph3 -> SetPointError(d+7,0,Perrsum/1000);
        }
    }
    graph1 -> Draw("APL");
    graph2 -> Draw("PL");
    graph3 -> Draw("PL");
    /*rep(bin,20){
        signal1 -> SetBinContent(bin,F_sig_delta(f0+dnu*(bin-1),f0,1000,0.5,0));
        signal2 -> SetBinContent(bin,F_sig_delta(f0+dnu*(bin-1),f0,1000,0.5,delF));
        signal3 -> SetBinContent(bin,F_sig_delta(f0+dnu*(bin-1),f0,1000,0.5,-delF));
    }
    signal1 -> SetLineColor(kBlue);
    signal2 -> SetLineColor(kRed);
    signal3 -> SetLineColor(kGreen);
    st.Hist(signal1);
    signal1 -> Draw();
    signal2 -> Draw("same");
    signal3 -> Draw("same");
    TLegend *legend = new TLegend(0.55, 0.55, 0.75, 0.75); 
    legend->AddEntry(signal1, "#pm 0kHz", "l");
    legend->AddEntry(signal2, "+30kHz", "l");
    legend->AddEntry(signal3, "-30kHz", "l");
    legend -> Draw();
    //さらに3通りの周波数でこのプロットを作る
    
    //5000回まで100回ずつその収束性を確認していく
    prep(f,0,3){
        rep(p,6){
            double pvec[5000];
            double pevec[5000];
            double P = 0;
            double Perr = 0;
            TF1* pfunc = new TF1("pfunc","F_sig2(x,[0],0.95*[1],0.5)",candmass[f]-dnu,candmass[f]+19*dnu);
            pfunc -> FixParameter(0,candmass[f]);
            rep(ite,5000){
                TH1* fithist = new TH1D("fithist","",20,candmass[f]-dnu,candmass[f]+19*dnu);
                double pft,perr;
                rep(bin,20){
                    double freq = candmass[f]+(bin-1)*dnu;
                    double dp = normrand(mt);
                    fithist -> SetBinContent(bin,F_sig2(freq,candmass[f],candpow[p],0.5)+dp);
                    fithist -> SetBinError(bin,100);
                }
                fithist -> Fit(pfunc);
                pft = pfunc -> GetParameter(1);
                perr = pfunc -> GetParError(1);
                P += pft;
                Perr += perr;
                pvec[ite] = pft;
                pevec[ite] = perr;
                delete fithist;
            }
            P /= 5000;
            Perr /= 5000;
            double stdP = 0;
            rep(ite,5000)stdP += (pvec[ite]-P)*(pvec[ite]-P);
            stdP /= 5000;
            stdP = sqrt(stdP);
            if(f==0)graph1 -> SetPoint(p,candpow[p],stdP/Perr);
            else if(f==1)graph2 -> SetPoint(p,candpow[p],stdP/Perr);
            else if(f==2)graph3 -> SetPoint(p,candpow[p],stdP/Perr);
        }
        
        //else if(f==1) graph2 -> Draw("P");
        //else graph3 -> Draw("P");
        //delete graph;
    }
    rep(f,3){
        int pnum = 0;
        //TH1D * testhist = new TH1D("testhist",";;",100,200,1300);
        for(int ite=100;ite<5000;ite+=100){
            double P = 0;
            rep(k,ite){
                TH1* fithist = new TH1D("fithist","",20,candmass[f]-dnu,candmass[f]+19*dnu);
                TF1* pfunc = new TF1("pfunc","F_sig2(x,[0],0.94*[1],0.5)",candmass[f]-dnu,candmass[f]+19*dnu);
                pfunc -> FixParameter(0,candmass[f]);
                double pft,perr;
                //ここでフィットするヒストグラムを作る
                rep(bin,20){
                    double freq = candmass[f]+(bin-1)*dnu;
                    double dp = normrand(mt);
                    fithist -> SetBinContent(bin,F_sig2(freq,candmass[f],1000,0.5)+dp);
                    fithist -> SetBinError(bin,100);
                }
                fithist -> Fit(pfunc);
                pft = pfunc -> GetParameter(1);
                perr = pfunc -> GetParError(1);
                P += pft;
                delete fithist;
                //testhist -> Fill(pft);
            }
            if(f==0)graph1 -> SetPoint(pnum,ite,P/(ite*1000));
            else if(f==1)graph2 -> SetPoint(pnum,ite,P/(ite*1000));
            else if(f==2)graph3 -> SetPoint(pnum,ite,P/(ite*1000));
            pnum++;
        }
        //st.Hist(testhist);
        //testhist -> Draw();
    }
    
    TF1* gline = new TF1("gline","0.005*x",0,1100);
    gline -> SetLineColor(kBlack);
    gline -> SetLineStyle(kDashed);
    gline -> Draw("same");
    rep(bin,20){
        double freq = 220+(bin-1)*dnu;
        
        double perr = normrand(mt);
        testhist -> SetBinContent(bin,F_sig2(freq,220,1000,0.5)+perr);
        testhist -> SetBinError(bin,100);
        //testhist2 -> SetBinContent(bin,F_sig2(freq,220-1*pow(10,-6),1000,0.5));
    }
    testhist -> SetFillColor(kBlue);
    st.Hist(testhist);
    testhist -> Draw("E");
    TF1* pfunc = new TF1("pfunc","F_sig2(x,220,[0],0.5)",f0-dnu,f0+19*dnu);
    testhist -> Fit(pfunc);*/
    //このフィットを5000回回せばええの？
}