#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <TROOT.h>
#include "../headers/mask_map.h"
#include "../headers/setting.h"
using namespace std;
#define rep(i,n) for(int i=0;i<n;i++)
#define prep(i,s,f) for(int i=s;i<f;i++)
const int nbin=32767;
const int sb = 2585;//切り出してくるビンの最初
const int cb = 15725;//探索範囲の中点,ここを境にデータの振幅が変わっているデータがある
const int fb = 28865;//切り出してくるビンの最後
vector<string> XFFT = {"lsbo","lsbi","usbi","usbo"};
void read_map(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    Setting st;
    st.dot_size = 0.8;
    st.markerstyle = 20;
    st.color = kRed;
    Mask ms;
    vector<vector<int>> vec = ms.maskmap;
    bool hantei[4][nbin];//ここを前後15ビンまでマークするようにする
    rep(xfft,4){
        rep(bin,nbin)hantei[xfft][bin] = false;
    }
    rep(xfft,4){
        for(auto v:vec[xfft]){
            for(int d=-15;d<15;d++)hantei[xfft][v+d] = true;//探索できないところは全て潰した
        }
    }
    int x = 0;
    rep(bin,nbin){
        if(hantei[1][bin] && hantei[1][bin+512] && hantei[1][bin-512] && hantei[1][bin-1024]){
            cout << bin << endl;
            x++;
        }
    }
    cout << "x : " << x << endl;
    /*TGraph* graph1 = new TGraph;
    TGraph* graph2 = new TGraph;
    TGraph* graph3 = new TGraph;
    TGraph* graph4 = new TGraph;
    axrange ax = {0,nbin,-1,1,0,1,"LSBO;Bin;Judge"};
    st.Graph(graph1,ax);
    st.Graph(graph2,ax);
    st.Graph(graph3,ax);
    st.Graph(graph4,ax);
    graph1 -> SetMarkerColor(kBlue);
    graph2 -> SetMarkerColor(kRed);
    graph3 -> SetMarkerColor(kGreen);
    graph4 -> SetMarkerColor(kMagenta);
    int num = 0;
    for(auto v:vec[1]){
        graph1 -> SetPoint(num,v,0.75);
        graph2 -> SetPoint(num,v-512,1);
        graph3 -> SetPoint(num,v+512,0.5);
        graph4 -> SetPoint(num,v+1024,0.25);
        num++;
    }
    //被っている点を具体的に列挙してみる
    graph1 -> Draw("AP");
    graph2 -> Draw("P");
    graph3 -> Draw("P");
    graph4 -> Draw("P");*/
}