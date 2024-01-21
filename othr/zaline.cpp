#include <iostream>
#include "../headers/setting.h"
using namespace std;
double pos[4] = {0,40,390,-233};
double xpower[4] = {257.5,334.8,336.2,276.4};
double ypower[4] = {278.8,329.4,335.5,265.5};
void zaline(){
    Setting st;
    st.dot_size = 0.6;
    st.markerstyle = 20;
    st.color = kBlue;
    //z軸に沿ったプロットを作る(x-yに関して座標固定or ピークの値)
    /*
    ※最初は最大値で議論
    初期位置：(269.2,270.8)
    日を跨いで
    初期位置：(257.5,278.8)
    40cm前：(334.8,329.4)
    390mm後：(336.2,335.5)
    233mm前：(276.4,265.5)
    */ 
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    axrange ax = {-300,500,200,400,0,1,";pos[mm];Power[uW]"};
    TGraph* graph = new TGraph;
    rep(i,4)graph -> SetPoint(i,pos[i],xpower[i]);
    st.Graph(graph,ax);
    graph -> Draw("AP");

}