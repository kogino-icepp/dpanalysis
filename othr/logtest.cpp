#include <iostream>
#include "../headers/fitter.h"

using namespace std;
void logtest(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.12,0.14,0.2,0.1);
    Setting st;
    st.dot_size=0.8;
    st.markerstyle=20;
    st.color = kGreen;
    TGraph* graph = new TGraph;
    TF1* area = new TF1("area","0.5",216,264);
    axrange axg = {1,200,0,1,0,1,"check_log;m_{DP};#chi"};
    st.Graph(graph,axg);
    graph -> SetPoint(0,250,0.5);
    graph -> SetPoint(1,1000,0.5);
    graph -> SetPoint(2,10,0.5);
    c1 -> SetLogx();
    graph -> Draw("AP");
    area -> Draw("same");
}