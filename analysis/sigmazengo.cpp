#include <iostream>
#include "../headers/setting.h"

void sigmazengo(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    Setting st;
    st.dot_size = 0.8;
    st.markerstyle = 20;
    st.color = kBlue;
    vector<double> sigzen = {0.0680209,0.0624807,0.0581624,0.0804053,0.0617756,0.100422,0.0658984,0.0629864,0.0755145,
    0.0545146,0.107824,0.0586616,0.0782594,0.146514,0.069492,0.0713242,0.0789449,0.0836776,0.137635,0.0752521,0.0726951,0.0664835,0.0805179,0.0874365};
    vector<double> sigkou = {0.0684595,0.078938,0.0775487,0.070623,0.0694093,0.0666186,0.0846847,0.0727667,0.059874,0.0639881,0.0642888,0.0663619,0.0963929,0.123084,
    0.0913639,0.120602,0.0896919,0.134783,0.109617,0.0942582,0.104724,0.0815718,0.124269,0.099359};

    cout << (int)sigzen.size() << " : " << (int)sigkou.size() << endl;
    TGraph* graph = new TGraph;
    axrange ax = {0,25,-0.1,0.1,0,1,"the_diff;measure;dsigma[K]"};
    rep(i,24){
        graph -> SetPoint(i,i+1,sigkou[i]-sigzen[i]);
    }
    TGraph* rgraph = new TGraph;
    st.Graph(graph,ax);
    graph -> Draw("AP");
}