#include <iostream>
#include "../headers/fitter.h"
using namespace std;

void apdist(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    //c1 -> SetLogy();
    int nvec[10] = {3,18,130,768,2382,4481,4881,2994,816,62};
    int nvecp[10] = {5,28,477,2864,7174,10645,10386,5175,976,33};
    TH1D* hist = new TH1D("hist",";Score;Count",10,0,100);
    rep(i,10){
        hist -> SetBinContent(i+1,nvecp[i]);
        cout << i << " " << nvecp[i] << endl;
    }

    st.Hist(hist);
    hist -> Draw();
    hist -> Fit("gaus","I");
}