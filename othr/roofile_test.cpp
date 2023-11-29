#include <iostream>
#include "../headers/fitter.h"
#include <random>
#define rep(i,n) for(int i=0;i<n;i++)
random_device rnd;     // 非決定的な乱数生成器を生成
mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
uniform_int_distribution<> rand100(0, 99);
normal_distribution<double> normrand(0,100);
void MakeTree(TFile* file,int i){
    file -> cd();
    string treename = "testtree" + to_string(i);
    TTree* tree = new TTree(treename.c_str(),treename.c_str());
    int testnum;
    tree -> Branch("a",&testnum,"a/I");
    rep(j,100){
        testnum = normrand(mt);
        tree -> Fill();
    }
}
void roofile_test(){
    TFile* file = new TFile("testfile.root","recreate");
    rep(i,10){
        MakeTree(file,i);
    }
    file -> Write();
    file -> Close();
}