//極めて簡便にrootfileの中身を確認するためだけのマクロ
#include <iostream>
#include "../headers/fitter.h"
using namespace std;

void checkroofile(){
    int i = 11;
    int j = 0;
    int p = 1;
    string filename = "peakfitdata" + to_string(i) + ".root";
    TFile* file = new TFile(filename.c_str());//読み取り専用なので鉄の意志でrecreateはしない
    string treename = "tree" + to_string(2*j+(p-1));
    TTree* tree = (TTree*)file->Get(treename.c_str());
    int bin;
    double pfit,dpfit,chindf;
    tree -> SetBranchAddress("bin",&bin);
    // tree -> SetBranchAdress("pfit",&pfit);
    // tree -> SetBranchAdress("dpfit",&dpfit);
    // tree -> SetBranchAdress("chiNDF",&chindf);
    // rep(num,100){
    //     tree -> GetEntry(num);
    //     cout << "======================" << endl;
    //     cout < "bin: " << bin << " pfit: " << pfit << endl;
    //     cout < "dpfit: " << dpfit << " chiNDF: " << chindf << endl;
    // }
}