#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include "../headers/setting.h"
using namespace std;
Double_t chiF_free(double x,double p0,double k,double p1){
    return p0*TMath::Gamma(k/2,x/(2*p1));
}
Double_t chiF_freefit(double x,double p0,double p1,double k,double bin){
    return (chiF_free((x+bin/2),p0,k,p1)-chiF_free((x-bin/2),p0,k,p1));
}
Setting st;
void PrintEntryInfo(const char* filename, const char* treeName, Int_t numEntriesToShow) {
    // Open the ROOT file
    TFile* file = TFile::Open(filename);
    if (!file) {
        cerr << "Error opening file " << filename << std::endl;
        return;
    }
    // Get the TTree
    TTree* tree = dynamic_cast<TTree*>(file->Get(treeName));
    if (!tree) {
        cerr << "Error getting tree " << treeName << " from file " << filename << std::endl;
        file->Close();
        return;
    }
    // Get the number of entries in the tree
    Int_t numEntries = tree->GetEntries();
    // Determine the number of entries to display
    Int_t numToShow = (numEntriesToShow <= 0) ? numEntries : std::min(numEntriesToShow, numEntries);
    // Loop over the entries and display information
    for (Int_t entry = 0; entry < numToShow; ++entry) {
        tree->GetEntry(entry);
        // Print entry information here
        // For example, print the values of specific branches
        Double_t a, b, c, chi;
        tree->SetBranchAddress("a", &a);
        tree->SetBranchAddress("b", &b);
        tree->SetBranchAddress("c", &c);
        tree->SetBranchAddress("chi", &chi);
        //tree->SetBranchAddress("ndf", &ndf);
        // Print the values for this entry
        cout << "Entry " << entry << ": a = " << a << ", b = " << b << ", c = " << c
                << ", chi = " << chi << endl;
    }
    // Close the file
    file->Close();
}
void check_chi(const char* filename,const char* treeName){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.15,0.1,0.2,0.1);
    TFile* file = TFile::Open(filename);
    if (!file) {
        cerr << "Error opening file " << filename << std::endl;
        return;
    }
    // Get the TTree
    TTree* tree = dynamic_cast<TTree*>(file->Get(treeName));
    if (!tree) {
        cerr << "Error getting tree " << treeName << " from file " << filename << std::endl;
        file->Close();
        return;
    }
    Int_t numEntries = tree->GetEntries();
    int binnum = 100;
    int binmin = 0;
    int binmax = 10;
    double binhaba = (binmax-binmin)/binnum;
    TH1D* chi_hist = new TH1D("chi_hist","chi_hist;Chi2/NDF;Count",100,0,10);
    vector<double> vchi;
    vector<vector<double>> vpara(3);
    for (Int_t entry = 0; entry < numEntries; ++entry) {
        tree->GetEntry(entry);
        // Print entry information here
        // For example, print the values of specific branches
        Double_t a, b, c, chi;
        int bin;
        tree->SetBranchAddress("a", &a);
        tree->SetBranchAddress("b", &b);
        tree->SetBranchAddress("c", &c);
        tree->SetBranchAddress("chi", &chi);
        tree->SetBranchAddress("bin", &bin);
        cout << bin << endl;
        chi_hist -> Fill(chi);
    }
    //c1 -> SetLogy();
    st.Hist(chi_hist);
    chi_hist -> Draw();
    
    TF1* chifit = new TF1("chifit","chiF_freefit(x,[0],[1],27,0.1)");
    TF1* chifit2 = new TF1("chifit2","chiF_freefit(x,[0],[1],27,0.1)+chiF_freefit(x,[2],[3],27,0.1)");
    int maxbin = chi_hist -> GetMaximumBin();
    chifit -> FixParameter(0,numEntries);
    chifit -> SetParameter(1,maxbin*0.1/(25.0));
    cout << maxbin << endl;
    chi_hist -> Fit(chifit,"E","",0,10);
    chifit -> Draw("same");
    double p1 = chifit -> GetParameter(1);
    cout << fixed;
    cout << setprecision(30) << "alpha : " << 1-chiF_free(1.5,1,27,p1) << endl;
}

void chiroot_check() {
    filesystem::path path=filesystem::current_path();
    string root_dir = "/Users/oginokyousuke/data/root_file/";
    filesystem::current_path(root_dir);
    const char* filename = "fit_res2_2_2.root"; // Specify the path to your ROOT file
    const char* treeName = "tree2"; // Specify the name of the TTree
    Int_t numEntriesToShow = 10; // Specify the number of entries to show
    //PrintEntryInfo(filename, treeName, numEntriesToShow);
    TFile* file = TFile::Open(filename);
    TTree* tree = dynamic_cast<TTree*>(file->Get(treeName));
    check_chi(filename,treeName);
    
}