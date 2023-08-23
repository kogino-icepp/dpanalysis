#include <TFile.h>
#include <TTree.h>
#include <iostream>
using namespace std;
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

void chiroot_check() {
    filesystem::path path=filesystem::current_path();
    string root_dir = "/Users/oginokyousuke/data/root_file/";
    filesystem::current_path(root_dir);
    const char* filename = "fit_res1_0_1.root"; // Specify the path to your ROOT file
    const char* treeName = "tree1"; // Specify the name of the TTree
    Int_t numEntriesToShow = 10; // Specify the number of entries to show

    PrintEntryInfo(filename, treeName, numEntriesToShow);

    
}
