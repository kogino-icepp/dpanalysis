#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "../headers/mask_map.h"
using namespace std;

void read_map(){
    Mask ms;
    vector<vector<int>> vec = ms.maskmap;
    for(auto v:vec[1]){
        cout << v << endl;
    }
}