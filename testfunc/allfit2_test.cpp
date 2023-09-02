#include <iostream>
#include <vector>
#include <algorithm>
#define rep(i,n) for(int i=0;i<n;i++)

void allfit2_test(){
//0以上bnum+selnum未満まで走査してくれる
    int bnum = 27;
    int selnum = 3;
    vector<int> field;
    rep(i,bnum)field.push_back(0);
    rep(i,selnum)field.push_back(1);
    vector<double> minpara(3,0);
    do{
        vector<int> ans;
        int cand_num = 0;
        rep(i,field.size()){
            if(field[i]==0)cand_num++;
            if(field[i]==1){
                ans.push_back(cand_num);
                cand_num++;
            }
        }
        bool hantei = true;
        rep(j,3){
            if(ans[j]>=10 && ans[j]<20){
                hantei = false;
                break;
            }
        }
        if(!hantei)continue;
        cout << ans[0] << " " << ans[1] << " " << ans[2] << endl;
    }while (next_permutation(field.begin(), field.end()));
}