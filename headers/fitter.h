#include <iostream>
#include <random>
#include "p_func.h"
#include "setting.h"
#include <TMath.h>
#include <TROOT.h>
using namespace std;
#define rep(i,n) for(int i=0;i<n;i++)
#define prep(i,m,n) for(int i=m;i<n;i++)
#define rrep(i,n) for(int i=n;i>=0;i--)
#define rprep(i,m,n) for(int i=m;i>n;i--)
const int dbin = 30;
typedef long long ll;
string savedirt = "/Users/oginokyousuke/data/chi_hist/";
Setting st;
class Fitter{
    public:
    void GetThreeParameter(TGraphErrors* graph,TF1* f,int p0,int p1,int p2){
        double x0 = graph -> GetPointX(p0);
        double y0 = graph -> GetPointY(p0);
        double x1 = graph -> GetPointX(p1);
        double y1 = graph -> GetPointY(p1);
        double x2 = graph -> GetPointX(p2);
        double y2 = graph -> GetPointY(p2);
        
        double m = (y2-y0)/(x2-x0);
        double n = (y0*x2-y2*x0)/(x2-x0);
        double a = (y1-m*x1-n)/((x1-x0)*(x1-x2));
        double b = (x0+x2-(m/a))/2;
        double c = a*x0*x2+n-a*b*b;
        f -> SetParameter(0,a);
        f -> SetParameter(1,b);
        f -> SetParameter(2,c);
    }
    //点を与えたらパラメータを返してくれる関数
    vector<double> GetQuadPara(vector<pair<double,double>> points){
        vector<double> ans(3);
        double x0 = points[0].first;
        double y0 = points[0].second;
        double x1 = points[1].first;
        double y1 = points[1].second;
        double x2 = points[2].first;
        double y2 = points[2].second;
        double m = (y2-y0)/(x2-x0);
        double n = (y0*x2-y2*x0)/(x2-x0);
        double a = (y1-m*x1-n)/((x1-x0)*(x1-x2));
        double b = (x0+x2-(m/a))/2;
        double c = a*x0*x2+n-a*b*b;
        ans[0] = a;
        ans[1] = b;
        ans[2] = c;
        return ans;
    }
    void syoki_para(TGraphErrors* graph,TF1* f,int bin){
        //y=a(x-x0)(x-x2)+mx+nが(x1,y1)を通る
        //取得するパラメータはy=a(x-b)^2+cの方
        double x0 = graph -> GetPointX(bin);
        double y0 = graph -> GetPointY(bin);
        double x1 = graph -> GetPointX(bin+dbin/2);
        double y1 = graph -> GetPointY(bin+dbin/2);
        double x2 = graph -> GetPointX(bin+dbin-1);
        double y2 = graph -> GetPointY(bin+dbin-1);
        /*cout << x0 << " " << x1 << " " << x2 << endl;
        cout << y0 << " " << y1 << " " << y2 << endl;*/
        
        double m = (y2-y0)/(x2-x0);
        double n = (y0*x2-y2*x0)/(x2-x0);
        double a = (y1-m*x1-n)/((x1-x0)*(x1-x2));
        double b = (x0+x2-(m/a))/2;
        double c = a*x0*x2+n-a*b*b;
        /*cout << "a = " << a << endl;
        cout << "b = " << b << endl;
        cout << "c = " << c << endl;*/
        f -> SetParameter(0,a);
        f -> SetParameter(1,b);
        f -> SetParameter(2,c);
    }
    void syoki_para_d(TGraphErrors* graph,TF1* f,int bin){
        vector<double> ans(3);
        //y=a(x-x0)(x-x2)+mx+nが(x1,y1)を通る
        //取得するパラメータはy=a+bx+cx^の方
        double x0 = graph -> GetPointX(bin);
        double y0 = graph -> GetPointY(bin);
        double x1 = graph -> GetPointX(bin+dbin/2);
        double y1 = graph -> GetPointY(bin+dbin/2);
        double x2 = graph -> GetPointX(bin+dbin-1);
        double y2 = graph -> GetPointY(bin+dbin-1);
        
        double m = (y2-y0)/(x2-x0);
        double n = (y0*x2-y2*x0)/(x2-x0);
        double a = (y1-m*x1-n)/((x1-x0)*(x1-x2));
        double b = m-a*(x0+x2);
        double c = a*x0*x2+n;
        f -> SetParameter(0,c);
        f -> SetParameter(1,b);
        f -> SetParameter(2,a);
    }
    //初期値を、ランダムな3点のデータ点を取ることで与える
    void syoki_rand1(TGraphErrors* graph,TF1* f){
        //cout << "first" << endl;
        random_device rnd;
        mt19937 mt(rnd());
        uniform_int_distribution<> rand_bin(0,dbin-1);
        //cout << "second" << endl;
        int p0,p1,p2;
        vector<int> p(3);
        set<int> st;
        while((int)st.size()<3){
            int x = rand_bin(mt);
            st.insert(x);
        }
        int bnum = 0;
        for(auto v:st){
            p[bnum]=v;
            bnum++;
        }
        /*p0 = rand_bin(mt);
        cout << "third" << endl;
        while(p0==p1) p1 = rand_bin(mt);
        while(p2==p0 || p2==p1) p2 = rand_bin(mt);*/

        //cout << "forth" << endl;
        double x0 = graph -> GetPointX(p[0]);
        double y0 = graph -> GetPointY(p[0]);
        double x1 = graph -> GetPointX(p[1]);
        double y1 = graph -> GetPointY(p[1]);
        double x2 = graph -> GetPointX(p[2]);
        double y2 = graph -> GetPointY(p[2]);
        
        double m = (y2-y0)/(x2-x0);
        double n = (y0*x2-y2*x0)/(x2-x0);
        double a = (y1-m*x1-n)/((x1-x0)*(x1-x2));
        double b = (x0+x2-(m/a))/2;
        double c = a*x0*x2+n-a*b*b;
        //cout << "third" << endl;
        f -> SetParameter(0,a);
        f -> SetParameter(1,b);
        f -> SetParameter(2,c);
        //cout << a << " " << b << " " << c << endl;
    }
    /*
    初期値をランダムな3つのパラメータを振ることで最適化を目指す
    課題：どの程度探す？(乱数で回数調整or一通り走査)
    */
    void syoki_rand2(TF1*f){
        random_device rnd;
        mt19937 mt(rnd());
        uniform_real_distribution<> rand0(-1,1);
        uniform_real_distribution<> rand1(-2,2);
        uniform_real_distribution<> rand2(0,1);
        f -> SetParameter(0,rand0(mt));
        f -> SetParameter(1,rand1(mt));
        f -> SetParameter(2,rand2(mt));
    }
    void rand_fit(TGraphErrors* graph,TF1* f,int ite,int fite,double fm,double fM,double &res){
        syoki_para(graph,f,0);
        rep(i,fite)graph -> Fit(f,"MQN","",fm,fM);
        double chi2 = f -> GetChisquare();
        double ndf = f -> GetNDF();
        double chimin = chi2/ndf;
        double chimin0 = chimin;
        vector<double> good_para(3);
        rep(i,3)good_para[i] = f -> GetParameter(i);
        
        rep(i,ite){
            syoki_rand1(graph,f);
            //cout << i << endl;
            //syoki_rand2(f);
            rep(j,fite)graph -> Fit(f,"MQN","",fm,fM);
            chi2 = f -> GetChisquare();
            ndf = f -> GetNDF();
            chi2 /= ndf;
            //cout << chi2 << endl;
            if(chi2 < chimin){
                chimin = chi2;
                rep(k,3)good_para[k] = f -> GetParameter(k);
            }
        }
        rep(i,3) f -> SetParameter(i,good_para[i]);
        res = chimin;
    }
    
    void rand_fit_test(TGraphErrors* graph,TF1* f,int fite,double fm,double fM,double chi_test,int &num){
        syoki_para(graph,f,0);
        double chi2,ndf;
        rep(i,fite){
            graph -> Fit(f,"MQN","",fm,fM);
            chi2 = f -> GetChisquare();
            ndf = f -> GetNDF();
            
        }
        //cout << chi2/ndf << endl;
        chi2 = f -> GetChisquare();
        ndf = f -> GetNDF();
        double chimin = chi2/ndf;
        if(chi_test!=chimin){
            num++;
            cout << "Yes" << endl;
            cout << chi_test << " : " << chimin << endl;
        }
        
    }
    Double_t MyFunction(double x,double p0,double p1){
        return p0*TMath::Exp(-x/2)*pow(x,(p1/2)-1)/(pow(2,p1/2)*TMath::Gamma(p1/2));
    }
    //きちんと元の縮尺に戻してからヒストに入れてホワイトノイズの散らばりを見る
    void FillHist(TF1* f,TGraphErrors* graph,TH1D* hist,double yscale,double &dym){
        double x,y,dy,y1;
        double dymax = -100;
        rep(i,dbin){
            x = graph -> GetPointX(i);
            y = graph -> GetPointY(i);
            y1 = f -> Eval(x);
            dy = y1 - y;
            dy *= yscale;
            hist -> Fill(dy);
            if(abs(dy)>=dymax){
                //cout << dy << endl;
                dymax = abs(dy);
            }
        }
        if(dym<dymax)dym = dymax;
        //cout << "dymax = " << dymax << endl;

    }
    void exfit(TGraphErrors* graph,TF1*f,double &res){
        //4パターンでフィットする
        TF1* tf = new TF1("tf","[0]*(x-[1])*(x-[1])+[2]",0,1);
        vector<vector<double>> para(4,vector<double>(4,0));
        rep(i,4){
            if(i%2==0)tf -> SetParameter(0,0.1);
            else if(i%2==1)tf -> SetParameter(0,-0.1);
            if(i/2==0)tf -> SetParameter(2,1);
            else if(i/2==0)tf -> SetParameter(2,0);
            tf -> SetParameter(1,0.5);
            rep(j,100)graph -> Fit(tf,"MQE","",0,1);
            double chi = tf -> GetChisquare();
            int ndf = tf -> GetNDF();
            para[i][3] = chi/ndf;
            para[i][0] = tf -> GetParameter(0);
            para[i][1] = tf -> GetParameter(1);
            para[i][2] = tf -> GetParameter(2);
        }
        int index = -1;
        double resx = 100;
        rep(i,4){
            if(para[i][3]<resx){
                index = i;
                resx = para[i][3];
                rep(j,3)f -> SetParameter(j,para[i][j]);
            }
        }
        res = resx;
    }
    void make_scale(TGraphErrors* graph,TGraph* mgraph,int sbin,double &yscale){
        //走査範囲のレンジ調査
        double xmin,xmax,ymin,ymax;
        ymin = pow(10,30);
        ymax = -200;
        double x1 = mgraph -> GetPointX(sbin);
        double x2 = mgraph -> GetPointX(sbin+dbin-1);
        xmin = min(x1,x2);
        xmax = max(x1,x2);
        double x,y;
        prep(bin,sbin,sbin+dbin){
            y = mgraph -> GetPointY(bin);
            //cout << y << endl;
            if(y<ymin)ymin = y;
            if(y>ymax)ymax = y;
        }
        //レンジに合わせてセットポイントを変える
        //x -> (x-x0)/xrange, y ->  (y-ymin)/yrange ??
        //cout << ymin << " " << ymax << endl;
        double xrange = xmax-xmin;
        double yrange = ymax-ymin;
        prep(bin,sbin,sbin+dbin){
            x = mgraph -> GetPointX(bin);
            y = mgraph -> GetPointY(bin);
            x = (x-xmin)/xrange;
            y = (y-ymin)/yrange;
            graph -> SetPoint(bin-sbin,x,y);
            
            //cout << x << " " << y << endl;
            graph -> SetPointError(bin-sbin,0,0.1/yrange);
        }
        yscale = yrange;
    }
    //基本はmake_scaleとほぼ同じだが元のスケールでのパラメータを再現できるよう最小値と最大値、xscale,yscaleを持ってもらう
    //
    void make_scale2(TGraph* mgraph,TGraphErrors* graph,int sbin,double &yMin,double &yscale){
        //走査範囲のレンジ調査
        double xmin,xmax,ymin,ymax;
        ymin = 10000;
        ymax = -200;
        double x1 = mgraph -> GetPointX(sbin);
        double x2 = mgraph -> GetPointX(sbin+dbin-1);
        if(!x1 || !x2)return;
        if(x1<x2){
            xmin = x1;
            xmax = x2;
        }
        else{
            xmin = x2;
            xmax = x1;
        }
        double x,y;
        prep(bin,sbin,sbin+dbin){
            y = mgraph -> GetPointY(bin);
            //cout << y << endl;
            if(y<ymin)ymin = y;
            if(y>ymax)ymax = y;
        }
        double xrange = xmax-xmin;
        double yrange = ymax-ymin;
        prep(bin,sbin,sbin+dbin){
            y = mgraph -> GetPointY(bin);
            x = mgraph -> GetPointX(bin);
            x = (x-xmin)/xrange;
            y = (y-ymin)/yrange;
            graph -> SetPoint(bin-sbin,x,y);
            //cout << x << " " << y << endl;
            /*if(y<ymin)ymin = y;
            if(y>ymax)ymax = y;*/
        }
        //レンジに合わせてセットポイントを変える
        //x -> (x-x0)/xrange, y ->  (y-ymin)/yrange ??
        
        
        yMin = ymin;
        yscale = yrange;
    }
    //xmin,ymin,xscale,yscaleをもとにパラメータを元のスケールに直す関数
    void rescale_para(double &a,double &b,double &c,double xmin,double ymin,double xscale,double yscale,TF1*f){
        //c = ;
        //b = ;
        //a = ;
        f -> SetParameter(0,a*yscale/(xscale*xscale));
        f -> SetParameter(1,b*xscale+xmin);
        f -> SetParameter(2,c*yscale+ymin);
    }
    vector<double> fparameters(TF1* f,int pnum){
        vector<double> para;
        rep(i,pnum){
            double p;
            p = f -> GetParameter(i);
            para.push_back(p);
        }
        return para;
    }
    //100回のイテレーションがconvite回分毎回収束しているかなんかいかにバラけているかを確認
    void rand_conv1(TGraphErrors* graph,TF1* f,int convite){
        map<double,int> mp;
        rep(i,convite){
            //rand_fit(graph,f,ite,fite,fm,fM, &res)
            double res;
            rand_fit(graph,f,100,10,0,1,res);
            mp[res]++;
        }
        if((int)mp.size()==1){
            cout << "Random_Fitting is convergent" << endl;
        }
        else{
            for(auto v:mp){
                cout << fixed;
                cout << setprecision(10) <<v.first << " : " << v.second << endl;
            }
            cout << "-------" << endl;
            //cout << "Random_Fitting is not convergent" << endl;
        }
    }
    //必ず合う桁の確認、桁が緩いものに関してリストアップ
    void rand_conv2(TGraphErrors* graph,TF1* f,int convite,TH1D* hist,int bin){
        map<double,int> mp;
        rep(i,convite){
            //rand_fit(graph,f,ite,fite,fm,fM, &res)
            double res;
            rand_fit(graph,f,100,10,0,1,res);
            mp[res]++;
        }
        
        //10ずつ元のdouble keyにかけていきmpのサイズが1になるギリギリを見る
        int keta = 0;
        rep(i,20){
            
            map<ll,int> ketamp;
            for(auto v:mp){
                double res = v.first;
                rep(j,i)res *= 10;
                //cout << i << ": " << res << endl;
                ketamp[res] += v.second;
            }
            
            if(ketamp.size()==1){
                if(i==19){
                    cout << "correct keta : " << 20 << endl;
                    break;
                }
                else keta++;
            }
            else{
                cout << "correct keta : " << keta << endl;
                if(keta<=5){
                    cout << bin << endl;
                    for(auto v: mp){
                        cout << v.first << " : " << v.second << endl;
                    }
                }
                break;
            }
        }
        hist -> Fill(keta);
    }
    void section_fit(TGraphErrors* graph,TF1* f,double &res){
        double cand = 1000000;
        for(int p0=0;p0<10;p0++){
            for(int p1=10;p1<20;p1++){
                for(int p2=20;p2<30;p2++){
                    GetThreeParameter(graph,f,p0,p1,p2);
                    rep(i,5)graph -> Fit(f,"EQN","",0,1);
                    double chi = f -> GetChisquare();
                    double ndf = f -> GetNDF();
                    cand = min(cand,chi/ndf);
                }
            }
        }
        res = cand;
    }
    
    void all_fit(TGraphErrors* graph,TF1*f,int ite,double &res){
        //0以上bnum+selnum未満まで走査してくれる
        int bnum = 27;
        int selnum = 3;
        vector<int> field;
        rep(i,bnum)field.push_back(0);
        rep(i,selnum)field.push_back(1);
        res = 100000;
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
            //このパートでフィットをする
            vector<pair<double,double>> points;
            rep(j,selnum){
                double x = graph -> GetPointX(ans[j]);
                double y = graph -> GetPointY(ans[j]);
                points.push_back({x,y});
            }
            vector<double> paras = GetQuadPara(points);
            //cout << ans[0] << " " << ans[1] << " " << ans[2] << endl;
            //cout << paras[0] << " " << paras[1] << " " << paras[2] << endl;
            //cout << "---------------" << endl;
            rep(j,selnum)f -> SetParameter(j,paras[j]);
            rep(j,ite)graph -> Fit(f,"MQN","",0,1);
            double chi2 = f -> GetChisquare();
            double ndf = f -> GetNDF();
            if(res>chi2/ndf){
                res = chi2/ndf;
                rep(j,selnum)minpara[j]=f->GetParameter(j);
            }
            
        }while (next_permutation(field.begin(), field.end()));
        rep(i,selnum)f->SetParameter(i,minpara[i]);
    }
    //ダイナミックレンジが大きいor物理的な意味のない縦軸のグラフをフィットする際
    void calfit(TGraphErrors* graph,TF1*f,int ite,double &sigma){
        vector<pair<double,double>> points(3);
        double x0,x1,x2,y0,y1,y2;
        x0 = graph -> GetPointX(0);
        y0 = graph -> GetPointY(0);
        x1 = graph -> GetPointX(1);
        y1 = graph -> GetPointY(1);
        x2 = graph -> GetPointX(2);
        y2 = graph -> GetPointY(2);
        points[0] = {x0,y0};
        points[1] = {x1,y1};
        points[2] = {x2,y2};
        vector<double> paras = GetQuadPara(points);
        f -> SetParameter(0,paras[0]);
        f -> SetParameter(1,paras[1]);
        f -> SetParameter(2,paras[2]);
        rep(ite,5)graph -> Fit(f,"QN","",0,1);
        //rep(ite,5)graph -> Fit(f,"EQN","",0,1);
        //rep(ite,5)graph -> Fit(f,"MQN","",0,1);
        double x,y,yf;
        sigma = 0;
        rep(bin,30){
            x = graph -> GetPointX(bin);
            y = graph -> GetPointY(bin);
            yf = f -> Eval(x);
            sigma += (y-yf)*(y-yf);
        }
        sigma /= 27;
        sigma = sqrt(sigma);
        //cout << "sigma : " << sigma  << endl;
        rep(bin,30)graph -> SetPointError(bin,0,sigma);
        rep(ite,5)graph -> Fit(f,"QN","",0,1);
        //rep(ite,5)graph -> Fit(f,"EQN","",0,1);
        //rep(ite,5)graph -> Fit(f,"MQN","",0,1);
        double chi = f -> GetChisquare();
        int ndf = f -> GetNDF();
        //sigmaも引数として渡して評価するのはあり
        
    }
    void allfit2(TGraphErrors* graph,TF1*f,int ite,double &res){
    //0以上bnum+selnum未満まで走査してくれる
        int bnum = 27;
        int selnum = 3;
        vector<int> field;
        rep(i,bnum)field.push_back(0);
        rep(i,selnum)field.push_back(1);
        vector<double> minpara(3,0);
        double resx = pow(10,30);
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
            vector<pair<double,double>> points;
            rep(j,selnum){
                double x = graph -> GetPointX(ans[j]);
                double y = graph -> GetPointY(ans[j]);
                points.push_back({x,y});
            }
            vector<double> paras = GetQuadPara(points);
            //cout << ans[0] << " " << ans[1] << " " << ans[2] << endl;
            //cout << paras[0] << " " << paras[1] << " " << paras[2] << endl;
            //cout << "---------------" << endl;
            rep(j,selnum)f -> SetParameter(j,paras[j]);
            rep(j,ite){
                graph -> Fit(f,"Q","",0,1);
            }
            double chi2 = f -> GetChisquare();
            double ndf = f -> GetNDF();
            
            if(ndf<5){
                cout << "Less NDF" << endl;
                exit(1);
            }
            if(resx>chi2/ndf){
                //cout << chi2 << " " << ndf << endl;
                resx = chi2/ndf;
                rep(j,selnum)minpara[j]=f->GetParameter(j);
            }
            //cout << ans[0] << " " << ans[1] << " " << ans[2] << endl;
        }while (next_permutation(field.begin(), field.end()));
        rep(i,3)f -> SetParameter(i,minpara[i]);
        res = resx;
    }
    //spgraphと二次関数からホワイトノイズを導出する関数
    //フィットで得られた結果のステータスを確認するための関数を作りたい(結果が収束しているかなど)
};
class CheckData{
    public:
    void SaveTwoGraph(TGraphErrors* graph,TF1* f1,TF1* f2,string title,TCanvas* c1){
        axrange axstg = {0,1,0,1,0,1,"rand_fit vs section_fit;xscale[a.u];yscale[a.u]"};
        
        st.GraphErrors(graph,axstg);
        graph -> SetMarkerColor(kGreen);
        graph -> SetMarkerStyle(20);
        graph -> SetMarkerSize(0.8);
        graph -> Draw("AP");
        f1 -> SetLineColor(kBlue);
        f1 -> Draw("same");
        f2 -> SetLineColor(kRed);
        f2 -> Draw("same");
        double p1[3],p2[3];
        rep(i,3){
            p1[i] = f1 -> GetParameter(i);
            p2[i] = f2 -> GetParameter(i);
            cout << "p" << i << " : " << p1[i] << " <==> " << p2[i] << endl; 
        }
        filesystem::current_path(savedirt);
        c1 -> SaveAs(title.c_str());
    } 
};
class FitFunc{
    public:
    double v_conv(double f,double f0){
        double c = 3.0*pow(10,8);
        double rtn=c*sqrt(1-((f0/f)*(f0/f)));
        return rtn;
    }
    double F_nu(double f,double f0){
        double rtn;
        double v0=220000.0;
        double vE=v0;
        double vc=v0;
        double v=v_conv(f,f0);
        double p_kata=(v+vE)/v0;
        double m_kata=(v-vE)/v0;
        rtn = (vc/(2*sqrt(M_PI)*vE))*(exp(-(p_kata*p_kata))-exp(-(m_kata*m_kata)));
        rtn += 0.5*(erf(p_kata)+erf(m_kata));
        //rtn *= (c*f0*f0)/(f*f*f*sqrt(1-(f0/f)*(f0/f)));
        return rtn;
    }
    double F_sig2(double f,double f0,double P,double r){
        double dNu=0.0000885;//周波数分解能[GHz]
        double dnu=0.000076296;//ビン幅[GHz]
        if(f+dnu*r<=f0)return 0;
        else if(f+r*dNu>f0 && f-(1-r)*dNu<=f0){
            return P*(F_nu(f+r*dNu,f0)-F_nu(f0,f0));
        }
        else if(f-dnu*(1-r)>f0){
            return P*(F_nu(f+r*dNu,f0)-F_nu(f-(1-r)*dNu,f0));
        }
        else return 0;
    }
    double F_sig_delta(double f,double f0,double P,double r,double delF){
        double omega0 = f0+delF;
        double dNu=0.0000885;//周波数分解能[GHz]
        if((f-f0)-r*dNu-delF<=0 && f>omega0-r*dNu){
            return P*(F_nu(f+r*dNu,omega0)-F_nu(omega0,omega0));
        }
        else if((f-f0)-r*dNu-delF>0){
            return P*(F_nu(f+r*dNu,omega0)-F_nu(f-r*dNu,omega0));
        }
        else return 0;
    }
    double F_sigscale(double x,double P,double r,double fmin,double delF){
        double bin = x/0.0344827586;
        double dnu=0.000076296;//ビン幅[GHz]
        double f = fmin+bin*dnu;//ここをxとfminの関数に変える
        double f0 = fmin+10*dnu;//ここは自明ではあるがx0からf0へ
        return F_sig_delta(f,f0,P,r,delF);
    }
};