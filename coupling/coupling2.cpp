//2023/12/08追記　これ途中で諦めてるっぽい
#include <iostream>
#include <TROOT.h>
using namespace std;
double cv=3.0*pow(10,11); 
int freq=240; 
double lam=cv/(freq*pow(10,9)); 
double w0=2.67;
double z0=M_PI*w0*w0/lam;
double k=2*M_PI/lam;
double Rm=1500;
double west(double z){
    double w=w0*sqrt(1+pow(z/z0,2));
    return w;
}
double R(double z){
    double r=0;
    r=z+(1/z)*pow((M_PI*w0*w0/lam),2);
    return r;
}
double eff(double z){
    double e=0;
    e=sqrt(2/(M_PI*west(z)*west(z)));
    return e;
}
complex<double> gauss(double x,double y,double z){
    complex<double> kata(-(x*x+y*y)/(west(z)*west(z)),-k*z-M_PI*(x*x+y*y)/(lam*R(z)));
    complex<double> g=exp(kata);
    return g;
}
//メッシュ上の一点から出たビームがどこに帰ってくるのかを変換して計算,カーテシァンで表示(ミラー位置とメッシュ上での座標がパラメータ)
//ヒットした時の軸からの距離
double rhit(double x,double y,double dz){
    double rtn=0;
    double r0=sqrt(pow(x,2)+pow(y,2));
    double ratio=z0/r0;
    rtn=(dz*ratio)*dz+sqrt(pow(dz,2)*pow(ratio,2)-(1+pow(ratio,2))*(pow(dz,2)-pow(Rm,2)));
    return rtn;
}
//ヒットした光がz=0のどの位置に返ってくるか(x座標について)
double xref(double x,double y,double dz){
    double rtn=0;
    double r0=sqrt(pow(x,2)+pow(y,2));
    if (rhit>350){
        return 1000;
    }
    else{
        double xhit=rhit(x,y,dz)*(x/sqrt(pow(x,2)+pow(y,2)));
        double zhit=z0*rhit(x,y,dz)/r0;
        double dxb=xhit/zhit;
        double dxa=-(2/Rm)*xhit+dxb;
        rtn=xhit-dxa*zhit;
        return rtn;
    }
}
//返ってくる光のy座標
double yref(double x,double y,double dz){
    double rtn=0;
    double r0=sqrt(pow(x,2)+pow(y,2));
    if (rhit>350){
        return 1000;
    }
    else{
        double yhit=rhit(x,y,dz)*(y/sqrt(pow(x,2)+pow(y,2)));
        double zhit=z0*rhit(x,y,dz)/r0;
        double dyb=yhit/zhit;
        double dya=-(2/Rm)*yhit+dyb;
        rtn=yhit-dya*zhit;
        return rtn;
    }
}
/*ここまで定義したはいいけど波面がミラーでどう跳ね返るかを本当は計算すべきでは？詰んだンゴ*/
double beampower()
//受信機ビームウエスト位置で数値積分())
complex<double> etaz(double dz){
    complex<double> eta(0);
    for(int i=0;i<100;i++){
        for(int j=0;j<100;j++){
            double meshx=-10+0.2*i;
            double meshy=-10+0.2*j;
            double ex=conj(gauss());
        }
    }
}
void coupling2()
{
    TGraph *g=new TGraph;
    for(int i=0;i<200;i++){
        double z=-1000+10*i;
        double c=pow(abs(etaz(z)),2);
        graph->SetPoint(i,z,c);
    }
    graph->Draw();
}