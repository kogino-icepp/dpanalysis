#include <iostream>
#include <TROOT.h>
using namespace std;
double cv=3.0*pow(10,11); 
int freq=240; 
double lam=cv/(freq*pow(10,9)); 
double w0=2.67;
double z0=M_PI*w0*w0/lam;
double k=2*M_PI/lam;
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
//平行移動と角度移動が同じ方向の場合の変換
double xrot2(double x,double z,double theta,double dx){
    double rtn=0;
    rtn=-sin(theta)*z+cos(theta)*x-dx;
    return rtn;
}
double zrot2(double x,double z,double theta,double dz){
    double rtn=0;
    rtn=cos(theta)*z+sin(theta)*x-dz;
    return rtn;
}
//平行移動と角度移動が直交している場合の変換
double yrot(double y,double z,double theta){
    double rtn=0;
    rtn=-sin(theta)*z+cos(theta)*y;
    return rtn;
}

complex<double> gauss(double x,double y, double z){
    complex<double> kata(-(x*x+y*y)/(west(z)*west(z)),-k*z-M_PI*(x*x+y*y)/(lam*R(z)));
    complex<double> g=exp(kata);
    return g;
}
//平行にずれる方向と振れる角度が同じ場合のカップリングの計算
complex<double> eta_all(double dx,double dz,double theta){
    complex<double> eta(0);
    for(int i=-350;i<351;i++){
        for(int j=-350;j<351;j++){
            if(i*i+j+j<=350*350){
                double z=sqrt(1500*1500-i*i-j*j);
                double norm1=1500;
                double norm2=sqrt(xrot2(i,z,theta,dx)*xrot2(i,z,theta,dx)+j*j+zrot2(i,z,theta,dz)*zrot2(i,z,theta,dz));
                double inpro=(i*xrot2(i,z,theta,dx)+j*j+zrot2(i,z,theta,dz)*z)/(norm1*norm2);
                eta+=eff(z)*eff(zrot2(i,z,theta,dz))*conj(gauss(i,j,z))*gauss(xrot2(i,z,theta,dx),j,zrot2(i,z,theta,dz))*inpro;
        }
        }
    }
    return eta;
}
//平行移動の方向と角度方向が直交している場合
complex<double> eta_all2(double dx,double dz,double theta){
    complex<double> eta(0);
    for(int i=-350;i<351;i++){
        for(int j=-350;j<351;j++){
            double z=sqrt(1500*1500-i*i-j*j);
            double norm1=1500;
            double norm2=sqrt((i-dx)*(i-dx)+yrot(j,z,theta)*yrot(j,z,theta)+zrot2(j,z,theta,dz)*zrot2(j,z,theta,dz));
            double inpro=(i*(i-dx)+j*yrot(j,z,theta)+zrot2(j,z,theta,dz)*z)/(norm1*norm2);
            eta+=eff(z)*eff(zrot2(j,z,theta,dz))*conj(gauss(i,j,z))*gauss(i-dx,yrot(j,z,theta),zrot2(j,z,theta,dz))*inpro;
        }
    }
    return eta;
}
//xまたはy方向(多分xの方が意味がある)ずらした時にθ方向に関してどの程度の幅で推移するか
//本処理パート
void coupling_3d()
{
    //TGraph2D* graph=new TGraph2D;
    TGraph* graph1=new TGraph;
    TGraph* graph2=new TGraph;
    TGraph* graph3=new TGraph;
    TGraph* graph4=new TGraph;
    Double_t xlo=-0.2;
    Double_t xhi=0.2;
    Double_t ylo=0;
    Double_t yhi=1;
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,800,600);
    TH1F *frame=gPad->DrawFrame(xlo,ylo,xhi,yhi);
    const Int_t n=200;
    int count =0;
    //角度とxをずらしながらカップリングの推移を見ていく
    /*for(int i=0;i<n;i++){
        double x=0.1*i-10;
        for(int j=0;j<n;j++){
            double y=0.005*j-0.5;
            double z=pow(abs(eta_all(x,0,y)),8);
            graph->SetPoint(1000*i+j,x,y,z);
            cout << count << endl;
            count+=1;
        }
    }*/
    double wg=(2/pow(west(1500),2));
    //ウエストの違いがフェータルに効いてくるのでちゃんと調べないといけないのでは？
    double thetat=(lam/M_PI)*sqrt(wg);
    double thetaw=(lam/M_PI)*sqrt(2)/w0;
    for(int i=0;i<n;i++){
        double x=0.001*i-0.1;
        double y=pow(abs(eta_all(0,0,x)),2);
        double y2=exp(-2*x*x/pow(thetat,2));
        //double y3=exp(-2*(x*x/(2*w0*w0)));
        double y4=exp(-2*x*x/pow(thetaw,2));
        graph1->SetPoint(i,x,y);
        graph2->SetPoint(i,x,y2);
        //graph3->SetPoint(i,x,y3);
        graph4->SetPoint(i,x,y4);
        cout << count << endl;
        count+=1;
    }

    //graph->SetTitle("Coupling_x_theta;x[mm];theta[rad];coupling");
    graph1->SetTitle("coupling_theta;theta[rad];coupling");
    graph1->SetMarkerStyle(20);
    graph1->SetMarkerSize(0.7);
    graph1->SetMarkerColor(kBlue);
    graph2->SetLineColor(kRed);
    //graph3->SetLineColor(kMagenta);
    graph4->SetLineColor(kCyan);
    //graph->Draw("surf1");
    graph1->Draw("AP");
    graph2->Draw("CP");
    //graph3->Draw("CP");
    graph4->Draw("CP");
}