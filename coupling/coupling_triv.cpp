#include <iostream>
#include <TROOT.h>
using namespace std;
double cv=3.0*pow(10,11); 
int freq=240; 
double lam=cv/(freq*pow(10,9)); 
double w0=2.67;
double z0=M_PI*w0*w0/lam;
double k=2*M_PI/lam;
double west1(double z){
    double w=w0*sqrt(1+pow(z/z0,2));
    return w;
}
double west2(double z,double kai){
    double w=kai*w0*sqrt(1+pow(z/z0,2));
    return w;
}
double R1(double z){
    double r=0;
    r=z+(1/z)*pow((M_PI*w0*w0/lam),2);
    return r;
}
double R2(double z,double kai){
    double r=0;
    r=z+(1/z)*pow((M_PI*w0*w0*kai*kai/lam),2);
    return r;
}
double eff1(double z){
    double e=0;
    e=sqrt(2/(M_PI*west1(z)*west1(z)));
    return e;
}
double eff2(double z,double kai){
    double e=0;
    e=sqrt(2/(M_PI*west2(z,kai)*west2(z,kai)));
    return e;
}
complex<double> gauss1(double x,double y, double z){
    complex<double> kata(-(x*x+y*y)/(west1(z)*west1(z)),-k*z-M_PI*(x*x+y*y)/(lam*R1(z)));
    complex<double> g=exp(kata);
    return g;
}
complex<double> gauss2(double x,double y,double z,double kai){
    complex<double> kata(-(x*x+y*y)/(west2(z,kai)*west2(z,kai)),-k*z-M_PI*(x*x+y*y)/(lam*R2(z,kai)));
    complex<double> g=exp(kata);
    return g;
}
//平行にずれる方向と振れる角度が同じ場合のカップリングの計算
complex<double> eta_all(double kai){
    complex<double> eta(0);
    for(int i=-1000;i<1001;i++){
        for(int j=-1000;j<1001;j++){
            if(i*i+j*j<=1000*1000){
                double z=sqrt(1500*1500-i*i-j*j);
                //double inpro=(i*xrot2(i,z,theta,dx)+j*j+zrot2(i,z,theta,dz)*z)/(norm1*norm2);
                eta+=eff1(z)*eff2(z,kai)*conj(gauss1(i,j,z))*gauss2(i,j,z,kai);
            }
            
        }
    }
    return eta;
}

//xまたはy方向(多分xの方が意味がある)ずらした時にθ方向に関してどの程度の幅で推移するか
//本処理パート
void coupling_triv()
{
    Double_t xlo=-1;
    Double_t xhi=7;
    Double_t ylo=0;
    Double_t yhi=1;
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,800,600);
    TH1F *frame=gPad->DrawFrame(xlo,ylo,xhi,yhi);
    TGraph* graph1=new TGraph;
    TGraph* graph2=new TGraph;
    /*TGraph* graph3=new TGraph;
    TGraph* graph4=new TGraph;
    TGraph* graph5=new TGraph;*/
    TLegend *leg=new TLegend(0.8,0.68,0.99,0.78);
    leg->SetTextSize(0.02);
    leg->SetTextFont(42);
    leg->SetFillStyle(0);
    leg->AddEntry(graph1,"calc","lp");
    leg->AddEntry(graph2,"0.0","lp");
    /*leg->AddEntry(graph3,"0.5","lp");
    leg->AddEntry(graph4,"1.0","lp");
    leg->AddEntry(graph5,"2.5","lp");*/
    
    const Int_t n=50;
    int count =0;
    
    for(int i=0;i<n;i++){
        double x=0.1*i;
        double z=1500;
        double kai=(M_PI*west1(z)*west1(z)/lam)*((1/R1(z))-(1/R2(z,x)));
        cout << kai << endl;
        double y=pow(abs(eta_all(x)),2);
        double y2=4/((x+(1/x))*(x+(1/x))+x*x*kai*kai);
        /*double y3=4/((x+(1/x))*(x+(1/x))+0.5*0.5*x*x);
        double y4=4/((x+(1/x))*(x+(1/x))+x*x);
        double y5=4/((x+(1/x))*(x+(1/x))+2.5*2.5*x*x);*/
        if(x>0.3)graph1->SetPoint(i,x,y);
        graph2->SetPoint(i,x,y2);
        /*graph3->SetPoint(i,x,y3);
        graph4->SetPoint(i,x,y4);
        graph5->SetPoint(i,x,y5);*/
    }
    gStyle->SetPalette(1);
    graph1->SetMarkerStyle(20);
    graph1->SetMarkerSize(0.7);
    graph1->SetMarkerColor(kBlue);
    graph1->SetTitle("K;wa/wb;coupling");
    graph1->SetLineColor(kBlue);
    graph2->SetLineColor(kMagenta);
    /*graph3->SetLineColor(kGreen);
    graph4->SetLineColor(kRed);
    graph5->SetLineColor(kCyan);
    TF1 *f1 = new TF1("f1","[0]*exp(-pow((x-[1]),2)/(2*pow([2],2)))",-0.2,0.2);
    f1->SetParameter(0,0.9);
    f1->SetParameter(1,0);
    f1->SetParameter(2,2);
    graph1->Fit("f1");
    gStyle->SetOptFit(1111);*/
    graph1->Draw("AP");
    graph2->Draw("CP");
    /*graph3->Draw("CP");
    graph4->Draw("CP");*/
    leg->Draw("same");
    /*graph5->Draw("CP");*/
}