#include <TROOT.h>
#include <iostream>
using namespace std;
#define rep(i,n) for(int i=0;i<n;i++)
struct axrange{
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    string title;
};

class Setting{
public:
    double dot_size;
    double markerstyle;
    Color_t color;
    Color_t lcolor;
public:
    void GetLimit(TGraph2D* graph,axrange& ax){
        ax.xmax = graph -> GetXmax();
        ax.xmin = graph -> GetXmin();
        ax.ymax = graph -> GetYmax();
        ax.ymin = graph -> GetYmin();
        ax.zmax = graph -> GetZmax();
        ax.zmin = graph -> GetZmin();
    }
    void Graph(TGraph* graph,axrange& ax){
        graph -> SetMaximum(ax.ymax);
        graph -> SetMinimum(ax.ymin);
        graph -> GetXaxis() -> SetLimits(ax.xmin,ax.xmax);
        graph -> GetXaxis() -> SetLabelSize(0.08);
        graph -> GetXaxis() -> SetTitleOffset(0.9);
        graph -> GetXaxis() -> SetTitleSize(0.08);
        graph -> GetYaxis() -> SetTitleOffset(0.7);
        graph -> GetYaxis() -> SetLabelSize(0.08);
        graph -> GetYaxis() -> SetTitleSize(0.08);
        graph -> GetXaxis() -> SetNdivisions(505);
        graph -> GetYaxis() -> SetNdivisions(505);
        //graph -> SetLineColor(color);
        graph -> SetMarkerColor(color);
        graph -> SetMarkerStyle(markerstyle);
        graph -> SetMarkerSize(dot_size);
        graph -> SetLineColor(lcolor);
        graph -> SetTitle(ax.title.c_str());
    }
    void GraphErrors(TGraphErrors* graph,axrange ax){
        graph -> SetMaximum(ax.ymax);
        graph -> SetMinimum(ax.ymin);
        graph -> GetXaxis() -> SetLimits(ax.xmin,ax.xmax);
        graph -> GetXaxis() -> SetLabelSize(0.08);
        graph -> GetXaxis() -> SetTitleOffset(0.8);
        graph -> GetXaxis() -> SetTitleSize(0.08);
        graph -> GetYaxis() -> SetTitleOffset(0.8);
        graph -> GetYaxis() -> SetLabelSize(0.08);
        graph -> GetYaxis() -> SetTitleSize(0.08);
        graph -> GetXaxis() -> SetNdivisions(505);
        graph -> GetYaxis() -> SetNdivisions(505);
        graph -> SetMarkerStyle(markerstyle);
        graph -> SetMarkerSize(dot_size);
        graph -> SetMarkerColor(color);
        graph -> SetLineColor(lcolor);
        graph -> SetTitle(ax.title.c_str());
    }
    void GraphErrorsDiv(TGraphErrors* graph,axrange ax){
        graph -> SetMaximum(ax.ymax);
        graph -> SetMinimum(ax.ymin);
        graph -> GetXaxis() -> SetLimits(ax.xmin,ax.xmax);
        graph -> GetXaxis() -> SetLabelSize(0.08);
        graph -> GetXaxis() -> SetTitleOffset(0.8);
        graph -> GetXaxis() -> SetTitleSize(0.08);
        graph -> GetYaxis() -> SetTitleOffset(0.6);
        graph -> GetYaxis() -> SetLabelSize(0.08);
        graph -> GetYaxis() -> SetTitleSize(0.08);
        graph -> GetXaxis() -> SetNdivisions(505);
        graph -> GetYaxis() -> SetNdivisions(505);
        graph -> SetMarkerStyle(markerstyle);
        graph -> SetMarkerSize(dot_size);
        graph -> SetMarkerColor(color);
        graph -> SetLineColor(lcolor);
        graph -> SetTitle(ax.title.c_str());
    }
    void Hist(TH1D* hist){
        hist -> GetXaxis() -> SetLabelSize(0.07);
        hist -> GetXaxis() -> SetTitleOffset(0.8);
        hist -> GetXaxis() -> SetTitleSize(0.08);
        hist -> GetYaxis() -> SetTitleOffset(0.8);
        hist -> GetYaxis() -> SetLabelSize(0.07);
        hist -> GetYaxis() -> SetTitleSize(0.08);
        hist -> GetXaxis() -> SetNdivisions(505);
        hist -> GetYaxis() -> SetNdivisions(505);
        double xmax = hist -> GetMaximum();
        double xmin = hist -> GetMinimum();
    }
    void Graph2D(TGraph2D* graph,struct axrange ax){
        graph -> SetMaximum(ax.zmax);
        graph -> SetMinimum(ax.zmin);
        //x軸について
        graph -> GetXaxis() -> SetLimits(ax.xmin,ax.xmax);
        graph -> GetXaxis() -> SetLabelSize(0.07);
        graph -> GetXaxis() -> SetTitleOffset(1);
        graph -> GetXaxis() -> SetTitleSize(0.07);
        graph -> GetXaxis() -> SetNdivisions(505);
        //y軸について
        graph -> GetYaxis() -> SetLimits(ax.ymin,ax.ymax);
        graph -> GetYaxis() -> SetTitleOffset(1);
        graph -> GetYaxis() -> SetLabelSize(0.07);
        graph -> GetYaxis() -> SetTitleSize(0.07);
        graph -> GetYaxis() -> SetNdivisions(505);
        //z軸について
        graph -> GetZaxis() -> SetTitleOffset(1);
        graph -> GetZaxis() -> SetLabelSize(0.07);
        graph -> GetZaxis() -> SetTitleSize(0.07);
        graph -> GetZaxis() -> SetNdivisions(505);
        //Marker(必要に応じて)
        graph -> SetMarkerStyle(markerstyle);
        graph -> SetMarkerSize(dot_size);
        graph -> SetMarkerColor(color);
        graph -> SetLineColor(color);
        //タイトル付け
        graph -> SetTitle(ax.title.c_str());
    }
    double ChiValue(TGraphErrors* graph,TF1* f,int boffset){
        int npoints = graph -> GetN();
        double chi2 = 0.0;
        double x,y,ex,ey;
        rep(i,npoints){
            graph -> GetPoint(i+boffset,x,y);
            ex = graph -> GetErrorX(i+boffset);
            ey = graph -> GetErrorY(i+boffset);

            double expectedY = f -> Eval(x);
            double residual = (y-expectedY)/ey;

            chi2 += residual*residual;
        }
        return chi2;
    }
    
};