#include <TROOT.h>
#include <iostream>
using namespace std;
//
struct axrange{
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
};

class Setting{
public:
    double dot_size;
    double markerstyle;
    Color_t color;
public:
    void Graph(TGraph* graph,double xmin,double xmax,double ymin,double ymax){
        graph -> SetMaximum(ymax);
        graph -> SetMinimum(ymin);
        graph -> GetXaxis() -> SetLimits(xmin,xmax);
        graph -> GetXaxis() -> SetLabelSize(0.08);
        graph -> GetXaxis() -> SetTitleOffset(0.9);
        graph -> GetXaxis() -> SetTitleSize(0.08);
        graph -> GetYaxis() -> SetTitleOffset(0.7);
        graph -> GetYaxis() -> SetLabelSize(0.08);
        graph -> GetYaxis() -> SetTitleSize(0.08);
        graph -> GetXaxis() -> SetNdivisions(505);
        graph -> GetYaxis() -> SetNdivisions(505);
        graph -> SetLineColor(color);
        graph -> SetMarkerColor(color);
        graph -> SetMarkerStyle(markerstyle);
        graph -> SetMarkerSize(dot_size);
    }
    void GraphErrors(TGraphErrors* graph,double xmin,double xmax,double ymin,double ymax){
        graph -> SetMaximum(ymax);
        graph -> SetMinimum(ymin);
        graph -> GetXaxis() -> SetLimits(xmin,xmax);
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
        graph -> SetLineColor(color);
    }
    void GraphErrorsDiv(TGraphErrors* graph,double xmin,double xmax,double ymin,double ymax){
        graph -> SetMaximum(ymax);
        graph -> SetMinimum(ymin);
        graph -> GetXaxis() -> SetLimits(xmin,xmax);
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
        graph -> SetLineColor(color);
    }
    void Hist(TH1D* hist){
        hist -> GetXaxis() -> SetLabelSize(0.08);
        hist -> GetXaxis() -> SetTitleOffset(0.8);
        hist -> GetXaxis() -> SetTitleSize(0.08);
        hist -> GetYaxis() -> SetTitleOffset(0.7);
        hist -> GetYaxis() -> SetLabelSize(0.08);
        hist -> GetYaxis() -> SetTitleSize(0.08);
        hist -> GetXaxis() -> SetNdivisions(505);
        hist -> GetYaxis() -> SetNdivisions(505);
    }
    void Graph2D(TGraph2D* graph,struct axrange ax){
        graph -> SetMaximum(ax.zmax);
        graph -> SetMinimum(ax.zmin);
        //x軸について
        graph -> GetXaxis() -> SetLimits(ax.xmin,ax.xmax);
        graph -> GetXaxis() -> SetLabelSize(0.08);
        graph -> GetXaxis() -> SetTitleOffset(0.8);
        graph -> GetXaxis() -> SetTitleSize(0.08);
        graph -> GetXaxis() -> SetNdivisions(505);
        //y軸について
        graph -> GetYaxis() -> SetLimits(ax.ymin,ax.ymax);
        graph -> GetYaxis() -> SetTitleOffset(0.6);
        graph -> GetYaxis() -> SetLabelSize(0.08);
        graph -> GetYaxis() -> SetTitleSize(0.08);
        graph -> GetYaxis() -> SetNdivisions(505);
        //z軸について
        graph -> GetZaxis() -> SetTitleOffset(0.6);
        graph -> GetZaxis() -> SetLabelSize(0.08);
        graph -> GetZaxis() -> SetTitleSize(0.08);
        graph -> GetZaxis() -> SetNdivisions(505);
    }
    
};
