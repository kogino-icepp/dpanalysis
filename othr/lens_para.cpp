#include <iostream>
#include <cmath>
#include <vector>
#include "../headers/fitter.h"
using namespace std;

const double a = 4.307257;
const double Rh = 35.6176;
const double window = 76.819+171.47+57.32;
const double f_place[4]={0,57.32,window,window+1500};
const double focal_length = 56.87;
const double freqs[3] = {220,240,260};
double gauss_beam_width(double waist_position, double waist_size, double freq, double z) {
    double z_0 = waist_position;
    double w0 = waist_size;
    double c = 3e11; // [mm/s]
    double freq_Hz = freq * 1e9;
    double lamda = c / freq_Hz;
    double zc = M_PI * w0 * w0 / lamda;
    double w = w0 * sqrt(1 + pow((z - z_0), 2) / (zc * zc));
    return w;
}
double ztrue(double w,double freq){
    double c = 3e11; // [mm/s]
    double w0 = 2.51;
    double lambda = c/(freq*pow(10,9));
    double z0 = M_PI*w0*w0/lambda;
    return z0*sqrt(pow(w/w0,2)-1);
}
std::pair<double, double> lens(double din, double w0_in, double freq) {
    double f = focal_length;
    double c = 3e11; // [mm/s]
    double freq_Hz = freq * 1e9;
    double lamda = c / freq_Hz;
    double zc = M_PI * w0_in * w0_in / lamda;
    double A = din / f - 1;
    double dout = f * (1 + A / (A * A + zc * zc / (f * f)));
    double w0_out = w0_in / sqrt(A * A + zc * zc / (f * f));
    return {dout, w0_out};
}
double initial_w0(double x, double y, double freq) {
    double c = 3e11; // [mm/s]
    double freq_Hz = freq * 1e9;
    double wave_len = c / freq_Hz;
    double a = x;
    double Rh = y;
    double w = a * 0.644;
    return w / (1 + pow(M_PI * w * w / (wave_len * Rh), 2));
}
double initial_z0(double x, double y, double freq) {
    double c = 3e11; // [mm/s]
    double freq_Hz = freq * 1e9;
    double wave_len = c / freq_Hz;
    double a = x;
    double Rh = y;
    double w = a * 0.644;
    return Rh / (1 + pow(wave_len * Rh / (M_PI * w * w), 2));
}
void lens_para(){
    TCanvas *c1 = new TCanvas("c1","My Canvas",10,10,700,500);
    c1 -> SetMargin(0.14,0.11,0.2,0.1);
    Setting st;
    st.dot_size = 0.8;
    st.markerstyle = 20;
    st.color = kBlue;
    st.lcolor = kBlue;
    TGraph* graph1 = new TGraph;
    TGraph* graph2 = new TGraph;
    TGraph* graph3 = new TGraph;
    TGraph* graph = new TGraph;
    

    prep(freq,216,265){
        double beam_params[3][2];
        beam_params[0][0] = -initial_z0(a,Rh,freq);
        beam_params[0][1] = initial_w0(a,Rh,freq);
        //cout << initial_w0(a,Rh,freq) << endl;
        rep(i,2){
            double din = f_place[i+1]-(beam_params[i][0]+f_place[i]);
            beam_params[i+1][0] = lens(din,beam_params[i][1],freq).first;
            beam_params[i+1][1] = lens(din,beam_params[i][1],freq).second;
        }
        graph -> SetPoint(freq-216,freq,beam_params[2][0]);
        //cout << "freq: " << freq << " => " << beam_params[2][0] << endl;
    }
    axrange ax = {215,265,50,60,0,1,";Freq[GHz];waist_position[mm]"};
    st.Graph(graph,ax);
    graph -> Draw("Ap");
}