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
    return w / sqrt(1 + pow(M_PI * w * w / (wave_len * Rh), 2));
}
double initial_z0(double x, double y, double freq) {
    double c = 3e11; // [mm/s]
    double freq_Hz = freq * 1e9;
    double wave_len = c / freq_Hz;
    double a = x;
    double Rh = y;
    double w = a * 0.644;
    return Rh / sqrt(1 + pow(wave_len * Rh / (M_PI * w * w), 2));
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
    cout << "window: " << window << endl;
    for(int i=0;i<3;i++){
        vector<vector<double>> beam_params;
        vector<double> z_mesh,w_values;

        double initial_z = -initial_z0(a,Rh,freqs[i]);
        double initial_w = initial_w0(a,Rh,freqs[i]);
        beam_params.push_back({initial_z,initial_w});

        for(int j=0;j<2;j++){
            double din = f_place[j + 1] - (beam_params[j][0] + f_place[j]);
            auto lens_result = lens(din, beam_params[j][1], freqs[i]);
            beam_params.push_back({lens_result.first, lens_result.second});
        }
        int bin = 0;
        for (int j = 0; j < 3; ++j) {
            double z_start = (j == 0) ? beam_params[j][0] + f_place[j] : f_place[j];
            double z_end = f_place[j + 1];
            int num_points = 1000;
            double z_step = (z_end - z_start) / num_points;

            for (int k = 0; k < num_points; ++k) {
                double z_current = z_start + k * z_step;
                double w = gauss_beam_width(beam_params[j][0] + f_place[j], beam_params[j][1], freqs[i], z_current);
                z_mesh.push_back(z_current);
                w_values.push_back(w);
                if(i==0){
                    graph1 -> SetPoint(bin,z_current,1.4*w);
                    bin++;
                }
                else if(i==1){
                    graph2 -> SetPoint(bin,z_current,1.4*w);
                    bin++;
                }
                else{
                    graph3 -> SetPoint(bin,z_current,1.4*w);
                    bin++;
                }
                
            }
        }
    }
    axrange ax = {-50,1800,0,250,0,1,";distance from horn aperture[mm];beam radious[mm]"};
    st.Graph(graph1,ax);
    graph1 -> SetLineColor(kBlue);
    graph2 -> SetLineColor(kRed);
    graph3 -> SetLineColor(kGreen);
    graph1 -> Draw("AL");
    graph2 -> Draw("L");
    graph3 -> Draw("L");
    TLegend* legend = new TLegend(0.25,0.5,0.45,0.8);
    legend -> AddEntry(graph1,"220GHz","l");
    legend -> AddEntry(graph2,"240GHz","l");
    legend -> AddEntry(graph3,"260GHz","l");
    //legend -> Draw();
    cout << "ztrue: " << ztrue(180,240) << endl;

}