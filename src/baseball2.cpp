///
/// Starter template for second baseball problem
///
///  Do not change the interface for running the program
///  Fill in the value of vPitch in the print statement with your solution
///  at the end of main()
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TAxis.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>

using namespace std;

struct pitch_params {
    double g;
    double phi;
    double omega;
    double b_magnus;
};

double drag_f(double v) {
    return 0.0039 + 0.0058 / (1.0 + std::exp((v - 35.0) / 5.0));
}

double deriv_x(double t, const vector<double> &y, void *params_ptr) {
    (void)t;
    (void)params_ptr;
    return y[1];
}

double deriv_vx(double t, const vector<double> &y, void *params_ptr) {
    (void)t;
    pitch_params *pars = (pitch_params *)params_ptr;
    double vx = y[1];
    double vy = y[3];
    double vz = y[5];
    double v = std::sqrt(vx * vx + vy * vy + vz * vz);
    if (v == 0.0) {
        return 0.0;
    }
    double a_drag = -drag_f(v) * v;
    double sin_phi = std::sin(pars->phi);
    double cos_phi = std::cos(pars->phi);
    double a_mx = pars->b_magnus * pars->omega * (vz * sin_phi - vy * cos_phi);
    return a_drag * vx + a_mx;
}

double deriv_y(double t, const vector<double> &y, void *params_ptr) {
    (void)t;
    (void)params_ptr;
    return y[3];
}

double deriv_vy(double t, const vector<double> &y, void *params_ptr) {
    (void)t;
    pitch_params *pars = (pitch_params *)params_ptr;
    double vx = y[1];
    double vy = y[3];
    double vz = y[5];
    double v = std::sqrt(vx * vx + vy * vy + vz * vz);
    if (v == 0.0) {
        return 0.0;
    }
    double a_drag = -drag_f(v) * v;
    double sin_phi = std::sin(pars->phi);
    double cos_phi = std::cos(pars->phi);
    double a_my = pars->b_magnus * pars->omega * vx * cos_phi;
    (void)sin_phi;
    return a_drag * vy + a_my;
}

double deriv_z(double t, const vector<double> &y, void *params_ptr) {
    (void)t;
    (void)params_ptr;
    return y[5];
}

double deriv_vz(double t, const vector<double> &y, void *params_ptr) {
    (void)t;
    pitch_params *pars = (pitch_params *)params_ptr;
    double vx = y[1];
    double vy = y[3];
    double vz = y[5];
    double v = std::sqrt(vx * vx + vy * vy + vz * vz);
    double a_drag = 0.0;
    if (v != 0.0) {
        a_drag = -drag_f(v) * v;
    }
    double sin_phi = std::sin(pars->phi);
    double cos_phi = std::cos(pars->phi);
    (void)cos_phi;
    double a_mz = -pars->b_magnus * pars->omega * vx * sin_phi;
    return -pars->g + a_drag * vz + a_mz;
}

void setup_pitch(int ip,
                 vector<double> &y0,
                 pitch_params &pars,
                 double &dt) {
    const double mph_to_mps = 0.44704;
    const double deg_to_rad = std::acos(-1.0) / 180.0;
    double v0_mph = 85.0;
    double theta_deg = 1.0;
    double phi_deg = 0.0;
    if (ip == 0) {
        v0_mph = 85.0;
        theta_deg = 1.0;
        phi_deg = 0.0;
    } else if (ip == 1) {
        v0_mph = 85.0;
        theta_deg = 1.0;
        phi_deg = 45.0;
    } else if (ip == 2) {
        v0_mph = 85.0;
        theta_deg = 1.0;
        phi_deg = 135.0;
    } else {
        v0_mph = 95.0;
        theta_deg = 1.0;
        phi_deg = 225.0;
    }
    double v0 = v0_mph * mph_to_mps;
    double theta = theta_deg * deg_to_rad;
    double omega_rpm = 1800.0;
    double omega = omega_rpm * 2.0 * std::acos(-1.0) / 60.0;

    pars.g = 9.81;
    pars.phi = phi_deg * deg_to_rad;
    pars.omega = omega;
    pars.b_magnus = 4.1e-4;

    y0[0] = 0.0;
    y0[2] = 0.0;
    y0[4] = 0.0;
    y0[1] = v0 * std::cos(theta);
    y0[3] = 0.0;
    y0[5] = v0 * std::sin(theta);

    dt = 1.0e-4;
}

void simulate_pitch(const vector<double> &y0,
                    pitch_params *pars,
                    double dt,
                    double x_plate_ft,
                    double &xend,
                    double &yend,
                    double &zend,
                    double &vxend,
                    double &vyend,
                    double &vzend,
                    vector<double> &x_vals,
                    vector<double> &y_vals,
                    vector<double> &z_vals) {
    const double ft_to_m = 0.3048;
    double x_plate_m = x_plate_ft * ft_to_m;

    vector<pfunc_t> fn_list(6);
    fn_list[0] = deriv_x;
    fn_list[1] = deriv_vx;
    fn_list[2] = deriv_y;
    fn_list[3] = deriv_vy;
    fn_list[4] = deriv_z;
    fn_list[5] = deriv_vz;

    vector<double> y = y0;
    double t = 0.0;

    double x_prev = y[0];
    double y_prev = y[2];
    double z_prev = y[4];

    x_vals.clear();
    y_vals.clear();
    z_vals.clear();

    x_vals.push_back(y[0]);
    y_vals.push_back(y[2]);
    z_vals.push_back(y[4]);

    int max_steps = 200000;
    bool crossed = false;

    for (int i = 0; i < max_steps; i += 1) {
        vector<double> y_next = RK4StepN(fn_list, y, t, dt, (void *)pars);
        double x_new = y_next[0];
        double y_new = y_next[2];
        double z_new = y_next[4];

        x_vals.push_back(x_new);
        y_vals.push_back(y_new);
        z_vals.push_back(z_new);

        if (x_prev <= x_plate_m && x_new >= x_plate_m) {
            double frac = 0.0;
            if (x_new != x_prev) {
                frac = (x_plate_m - x_prev) / (x_new - x_prev);
            }
            double x_cross = x_plate_m;
            double y_cross = y_prev + frac * (y_new - y_prev);
            double z_cross = z_prev + frac * (z_new - z_prev);
            double vx_cross = y[1] + frac * (y_next[1] - y[1]);
            double vy_cross = y[3] + frac * (y_next[3] - y[3]);
            double vz_cross = y[5] + frac * (y_next[5] - y[5]);

            xend = x_cross / ft_to_m;
            yend = y_cross / ft_to_m;
            zend = z_cross / ft_to_m;
            vxend = vx_cross;
            vyend = vy_cross;
            vzend = vz_cross;

            crossed = true;
            break;
        }

        if (z_new < -5.0) {
            break;
        }

        y = y_next;
        t += dt;
        x_prev = x_new;
        y_prev = y_new;
        z_prev = z_new;
    }

    if (!crossed) {
        double x_last = y[0];
        double y_last = y[2];
        double z_last = y[4];
        xend = x_last / ft_to_m;
        yend = y_last / ft_to_m;
        zend = z_last / ft_to_m;
        vxend = y[1];
        vyend = y[3];
        vzend = y[5];
    }
}

int main(int argc, char **argv){

  vector<double> y0(6);

  bool showPlot=true;
  int ip=1;
  int c;
  while ((c = getopt (argc, argv, "p:n")) != -1)
    switch (c) {
    case 'p':
      ip = atoi(optarg);
      break;
    case 'n':
      showPlot=false;
      break;
    }

  TString title;
  if (ip==0){
    cout << "Setting up initial conditions for slider" << endl;
  }
  else if (ip==1){
    cout << "Setting up initial conditions for curveball" << endl;
  }
  else if (ip==2){
    cout << "Setting up initial conditions for screwball" << endl;
  }
  else {
    cout << "Setting up initial conditions for fastball" << endl;
  }

  TApplication theApp("App", &argc, argv);
  gROOT->SetBatch(kTRUE);

  double xend=60.6;
  double yend=0;
  double zend=0;
  double vxend=0;
  double vyend=0;
  double vzend=0;

  const int npitches = 4;
  double x_res[npitches];
  double y_res[npitches];
  double z_res[npitches];
  double vx_res[npitches];
  double vy_res[npitches];
  double vz_res[npitches];

  pitch_params pars[npitches];
  double dt[npitches];
  vector<double> x_vals[npitches];
  vector<double> y_vals[npitches];
  vector<double> z_vals[npitches];

  for (int k = 0; k < npitches; k += 1) {
      setup_pitch(k, y0, pars[k], dt[k]);
      simulate_pitch(y0, &pars[k], dt[k], 60.6,
                     x_res[k], y_res[k], z_res[k],
                     vx_res[k], vy_res[k], vz_res[k],
                     x_vals[k], y_vals[k], z_vals[k]);
  }

  int ip_sel = ip;
  if (ip_sel < 0) {
      ip_sel = 0;
  }
  if (ip_sel > 3) {
      ip_sel = 3;
  }

  xend = x_res[ip_sel];
  yend = y_res[ip_sel];
  zend = z_res[ip_sel];
  vxend = vx_res[ip_sel];
  vyend = vy_res[ip_sel];
  vzend = vz_res[ip_sel];

  printf("********************************\n");
  printf("Coordinates when x=60 feet\n");
  printf("(x,y,x) = (%lf,%lf,%lf)\n",xend,yend,zend);
  printf("(vx,vy,vz) = (%lf,%lf,%lf)\n",vxend,vyend,vzend);
  printf("********************************\n");

  if (showPlot){
    TCanvas *canvas = new TCanvas("c_pitches",
                                  "Baseball pitch trajectories",
                                  900, 900);
    canvas->Divide(2, 2);

    const double ft_to_m = 0.3048;
    const char *ptitle[npitches] = {
        "Slider", "Curveball", "Screwball", "Fastball"
    };

    for (int k = 0; k < npitches; k += 1) {
        canvas->cd(k + 1);
        int n = static_cast<int>(x_vals[k].size());
        TGraph *g_z = new TGraph(n);
        TGraph *g_y = new TGraph(n);
        for (int i = 0; i < n; i += 1) {
            double x_ft = x_vals[k][i] / ft_to_m;
            double y_ft = y_vals[k][i] / ft_to_m;
            double z_ft = z_vals[k][i] / ft_to_m;
            g_z->SetPoint(i, x_ft, z_ft);
            g_y->SetPoint(i, x_ft, y_ft);
        }
        TString full_title = TString(ptitle[k]) + ";x(ft);y(ft)/z(ft)";
        g_z->SetTitle(full_title);
        g_z->SetLineWidth(2);
        g_y->SetLineStyle(2);
        g_y->SetLineWidth(2);
        g_z->Draw("AL");
        g_z->GetXaxis()->SetLimits(0.0, 60.6);
        g_z->SetMinimum(-4.0);
        g_z->SetMaximum(2.0);
        g_y->Draw("L SAME");
    }

    canvas->Print("pitches.pdf");
  }
  
  return 0;
}
