///
/// Starter template for first baseball problem
/// Solve for the initial speed of the pitch given the initial parameters
/// xend : distance to home plate [18.5] m
/// z0 : height of release of ball [1.4] m
/// theta0 : angle of release above horizontal [1] degree
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
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

struct Params {
    double g;
    double m;
    double d;
    double b;
    double c;
};

double deriv_x(double t, const vector<double> &y, void *params_ptr) {
    (void)t;
    (void)params_ptr;
    return y[2];
}

double deriv_z(double t, const vector<double> &y, void *params_ptr) {
    (void)t;
    (void)params_ptr;
    return y[3];
}

double deriv_vx(double t, const vector<double> &y, void *params_ptr) {
    (void)t;
    Params *pars = (Params *)params_ptr;
    double vx = y[2];
    double vz = y[3];
    double v = std::sqrt(vx * vx + vz * vz);
    if (v == 0.0) {
        return 0.0;
    }
    double b_eff = pars->b * pars->d;
    double c_eff = pars->c * pars->d * pars->d;
    double force_mag = b_eff * v + c_eff * v * v;
    double ax = -force_mag * vx / (pars->m * v);
    return ax;
}

double deriv_vz(double t, const vector<double> &y, void *params_ptr) {
    (void)t;
    Params *pars = (Params *)params_ptr;
    double vx = y[2];
    double vz = y[3];
    double v = std::sqrt(vx * vx + vz * vz);
    if (v == 0.0) {
        return -pars->g;
    }
    double b_eff = pars->b * pars->d;
    double c_eff = pars->c * pars->d * pars->d;
    double force_mag = b_eff * v + c_eff * v * v;
    double az_drag = -force_mag * vz / (pars->m * v);
    double az = az_drag - pars->g;
    return az;
}

double height_at_plate(double v0,
                       Params *pars,
                       double xend,
                       double z0,
                       double theta0_deg) {
    vector<pfunc_t> fn_list(4);
    fn_list[0] = deriv_x;
    fn_list[1] = deriv_z;
    fn_list[2] = deriv_vx;
    fn_list[3] = deriv_vz;

    vector<double> y(4);
    double theta = theta0_deg * std::acos(-1.0) / 180.0;
    y[0] = 0.0;
    y[1] = z0;
    y[2] = v0 * std::cos(theta);
    y[3] = v0 * std::sin(theta);

    double t = 0.0;
    double dt = 0.001;
    int max_steps = 50000;

    double x_prev = y[0];
    double z_prev = y[1];

    for (int i = 0; i < max_steps; i += 1) {
        vector<double> y_next = RK4StepN(fn_list, y, t, dt, (void *)pars);
        double x_new = y_next[0];
        double z_new = y_next[1];

        if (x_prev <= xend && x_new >= xend) {
            double frac = 0.0;
            if (x_new != x_prev) {
                frac = (xend - x_prev) / (x_new - x_prev);
            }
            double z_cross = z_prev + frac * (z_new - z_prev);
            return z_cross;
        }

        if (z_new < 0.0 && x_new < xend) {
            break;
        }

        y = y_next;
        t += dt;
        x_prev = x_new;
        z_prev = z_new;
    }

    return -1.0;
}

double find_v_pitch(Params *pars,
                    double xend,
                    double z0,
                    double theta0_deg,
                    double zend) {
    double v_low = 20.0;
    double v_high = 60.0;

    double f_low = height_at_plate(v_low, pars, xend, z0, theta0_deg) - zend;
    double f_high = height_at_plate(v_high, pars, xend, z0, theta0_deg) - zend;

    for (int iter = 0; iter < 40; iter += 1) {
        double v_mid = 0.5 * (v_low + v_high);
        double f_mid = height_at_plate(v_mid, pars, xend, z0, theta0_deg) - zend;

        if (f_mid == 0.0 || 0.5 * (v_high - v_low) < 1e-4) {
            return v_mid;
        }

        if (f_low * f_mid <= 0.0) {
            v_high = v_mid;
            f_high = f_mid;
        } else {
            v_low = v_mid;
            f_low = f_mid;
        }
    }

    return 0.5 * (v_low + v_high);
}

int main(int argc, char **argv){

  // examples of parameters
  Params pars;
  pars.g=9.81;
  pars.m=0.145;    
  pars.d=0.075;   
  pars.b=1.6e-4;  
  pars.c=0.25;

  double xend=18.5;       // meters to plate
  double z0=1.4;             // height of release [m]
  double theta0=1;         // angle of velocity at release (degrees)
                                      // convert to radians before using!
  bool showPlot=false;    // keep this flag false by default
  
  // allow changing the parameters from the command line
  int c;
  while ((c = getopt (argc, argv, "x:z:t:p")) != -1)
    switch (c) {
    case 'x':
      xend = atof(optarg);
      break;
    case 'z':
      z0 = atof(optarg);
      break;
    case 't':
      theta0 = atof(optarg);
      break;
    case 'p':
      showPlot=true;
      break;
    case '?':
      fprintf (stderr, "Unknown option `%c'.\n", optopt);
    }
  TApplication theApp("App", &argc, argv); // init ROOT App for displays


  double zend = 0.9;
  double vPitch = find_v_pitch(&pars, xend, z0, theta0, zend);

  // do not change these lines
  printf("********************************\n");
  printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n",xend,z0,theta0);
  printf("v_pitch = %lf m/s\n",vPitch);
  printf("********************************\n");

  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}
