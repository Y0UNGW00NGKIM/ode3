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
#include <unistd.h>
#include <cmath>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

typedef struct projectile_params {
    double g;
    double mass;
    double drag_k;
} projectile_params_t;

double deriv_z(double t, const vector<double> &y, void *params_ptr) {
    (void)t;
    (void)params_ptr;
    return y[1];
}

double deriv_v(double t, const vector<double> &y, void *params_ptr) {
    (void)t;
    projectile_params_t *params = (projectile_params_t *)params_ptr;
    double v = y[1];

    double drag_acc = 0.0;
    if (params->drag_k != 0.0) {
        drag_acc = -(params->drag_k / params->mass) * std::fabs(v) * v;
    }

    return -params->g + drag_acc;
}

double total_energy(const projectile_params_t &params, double z, double v) {
    double kinetic = 0.5 * params.mass * v * v;
    double potential = params.mass * params.g * z;
    return kinetic + potential;
}

void run_energy_test(unsigned int canvas_width,
                     unsigned int canvas_height,
                     int nsteps,
                     double tmax) {
    projectile_params_t params;
    params.g = 9.81;
    params.mass = 1.0;
    params.drag_k = 0.0;

    vector<pfunc_t> fn_list(2);
    fn_list[0] = deriv_z;
    fn_list[1] = deriv_v;

    vector<double> y(2);
    y[0] = 0.0;
    y[1] = 20.0;

    double t0 = 0.0;
    double t = t0;

    vector<TGraph> graphs = RK4SolveN(fn_list, y, nsteps, t, tmax, &params, 0);

    int n_points = graphs[0].GetN();
    TGraph energy_vs_t;

    double energy_min = 0.0;
    double energy_max = 0.0;
    double energy_sum = 0.0;

    for (int i = 0; i < n_points; i += 1) {
        double t_val = 0.0;
        double z_val = 0.0;
        double v_val = 0.0;

        graphs[0].GetPoint(i, t_val, z_val);
        graphs[1].GetPoint(i, t_val, v_val);

        double energy = total_energy(params, z_val, v_val);
        energy_vs_t.SetPoint(i, t_val, energy);

        if (i == 0) {
            energy_min = energy;
            energy_max = energy;
        } else {
            if (energy < energy_min) {
                energy_min = energy;
            }
            if (energy > energy_max) {
                energy_max = energy;
            }
        }

        energy_sum += energy;
    }

    double energy_avg = energy_sum / static_cast<double>(n_points);
    double relative_variation = (energy_max - energy_min) / energy_avg;

    cout << "Energy conservation test (no air resistance)" << endl;
    cout << "  nsteps             = " << nsteps << endl;
    cout << "  t_max              = " << tmax << " s" << endl;
    cout << "  E_min              = " << energy_min << " J" << endl;
    cout << "  E_max              = " << energy_max << " J" << endl;
    cout << "  relative variation = " << relative_variation << endl;

    graphs[0].SetTitle("Height vs time;time t [s];z(t) [m]");
    energy_vs_t.SetTitle("Total energy vs time;time t [s];E(t) [J]");

    TCanvas *canvas = new TCanvas("c_vterm_energy",
                                  "Energy conservation",
                                  canvas_width,
                                  canvas_height);
    canvas->Divide(1, 2);

    canvas->cd(1);
    graphs[0].Draw("AL");

    canvas->cd(2);
    energy_vs_t.Draw("AL");

    canvas->Draw();
}

void run_terminal_velocity_scan(unsigned int canvas_width,
                                unsigned int canvas_height) {
    projectile_params_t params;
    params.g = 9.81;
    params.drag_k = 0.1;

    vector<pfunc_t> fn_list(2);
    fn_list[0] = deriv_z;
    fn_list[1] = deriv_v;

    double mass_min = 0.001;
    double mass_max = 10.0;
    int mass_steps = 40;

    double t0 = 0.0;
    double tmax = 20.0;
    int nsteps = 2000;

    TGraph vt_vs_mass_num;
    TGraph vt_vs_mass_ana;
    TGraph rel_err_vs_mass;

    for (int i = 0; i < mass_steps; i += 1) {
        double mass = mass_min
                      + (mass_max - mass_min)
                        * static_cast<double>(i)
                        / static_cast<double>(mass_steps - 1);

        params.mass = mass;

        vector<double> y(2);
        y[0] = 0.0;
        y[1] = 0.0;

        double t = t0;
        vector<TGraph> graphs = RK4SolveN(fn_list, y, nsteps, t, tmax, &params, 0);
        (void)graphs;

        double v_final = y[1];
        double vt_numeric = std::fabs(v_final);
        double vt_analytic = std::sqrt(params.mass * params.g / params.drag_k);
        double rel_err = 0.0;
        if (vt_analytic != 0.0) {
            rel_err = (vt_numeric - vt_analytic) / vt_analytic;
        }

        cout << "m = " << mass
             << " kg, v_t (numeric)  = " << vt_numeric
             << " m/s, v_t (analytic) = " << vt_analytic
             << " m/s" << endl;

        vt_vs_mass_num.SetPoint(i, mass, vt_numeric);
        vt_vs_mass_ana.SetPoint(i, mass, vt_analytic);
        rel_err_vs_mass.SetPoint(i, mass, rel_err);
    }

    vt_vs_mass_num.SetTitle("Terminal speed vs mass;mass m [kg];v_{t}(m) [m/s]");
    rel_err_vs_mass.SetTitle("Relative error vs mass;mass m [kg];(v_{num}-v_{ana})/v_{ana}");

    TCanvas *canvas = new TCanvas("c_vterm_vt",
                                  "Terminal velocity vs mass",
                                  canvas_width,
                                  canvas_height);
    canvas->Divide(1, 2);

    canvas->cd(1);
    vt_vs_mass_num.Draw("AL*");
    vt_vs_mass_ana.Draw("L");

    TLegend *legend = new TLegend(0.15, 0.7, 0.45, 0.88);
    legend->AddEntry(&vt_vs_mass_num, "numeric v_{t}", "lp");
    legend->AddEntry(&vt_vs_mass_ana, "analytic v_{t}", "l");
    legend->Draw();

    canvas->cd(2);
    rel_err_vs_mass.Draw("AL*");

    canvas->Print("vterm.pdf");
}


int main(int argc, char **argv) {
    bool do_energy_test = false;
    int energy_nsteps = 500;
    double energy_tmax = 5.0;

    int opt = 0;
    while ((opt = getopt(argc, argv, "En:t:")) != -1) {
        if (opt == 'E') {
            do_energy_test = true;
        } else if (opt == 'n') {
            energy_nsteps = std::atoi(optarg);
        } else if (opt == 't') {
            energy_tmax = std::atof(optarg);
        } else {
            std::fprintf(stderr,
                         "Usage: %s [-E] [-n nsteps] [-t tmax]\n",
                         argv[0]);
            return 1;
        }
    }

    TApplication app("vterm", &argc, argv);
    gROOT->SetBatch(kTRUE);

    unsigned int canvas_width = 800;
    unsigned int canvas_height = 600;

    if (do_energy_test) {
        run_energy_test(canvas_width, canvas_height,
                        energy_nsteps, energy_tmax);
    } else {
        run_terminal_velocity_scan(canvas_width, canvas_height);
    }

    return 0;
}

