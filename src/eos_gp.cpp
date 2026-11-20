// Tabulated equation of state at mu_B = 0 (EOS_to_use 24), see eos_gp.h.
// Ported from luizafperin/MUSIC (branch main, ad67922, Luiza Perin).

#include "eos_gp.h"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

#include "util.h"

using std::string;

EOS_gp::EOS_gp(const int eos_id_in, int l_in, int sigma_in, int T_sw_in,
               const std::string &eos_type_in, const std::string &sample_in,
               const std::string &gp_type_in, const std::string &mean_in)
    : eos_id(eos_id_in), l(l_in), sigma(sigma_in), T_sw(T_sw_in),
      eos_type(eos_type_in), sample(sample_in), gp_type(gp_type_in),
      mean(mean_in) {
    set_EOS_id(eos_id);
    set_number_of_tables(0);
    set_eps_max(1e5);
    set_flag_muB(false);
    set_flag_muS(false);
    set_flag_muQ(false);
}

EOS_gp::~EOS_gp() {
    for (int itable = 0; itable < get_number_of_tables(); itable++) {
        Util::mtx_free(cs2_tb[itable], nb_length[itable], e_length[itable]);
        Util::mtx_free(entropy_tb[itable], nb_length[itable],
                       e_length[itable]);
        Util::mtx_free(energy_tb[itable], nb_length[itable],
                       e_length[itable]);
    }
    if (get_number_of_tables() > 0) {
        delete[] cs2_tb;
        delete[] entropy_tb;
        delete[] energy_tb;
    }
}

std::string EOS_gp::table_path() const {
    std::ostringstream path;
    path << get_hydro_env_path() << "/EOS/EOS-gp/Tsw_" << T_sw << "/"
         << gp_type << "/";
    if (gp_type == "unconstrained") {
        path << mean << "/";
    }
    path << "l" << l << "_s" << sigma << "/eos_" << eos_type << "_" << sample
         << "_l" << l << "_s" << sigma << ".dat";
    return path.str();
}

void EOS_gp::initialize_eos() {
    const std::string eos_file_name = table_path();
    music_message << "reading EOS_gp from " << eos_file_name;
    music_message.flush("info");

    std::ifstream eos_file(eos_file_name);
    if (!eos_file.is_open()) {
        music_message << "Could not open EOS file: " << eos_file_name;
        music_message.flush("error");
        exit(1);
    }

    set_number_of_tables(1);
    resize_table_info_arrays();
    const int ntables = get_number_of_tables();
    pressure_tb = new double **[ntables];
    temperature_tb = new double **[ntables];
    cs2_tb = new double **[ntables];
    entropy_tb = new double **[ntables];
    energy_tb = new double **[ntables];

    for (int itable = 0; itable < ntables; itable++) {
        nb_length[itable] = 1;  // mu_B = 0: one row in rho_B
        eos_file >> e_length[itable];
        string dummy;
        std::getline(eos_file, dummy);  // rest of the first line
        pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
        cs2_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                          e_length[itable]);
        entropy_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                              e_length[itable]);
        energy_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                             e_length[itable]);
        string line;
        for (int ii = 0; ii < e_length[itable]; ii++) {
            std::getline(eos_file, line);
            std::istringstream iss(line);
            double energy, temperature, pressure, cs2, entropy;
            if (!(iss >> energy >> temperature >> pressure >> cs2
                  >> entropy)) {
                music_message << "EOS_gp: cannot parse line " << ii + 2
                              << " of " << eos_file_name << ": " << line;
                music_message.flush("error");
                exit(1);
            }
            energy_tb[itable][0][ii] = energy / Util::hbarc;  // 1/fm^4
            temperature_tb[itable][0][ii] = temperature / Util::hbarc;  // 1/fm
            pressure_tb[itable][0][ii] = pressure / Util::hbarc;  // 1/fm^4
            cs2_tb[itable][0][ii] = cs2;
            entropy_tb[itable][0][ii] = entropy;  // 1/fm^3
        }
        music_message << "EOS_gp: " << e_length[itable]
                      << " rows, e up to "
                      << energy_tb[itable][0][e_length[itable] - 1]
                             * Util::hbarc
                      << " GeV/fm^3";
        music_message.flush("info");
    }
}

double EOS_gp::p_e_func(double e, double rhob, double rhoq,
                        double rhos) const {
    return get_dpOverde3(e, rhob);
}

void EOS_gp::get_pressure_with_gradients(
    double epsilon, double rhob, double rhoq, double rhos, double &p,
    double &dpde, double &dpdrhob, double &dpdrhoq, double &dpdrhos) const {
    p = get_pressure(epsilon, rhob);
    dpde = p_e_func(epsilon, rhob);
    dpdrhob = 0.0;  // mu_B = 0 EoS
    dpdrhoq = 0.0;
    dpdrhos = 0.0;
}

void EOS_gp::get_pressure_with_gradients_and_cs2(
    double epsilon, double rhob, double rhoq, double rhos, double &p,
    double &dpde, double &dpdrhob, double &dpdrhoq, double &dpdrhos,
    double &cs2) const {
    // the tabulated cs^2, as in the original implementation
    get_pressure_with_gradients(epsilon, rhob, rhoq, rhos, p, dpde, dpdrhob,
                                dpdrhoq, dpdrhos);
    cs2 = get_cs2(epsilon, rhob);
}

//! local temperature [1/fm] from e [1/fm^4]
double EOS_gp::get_temperature(double e, double rhob, double rhoq,
                               double rhos) const {
    double T = interpolate1D_nonuniform(e, 0, energy_tb, temperature_tb,
                                        7000.);
    return std::max(Util::small_eps, T);
}

//! local pressure [1/fm^4] from e [1/fm^4]
double EOS_gp::get_pressure(double e, double rhob, double rhoq,
                            double rhos) const {
    double f = interpolate1D_nonuniform(e, 0, energy_tb, pressure_tb,
                                        200. / Util::hbarc);
    return std::max(Util::small_eps, f);
}

double EOS_gp::get_cs2(double e, double rhob, double rhoq,
                       double rhos) const {
    double f = interpolate1D_nonuniform(e, 0, energy_tb, cs2_tb, 0.333);
    return std::max(Util::small_eps, f);
}

double EOS_gp::get_s2e(double s, double rhob, double rhoq,
                       double rhos) const {
    return get_s2e_finite_rhob(s, 0.0);
}

double EOS_gp::get_T2e(double T_in_GeV, double rhob, double rhoq,
                       double rhos) const {
    return get_T2e_finite_rhob(T_in_GeV, 0.0);
}
