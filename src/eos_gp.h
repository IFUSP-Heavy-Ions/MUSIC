// Tabulated equation of state at mu_B = 0 (EOS_to_use 24): one table per
// sample of a Gaussian-process EoS ensemble, read from
//   $HYDROPROGRAMPATH/EOS/EOS-gp/Tsw_<T_sw>/<gp_type>/[<EOS_gp_mean>/]
//       l<EOS_gp_l>_s<EOS_gp_sigma>/eos_<EOS_gp_type>_<EOS_gp_sample>_l<l>_s<sigma>.dat
// File: first line the number of rows N, then N rows of
//   e [GeV/fm^3]  T [GeV]  P [GeV/fm^3]  cs^2  s [1/fm^3]
// on a non-uniform grid in e, increasing.  Linear interpolation inside the
// table; above it an exponential approach to fixed asymptotic values.
// Ported from luizafperin/MUSIC (branch main, ad67922, Luiza Perin).

#ifndef SRC_EOS_GP_H_
#define SRC_EOS_GP_H_

#include <string>

#include "eos_base.h"

class EOS_gp : public EOS_base {
 private:
    const int eos_id;
    int l;
    int sigma;
    int T_sw;
    std::string eos_type;
    std::string sample;
    std::string gp_type;
    std::string mean;

    // tables of this EoS only (the base class frees pressure_tb and
    // temperature_tb; these three are freed here)
    double ***cs2_tb = nullptr;
    double ***entropy_tb = nullptr;
    double ***energy_tb = nullptr;

 public:
    EOS_gp(const int eos_id_in, int l_in, int sigma_in, int T_sw_in,
           const std::string &eos_type_in, const std::string &sample_in,
           const std::string &gp_type_in, const std::string &mean_in);
    ~EOS_gp();

    std::string table_path() const;

    void initialize_eos();
    double p_e_func(double e, double rhob, double rhoq = 0.0,
                    double rhos = 0.0) const;
    double get_temperature(double e, double rhob, double rhoq = 0.0,
                           double rhos = 0.0) const;
    double get_pressure(double e, double rhob, double rhoq = 0.0,
                        double rhos = 0.0) const;
    double get_cs2(double e, double rhob, double rhoq = 0.0,
                   double rhos = 0.0) const;
    double get_s2e(double s, double rhob, double rhoq = 0.0,
                   double rhos = 0.0) const;
    double get_T2e(double T, double rhob, double rhoq = 0.0,
                   double rhos = 0.0) const;

    void get_pressure_with_gradients(
        double epsilon, double rhob, double rhoq, double rhos, double &p,
        double &dpde, double &dpdrhob, double &dpdrhoq,
        double &dpdrhos) const;
    void get_pressure_with_gradients_and_cs2(
        double epsilon, double rhob, double rhoq, double rhos, double &p,
        double &dpde, double &dpdrhob, double &dpdrhoq, double &dpdrhos,
        double &cs2) const;

    void check_eos() const { check_eos_no_muB(); }
};

#endif  // SRC_EOS_GP_H_
