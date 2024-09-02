#pragma once

#include    "system.h"
#include     "timer.h"
#include "transform.h"

class CoupledCluster {
public:
    CoupledCluster(const Input::HartreeFock::CoupledCluster& input) : input(input) {}

    static torch::Tensor ccd_amplitudes(const System& system, const torch::Tensor& J_MS_AP, const torch::Tensor& E_MS_D, const torch::Tensor& T2);
    static double ccd_energy(const System& system, const torch::Tensor& J_MS_AP, const torch::Tensor& T2);
    static double ccsd_energy(const System& system, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP, const torch::Tensor& T1, const torch::Tensor& T2);
    static std::tuple<torch::Tensor, torch::Tensor> ccsd_amplitudes(const System& system, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP, const torch::Tensor& E_MS_S, const torch::Tensor& E_MS_D, const torch::Tensor& T1, const torch::Tensor& T2);
    static double perturbation_triple(const System& system, const torch::Tensor& J_MS_AP, const torch::Tensor& E_MS_T, const torch::Tensor& T1, const torch::Tensor& T2);

    std::tuple<torch::Tensor, double> ccd_scf(const System& system, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP) const;
    std::tuple<torch::Tensor, torch::Tensor, double> ccsd_scf(const System& system, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP) const;
    std::string get_name() const;
    double run(const System& system, const torch::Tensor& F_MS, const torch::Tensor& J_MS_AP) const;

private:
    Input::HartreeFock::CoupledCluster input;
};

