#pragma once

#include "tensor.h"

namespace Acorn {
    namespace CDYN {
        class LandauZener {
        public:
            LandauZener(int nstate, int points, bool adiabatic); LandauZener();

            // function to perform the Landau-Zener jump
            std::vector<std::tuple<int, double, bool>> jump(const Eigen::MatrixXd& U, int state, int i, double tstep);

            // matrix getters
            const Eigen::MatrixXd& getEd() const {return ed;} const Eigen::MatrixXd& getDed() const {return ded;} const Eigen::MatrixXd& getDded() const {return dded;}

        private:
            Eigen::MatrixXd ed, ded, dded; std::vector<std::vector<int>> combs; bool adiabatic;
        };
    }
}
