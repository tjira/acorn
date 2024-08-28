#pragma once

#include "linalg.h"

namespace Acorn {
    namespace CDYN {
        class LandauZener {
        public:
            LandauZener(int nstate, int points, bool adiabatic); LandauZener();

            // function to perform the Landau-Zener jump
            std::vector<std::tuple<int, double, bool>> jump(const Matrix& U, int state, int i, double tstep);

            // matrix getters
            const Matrix& getEd() const {return ed;} const Matrix& getDed() const {return ded;} const Matrix& getDded() const {return dded;}

        private:
            Matrix ed, ded, dded; std::vector<std::vector<int>> combs; bool adiabatic;
        };
    }
}
