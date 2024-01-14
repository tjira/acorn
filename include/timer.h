#pragma once

#include <libint2/diis.h>
#include <bits/stdc++.h>

namespace Timer {
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> Timepoint;
    typedef std::chrono::milliseconds Millis;

    // methods
    long Elapsed(Timepoint start); std::string Local();
    Timepoint Now(); std::string Format(long millis);
};
