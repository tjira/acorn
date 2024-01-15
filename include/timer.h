#pragma once

#include <libint2/diis.h>
#include <bits/stdc++.h>

namespace Timer {
    // typedefs of the timepoint and time intervals
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> Timepoint;
    typedef std::chrono::milliseconds Millis; typedef std::chrono::hours Hours;
    typedef std::chrono::seconds Seconds; typedef std::chrono::minutes Minutes;

    // methods
    long Elapsed(Timepoint start); std::string Local();
    Timepoint Now(); std::string Format(long millis);
};
