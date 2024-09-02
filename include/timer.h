#pragma once

#include  <chrono>
#include <iomanip>
#include <sstream>
#include  <string>

typedef std::chrono::time_point<std::chrono::high_resolution_clock> Timepoint;
typedef std::chrono::milliseconds Millis; typedef std::chrono::hours Hours;
typedef std::chrono::seconds Seconds; typedef std::chrono::minutes Minutes;

namespace Timer {
    long Elapsed(const Timepoint& start);
    std::string Format(long millis);
    std::string Local();
    Timepoint Now();
};
