#include "timer.h"

long Timer::Elapsed(Timepoint start) {
    return std::chrono::duration_cast<Millis>(Now() - start).count();
};

std::string Timer::Format(long ms) {
    long hours = ms / 3600000, mins = ms % 3600000 / 60000, secs = ms % 60000 / 1000; ms = ms % 1000;
    std::stringstream strstream; strstream << std::setfill('0'); strstream << std::setw(2);
    strstream << hours <<  ":" << std::setw(2) << mins << ":" << std::setw(2);
    strstream << secs << "." << std::setw(3) << ms; return strstream.str();
}

std::string Timer::Local() {
    long t = std::time(0); auto tm = *std::localtime(&t); std::stringstream strstream;
    strstream << std::put_time(&tm, "%a %b %d %H:%M:%S %Y"); return strstream.str();
}

Timer::Timepoint Timer::Now() {
    return std::chrono::high_resolution_clock().now();
};
