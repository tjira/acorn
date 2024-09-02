#include "timer.h"

long Timer::Elapsed(const Timepoint& start) {
    return std::chrono::duration_cast<Millis>(Now() - start).count();
};

std::string Timer::Format(long ms) {
    // convert milliseconds to hours, minutes, seconds, and milliseconds
    long hours = ms / 3600000, mins = ms % 3600000 / 60000, secs = ms % 60000 / 1000; ms = ms % 1000;

    // define a string stream with zero padding
    std::stringstream strstream; strstream << std::setfill('0'); strstream << std::setw(2);

    // format the time as "hours:minutes:seconds.milliseconds"
    strstream << hours <<  ":" << std::setw(2) << mins << ":" << std::setw(2) << secs << "." << std::setw(3) << ms; return strstream.str();
}

std::string Timer::Local() {
    // get the current local time
    std::time_t t = std::time(0); auto tm = *std::localtime(&t); std::stringstream strstream;

    // format the time as "day month data hour:minute:second year"
    strstream << std::put_time(&tm, "%a %b %d %H:%M:%S %Y"); return strstream.str();
}

Timepoint Timer::Now() {
    return std::chrono::high_resolution_clock().now();
};
