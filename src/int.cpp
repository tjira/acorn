#include <libint2.hpp>

extern "C" void initialize(void) {
    libint2::initialize();
}

extern "C" void finalize(void) {
    libint2::finalize();
}
