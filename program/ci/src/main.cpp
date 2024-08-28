#include "ci.h"

int main(int argc, char** argv) {
    auto[opt, timers] = Acorn::CI::initialize(argc, argv); Acorn::CI::run(opt, timers);
}
