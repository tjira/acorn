#include "qdyn.h"

int main(int argc, char** argv) {
    auto[opt, timers] = Acorn::QDYN::initialize(argc, argv); Acorn::QDYN::run(opt, timers);
}
