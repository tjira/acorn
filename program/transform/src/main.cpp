#include "transform.h"

int main(int argc, char** argv) {
    auto [opt, timers] = Acorn::Transform::initialize(argc, argv); Acorn::Transform::run(opt, timers);
}
