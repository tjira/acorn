const std = @import("std");

const cdn = @import("acorn").cdn;
const vec = @import("acorn").vec;

const Vector = @import("acorn").Vector;

test "cdyn_fssh_doubleState1D_1" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.72; (try pop.ptr(1)).* = 0.28;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "doubleState1D_1"; opt.fewest_switches = .{}; opt.initial_conditions.state = 1;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_doubleState1D_1" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.67; (try pop.ptr(1)).* = 0.33;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "doubleState1D_1"; opt.landau_zener = .{}; opt.initial_conditions.state = 1;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_doubleState1D_2" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.96; (try pop.ptr(1)).* = 0.04;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "doubleState1D_2"; opt.fewest_switches = .{}; opt.initial_conditions.state = 1;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_doubleState1D_2" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.92; (try pop.ptr(1)).* = 0.08;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "doubleState1D_2"; opt.landau_zener = .{}; opt.initial_conditions.state = 1;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_tripleState1D_1" {
    var pop = try Vector(f64).init(3, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.1; (try pop.ptr(1)).* = 0.11; (try pop.ptr(2)).* = 0.79;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "tripleState1D_1"; opt.fewest_switches = .{}; opt.initial_conditions.state = 2;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_tripleState1D_1" {
    var pop = try Vector(f64).init(3, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.08; (try pop.ptr(1)).* = 0.03; (try pop.ptr(2)).* = 0.89;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "tripleState1D_1"; opt.landau_zener = .{}; opt.initial_conditions.state = 2;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_tripleState1D_2" {
    var pop = try Vector(f64).init(3, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.74; (try pop.ptr(1)).* = 0.06; (try pop.ptr(2)).* = 0.2;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "tripleState1D_2"; opt.fewest_switches = .{}; opt.initial_conditions.state = 2;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_tripleState1D_2" {
    var pop = try Vector(f64).init(3, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.79; (try pop.ptr(1)).* = 0.13; (try pop.ptr(2)).* = 0.08;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "tripleState1D_2"; opt.landau_zener = .{}; opt.initial_conditions.state = 2;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_tripleState1D_3" {
    var pop = try Vector(f64).init(3, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.14; (try pop.ptr(1)).* = 0.82; (try pop.ptr(2)).* = 0.04;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "tripleState1D_3"; opt.fewest_switches = .{}; opt.initial_conditions.state = 1;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_tripleState1D_3" {
    var pop = try Vector(f64).init(3, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.12; (try pop.ptr(1)).* = 0.83; (try pop.ptr(2)).* = 0.05;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "tripleState1D_3"; opt.landau_zener = .{}; opt.initial_conditions.state = 1;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_tully1D_1" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.41; (try pop.ptr(1)).* = 0.59;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "tully1D_1"; opt.fewest_switches = .{}; opt.initial_conditions.state = 1;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_tully1D_1" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.51; (try pop.ptr(1)).* = 0.49;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "tully1D_1"; opt.landau_zener = .{}; opt.initial_conditions.state = 1;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_tully1D_2" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.15; (try pop.ptr(1)).* = 0.85;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "tully1D_2"; opt.fewest_switches = .{}; opt.initial_conditions.state = 1;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_tully1D_2" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.13; (try pop.ptr(1)).* = 0.87;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "tully1D_2"; opt.landau_zener = .{}; opt.initial_conditions.state = 1;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_fssh_tully1D_3" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0.58; (try pop.ptr(1)).* = 0.42;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "tully1D_3"; opt.fewest_switches = .{}; opt.initial_conditions.state = 1;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}

test "cdyn_lzsh_tully1D_3" {
    var pop = try Vector(f64).init(2, std.testing.allocator); defer pop.deinit();

    (try pop.ptr(0)).* = 0; (try pop.ptr(1)).* = 1;

    var opt = cdn.ClassicalDynamicsOptions(f64){}; opt.potential = "tully1D_3"; opt.landau_zener = .{}; opt.initial_conditions.state = 1;

    const output = try cdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(vec.eq(f64, output.pop, pop, 1e-12));
}
