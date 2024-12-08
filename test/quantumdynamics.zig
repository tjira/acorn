const std = @import("std");

const mat = @import("acorn").mat;
const qdn = @import("acorn").qdn;

const Matrix = @import("acorn").Matrix;

test "qdyn_doubleState1D_1" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.68939163804408; P.ptr(0, 1).* = 0.00007484823677;
    P.ptr(1, 0).* = 0.00007484823677; P.ptr(1, 1).* = 0.31060836195428;

    var opt = qdn.QuantumDynamicsOptions(f64){}; opt.potential = "doubleState1D_1"; opt.initial_conditions.state = 1;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_doubleState1D_2" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.96070275892968; P.ptr(0, 1).* = 0.00727393248796;
    P.ptr(1, 0).* = 0.00727393248796; P.ptr(1, 1).* = 0.03929724106854;

    var opt = qdn.QuantumDynamicsOptions(f64){}; opt.potential = "doubleState1D_2"; opt.initial_conditions.state = 1;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_tripleState1D_1" {
    var P = try Matrix(f64).init(3, 3, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* =  0.09406137151353; P.ptr(0, 1).* = -0.01173494963317; P.ptr(0, 2).* =  0.00173372711838;
    P.ptr(1, 0).* = -0.01173494963317; P.ptr(1, 1).* =  0.09966537676404; P.ptr(1, 2).* =  0.02130317470369;
    P.ptr(2, 0).* =  0.00173372711838; P.ptr(2, 1).* =  0.02130317470369; P.ptr(2, 2).* =  0.80627325172103;

    var opt = qdn.QuantumDynamicsOptions(f64){}; opt.potential = "tripleState1D_1"; opt.initial_conditions.state = 2;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_tripleState1D_2" {
    var P = try Matrix(f64).init(3, 3, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* =  0.80243373396785; P.ptr(0, 1).* = -0.04371707032137; P.ptr(0, 2).* =  0.01971557113131;
    P.ptr(1, 0).* = -0.04371707032137; P.ptr(1, 1).* =  0.04418760360549; P.ptr(1, 2).* = -0.02345042333257;
    P.ptr(2, 0).* =  0.01971557113131; P.ptr(2, 1).* = -0.02345042333257; P.ptr(2, 2).* =  0.15337866242496;

    var opt = qdn.QuantumDynamicsOptions(f64){}; opt.potential = "tripleState1D_2"; opt.initial_conditions.state = 2;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_tripleState1D_3" {
    var P = try Matrix(f64).init(3, 3, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* =  0.72007774176220; P.ptr(0, 1).* = -0.00000000057348; P.ptr(0, 2).* =  0.00000009131136;
    P.ptr(1, 0).* = -0.00000000057348; P.ptr(1, 1).* =  0.11247117674074; P.ptr(1, 2).* = -0.00000000354249;
    P.ptr(2, 0).* =  0.00000009131136; P.ptr(2, 1).* = -0.00000000354249; P.ptr(2, 2).* =  0.16745108149516;

    var opt = qdn.QuantumDynamicsOptions(f64){}; opt.potential = "tripleState1D_3"; opt.initial_conditions.state = 1;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_tully1D_1" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* =  0.41050635034636; P.ptr(0, 1).* = -0.02165467199410;
    P.ptr(1, 0).* = -0.02165467199410; P.ptr(1, 1).* =  0.58949364965213;

    var opt = qdn.QuantumDynamicsOptions(f64){}; opt.potential = "tully1D_1"; opt.initial_conditions.state = 1;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_tully1D_2" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* =  0.02483440861074; P.ptr(0, 1).* = -0.00009186833299;
    P.ptr(1, 0).* = -0.00009186833299; P.ptr(1, 1).* =  0.97516559138780;

    var opt = qdn.QuantumDynamicsOptions(f64){}; opt.potential = "tully1D_2"; opt.initial_conditions.state = 1;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_tully1D_3" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.76438466406603; P.ptr(0, 1).* = 0.09347537494805;
    P.ptr(1, 0).* = 0.09347537494805; P.ptr(1, 1).* = 0.23561533593184;

    var opt = qdn.QuantumDynamicsOptions(f64){}; opt.potential = "tully1D_3"; opt.initial_conditions.state = 1;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}
