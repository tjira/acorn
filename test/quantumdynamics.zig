const std = @import("std");

const mat = @import("acorn").mat;
const qdn = @import("acorn").qdn;

const Matrix = @import("acorn").Matrix;

test "qdyn_imag_harmonic1D_1" {
    var opt = qdn.QuantumDynamicsOptions(f64){};

    opt.potential = "harmonic1D_1"; opt.grid.limits = &[_]f64{-8, 8}; opt.grid.points = 32;

    opt.initial_conditions.state = 0; opt.initial_conditions.position = &[_]f64{1}; opt.initial_conditions.momentum = &[_]f64{0}; opt.initial_conditions.mass = 1;

    opt.adiabatic = false; opt.iterations = 300; opt.mode = &[_]u32{1, 0}; opt.time_step = 0.1;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(@abs(output.Ekin + output.Epot - 0.50000038965072) < 1e-12);
}

test "qdyn_imag_harmonic2D_1" {
    var opt = qdn.QuantumDynamicsOptions(f64){};

    opt.potential = "harmonic2D_1"; opt.grid.limits = &[_]f64{-8, 8}; opt.grid.points = 32;

    opt.initial_conditions.state = 0; opt.initial_conditions.position = &[_]f64{1, 1}; opt.initial_conditions.momentum = &[_]f64{0, 0}; opt.initial_conditions.mass = 1;

    opt.adiabatic = false; opt.iterations = 300; opt.mode = &[_]u32{1, 0}; opt.time_step = 0.1;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(@abs(output.Ekin + output.Epot - 1.00000077930144) < 1e-12);
}

test "qdyn_imag_harmonic3D_1" {
    var opt = qdn.QuantumDynamicsOptions(f64){};

    opt.potential = "harmonic3D_1"; opt.grid.limits = &[_]f64{-8, 8}; opt.grid.points = 32;

    opt.initial_conditions.state = 0; opt.initial_conditions.position = &[_]f64{1, 1, 1}; opt.initial_conditions.momentum = &[_]f64{0, 0, 0}; opt.initial_conditions.mass = 1;

    opt.adiabatic = false; opt.iterations = 300; opt.mode = &[_]u32{1, 0}; opt.time_step = 0.1;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(@abs(output.Ekin + output.Epot - 1.50000116895219) < 1e-12);
}

test "qdyn_real_doubleState1D_1" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* =  0.68939163804408; P.ptr(0, 1).* = -0.00007484823676;
    P.ptr(1, 0).* = -0.00007484823676; P.ptr(1, 1).* =  0.31060836195428;

    var opt = qdn.QuantumDynamicsOptions(f64){};

    opt.potential = "doubleState1D_1"; opt.grid.limits = &[_]f64{-16, 32}; opt.grid.points = 512;

    opt.initial_conditions.state = 1; opt.initial_conditions.position = &[_]f64{-10}; opt.initial_conditions.momentum = &[_]f64{15}; opt.initial_conditions.mass = 2000;

    opt.adiabatic = true; opt.iterations = 300; opt.mode = &[_]u32{0, 1}; opt.time_step = 10;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_real_doubleState1D_2" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* =  0.96070275892968; P.ptr(0, 1).* = -0.00727393248796;
    P.ptr(1, 0).* = -0.00727393248796; P.ptr(1, 1).* =  0.03929724106854;

    var opt = qdn.QuantumDynamicsOptions(f64){};

    opt.potential = "doubleState1D_2"; opt.grid.limits = &[_]f64{-16, 32}; opt.grid.points = 512;

    opt.initial_conditions.state = 1; opt.initial_conditions.position = &[_]f64{-10}; opt.initial_conditions.momentum = &[_]f64{15}; opt.initial_conditions.mass = 2000;

    opt.adiabatic = true; opt.iterations = 300; opt.mode = &[_]u32{0, 1}; opt.time_step = 10;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_real_tripleState1D_1" {
    var P = try Matrix(f64).init(3, 3, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* =  0.47168386547338; P.ptr(0, 1).* = -0.05402066588886; P.ptr(0, 2).* =  0.00017180750094;
    P.ptr(1, 0).* = -0.05402066588886; P.ptr(1, 1).* =  0.42854789634051; P.ptr(1, 2).* = -0.03148902264056;
    P.ptr(2, 0).* =  0.00017180750094; P.ptr(2, 1).* = -0.03148902264056; P.ptr(2, 2).* =  0.09976823818447;

    var opt = qdn.QuantumDynamicsOptions(f64){};

    opt.potential = "tripleState1D_1"; opt.grid.limits = &[_]f64{-16, 32}; opt.grid.points = 512;

    opt.initial_conditions.state = 2; opt.initial_conditions.position = &[_]f64{-10}; opt.initial_conditions.momentum = &[_]f64{15}; opt.initial_conditions.mass = 2000;

    opt.adiabatic = true; opt.iterations = 300; opt.mode = &[_]u32{0, 1}; opt.time_step = 10;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_real_tripleState1D_2" {
    var P = try Matrix(f64).init(3, 3, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* =  0.95479039377455; P.ptr(0, 1).* =  0.06561014255438; P.ptr(0, 2).* = -0.00059884544352;
    P.ptr(1, 0).* =  0.06561014255438; P.ptr(1, 1).* =  0.04466241929869; P.ptr(1, 2).* =  0.00181568734195;
    P.ptr(2, 0).* = -0.00059884544352; P.ptr(2, 1).* =  0.00181568734195; P.ptr(2, 2).* =  0.00054718692497;

    var opt = qdn.QuantumDynamicsOptions(f64){};

    opt.potential = "tripleState1D_2"; opt.grid.limits = &[_]f64{-16, 32}; opt.grid.points = 512;

    opt.initial_conditions.state = 2; opt.initial_conditions.position = &[_]f64{-10}; opt.initial_conditions.momentum = &[_]f64{15}; opt.initial_conditions.mass = 2000;

    opt.adiabatic = true; opt.iterations = 300; opt.mode = &[_]u32{0, 1}; opt.time_step = 10;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_real_tripleState1D_3" {
    var P = try Matrix(f64).init(3, 3, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* =  0.72007774176220; P.ptr(0, 1).* = -0.00000000050304; P.ptr(0, 2).* = 0.00000004962809;
    P.ptr(1, 0).* = -0.00000000050304; P.ptr(1, 1).* =  0.11247117674074; P.ptr(1, 2).* = 0.00000000354232;
    P.ptr(2, 0).* =  0.00000004962809; P.ptr(2, 1).* =  0.00000000354232; P.ptr(2, 2).* = 0.16745108149516;

    var opt = qdn.QuantumDynamicsOptions(f64){};

    opt.potential = "tripleState1D_3"; opt.grid.limits = &[_]f64{-16, 32}; opt.grid.points = 512;

    opt.initial_conditions.state = 1; opt.initial_conditions.position = &[_]f64{-10}; opt.initial_conditions.momentum = &[_]f64{15}; opt.initial_conditions.mass = 2000;

    opt.adiabatic = true; opt.iterations = 300; opt.mode = &[_]u32{0, 1}; opt.time_step = 10;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_real_tully1D_1" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.41050635034636; P.ptr(0, 1).* = 0.02165466513650;
    P.ptr(1, 0).* = 0.02165466513650; P.ptr(1, 1).* = 0.58949364965213;

    var opt = qdn.QuantumDynamicsOptions(f64){};

    opt.potential = "tully1D_1"; opt.grid.limits = &[_]f64{-16, 32}; opt.grid.points = 512;

    opt.initial_conditions.state = 1; opt.initial_conditions.position = &[_]f64{-10}; opt.initial_conditions.momentum = &[_]f64{15}; opt.initial_conditions.mass = 2000;

    opt.adiabatic = true; opt.iterations = 300; opt.mode = &[_]u32{0, 1}; opt.time_step = 10;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_real_tully1D_2" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* =  0.02483440861074; P.ptr(0, 1).* = -0.00009186830172;
    P.ptr(1, 0).* = -0.00009186830172; P.ptr(1, 1).* =  0.97516559138780;

    var opt = qdn.QuantumDynamicsOptions(f64){};

    opt.potential = "tully1D_2"; opt.grid.limits = &[_]f64{-16, 32}; opt.grid.points = 512;

    opt.initial_conditions.state = 1; opt.initial_conditions.position = &[_]f64{-10}; opt.initial_conditions.momentum = &[_]f64{15}; opt.initial_conditions.mass = 2000;

    opt.adiabatic = true; opt.iterations = 300; opt.mode = &[_]u32{0, 1}; opt.time_step = 10;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}

test "qdyn_real_tully1D_3" {
    var P = try Matrix(f64).init(2, 2, std.testing.allocator); defer P.deinit();

    P.ptr(0, 0).* = 0.76438466406603; P.ptr(0, 1).* = 0.09347537494805;
    P.ptr(1, 0).* = 0.09347537494805; P.ptr(1, 1).* = 0.23561533593184;

    var opt = qdn.QuantumDynamicsOptions(f64){};

    opt.potential = "tully1D_3"; opt.grid.limits = &[_]f64{-16, 32}; opt.grid.points = 512;

    opt.initial_conditions.state = 1; opt.initial_conditions.position = &[_]f64{-10}; opt.initial_conditions.momentum = &[_]f64{15}; opt.initial_conditions.mass = 2000;

    opt.adiabatic = true; opt.iterations = 300; opt.mode = &[_]u32{0, 1}; opt.time_step = 10;

    const output = try qdn.run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(mat.eq(f64, output.P, P, 1e-12));
}
