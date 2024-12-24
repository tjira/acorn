const std = @import("std"); const Complex = std.math.Complex;

const mat = @import("matrix.zig");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

pub fn ModelPotentialOptions(comptime T: type) type {
    return struct {
        const ValuePair = struct {
            index: u32 = 0, value: T = 0
        };

        adiabatic: bool = false,
        limits: []const T = &[_]T{-16, 16},
        output: []const u8 = "POTENTIAL.mat",
        points: u32 = 1024,
        potential: []const u8 = "tully1D_1",

        constant: []const ValuePair = &[_]ValuePair{}
    };
}

pub fn dims(potential: []const u8) !u32 {
    if (std.mem.eql(u8, potential,    "harmonic1D_1")) return 1;
    if (std.mem.eql(u8, potential,    "harmonic2D_1")) return 2;
    if (std.mem.eql(u8, potential,    "harmonic3D_1")) return 3;
    if (std.mem.eql(u8, potential, "doubleState1D_1")) return 1;
    if (std.mem.eql(u8, potential, "doubleState1D_2")) return 1;
    if (std.mem.eql(u8, potential, "tripleState1D_1")) return 1;
    if (std.mem.eql(u8, potential, "tripleState1D_2")) return 1;
    if (std.mem.eql(u8, potential, "tripleState1D_3")) return 1;
    if (std.mem.eql(u8, potential,       "tully1D_1")) return 1;
    if (std.mem.eql(u8, potential,       "tully1D_2")) return 1;
    if (std.mem.eql(u8, potential,       "tully1D_3")) return 1;

    if (std.mem.eql(u8, potential,      "uracil1D_1")) return 8;

    return error.UnknownPotential;
}

pub fn eval(comptime T: type, U: *Matrix(T), potential: []const u8, r: Vector(T)) !void {
    if (std.mem.eql(u8, potential,    "harmonic1D_1"))    return harmonic1D_1(T, U, r);
    if (std.mem.eql(u8, potential,    "harmonic2D_1"))    return harmonic2D_1(T, U, r);
    if (std.mem.eql(u8, potential,    "harmonic3D_1"))    return harmonic3D_1(T, U, r);
    if (std.mem.eql(u8, potential, "doubleState1D_1")) return doubleState1D_1(T, U, r);
    if (std.mem.eql(u8, potential, "doubleState1D_2")) return doubleState1D_2(T, U, r);
    if (std.mem.eql(u8, potential, "tripleState1D_1")) return tripleState1D_1(T, U, r);
    if (std.mem.eql(u8, potential, "tripleState1D_2")) return tripleState1D_2(T, U, r);
    if (std.mem.eql(u8, potential, "tripleState1D_3")) return tripleState1D_3(T, U, r);
    if (std.mem.eql(u8, potential,       "tully1D_1"))       return tully1D_1(T, U, r);
    if (std.mem.eql(u8, potential,       "tully1D_2"))       return tully1D_2(T, U, r);
    if (std.mem.eql(u8, potential,       "tully1D_3"))       return tully1D_3(T, U, r);

    if (std.mem.eql(u8, potential,      "uracil1D_1"))      return uracil1D_1(T, U, r);

    return error.UnknownPotential;
}

pub fn states(potential: []const u8) !u32 {
    if (std.mem.eql(u8, potential,    "harmonic1D_1")) return 1;
    if (std.mem.eql(u8, potential,    "harmonic2D_1")) return 1;
    if (std.mem.eql(u8, potential,    "harmonic3D_1")) return 1;
    if (std.mem.eql(u8, potential, "doubleState1D_1")) return 2;
    if (std.mem.eql(u8, potential, "doubleState1D_2")) return 2;
    if (std.mem.eql(u8, potential, "tripleState1D_1")) return 3;
    if (std.mem.eql(u8, potential, "tripleState1D_2")) return 3;
    if (std.mem.eql(u8, potential, "tripleState1D_3")) return 3;
    if (std.mem.eql(u8, potential,       "tully1D_1")) return 2;
    if (std.mem.eql(u8, potential,       "tully1D_2")) return 2;
    if (std.mem.eql(u8, potential,       "tully1D_3")) return 2;

    if (std.mem.eql(u8, potential,      "uracil1D_1")) return 4;

    return error.UnknownPotential;
}

pub fn harmonic1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.5 * r.at(0) * r.at(0);
}

pub fn harmonic2D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.5 * (r.at(0) * r.at(0) + r.at(1) * r.at(1));
}

pub fn harmonic3D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.5 * (r.at(0) * r.at(0) + r.at(1) * r.at(1) + r.at(2) * r.at(2));
}

pub fn doubleState1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.001 * r.at(0);
    U.ptr(0, 1).* = 0.001 * std.math.exp(-0.05 * r.at(0) * r.at(0));

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -0.001 * r.at(0);
}

pub fn doubleState1D_2(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.01 * std.math.tanh(0.6 * r.at(0));
    U.ptr(0, 1).* = 0.001 * std.math.exp(-r.at(0) * r.at(0));

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -0.01 * std.math.tanh(0.6 * r.at(0));
}

pub fn tripleState1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.001 * r.at(0);
    U.ptr(0, 1).* = 0.001 * std.math.exp(-0.01 * r.at(0) * r.at(0));
    U.ptr(0, 2).* = 0.002 * std.math.exp(-0.01 * r.at(0) * r.at(0));

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = 0;
    U.ptr(1, 2).* = 0.001 * std.math.exp(-0.01 * r.at(0) * r.at(0));

    U.ptr(2, 0).* = U.at(0, 2);
    U.ptr(2, 1).* = U.at(1, 2);
    U.ptr(2, 2).* = -0.001 * r.at(0);
}

pub fn tripleState1D_2(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.01 * std.math.tanh(0.5 * r.at(0));
    U.ptr(0, 1).* = 0.001 * std.math.exp(-r.at(0) * r.at(0));
    U.ptr(0, 2).* = 0.002 * std.math.exp(-r.at(0) * r.at(0));

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = 0;
    U.ptr(1, 2).* = 0.001 * std.math.exp(-r.at(0) * r.at(0));

    U.ptr(2, 0).* = U.at(0, 2);
    U.ptr(2, 1).* = U.at(1, 2);
    U.ptr(2, 2).* = -0.01 * std.math.tanh(0.5 * r.at(0));
}

pub fn tripleState1D_3(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.03 * (std.math.tanh(1.6 * r.at(0)) + std.math.tanh(1.6 * (r.at(0) + 7)));
    U.ptr(0, 1).* = 0.005 * std.math.exp(-r.at(0) * r.at(0));
    U.ptr(0, 2).* = 0.005 * std.math.exp(-(r.at(0) + 7) * (r.at(0) + 7));

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -0.03 * (std.math.tanh(1.6 * r.at(0)) + std.math.tanh(1.6 * (r.at(0) - 7)));
    U.ptr(1, 2).* = 0.005 * std.math.exp(-(r.at(0) - 7) * (r.at(0) - 7));

    U.ptr(2, 0).* = U.at(0, 2);
    U.ptr(2, 1).* = U.at(1, 2);
    U.ptr(2, 2).* = -0.03 * (std.math.tanh(1.6 * (r.at(0) + 7)) - std.math.tanh(1.6 * (r.at(0) - 7)));
}

pub fn tully1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = if (r.at(0) > 0) 0.01 * (1 - std.math.exp(-1.6 * r.at(0))) else -0.01 * (1 - std.math.exp(1.6 * r.at(0)));
    U.ptr(0, 1).* = 0.005 * std.math.exp(-r.at(0) * r.at(0));

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -U.at(0, 0);
}

pub fn tully1D_2(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0;
    U.ptr(0, 1).* = 0.0150 * std.math.exp(-0.06 * r.at(0) * r.at(0));

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -0.1 * std.math.exp(-0.28 * r.at(0) * r.at(0)) + 0.05;
}

pub fn tully1D_3(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 6e-4;
    U.ptr(0, 1).* = if (r.at(0) > 0) 0.1 * (2 - std.math.exp(-0.9 * r.at(0))) else 0.1 * std.math.exp(0.9 * r.at(0));

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -U.at(0, 0);
}

pub fn uracil1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    const au2cm = 219474.63068; const au2ev = 27.21138602;

    const E_0 = 9.59; const E_1 = 10.11; const E_2 = 10.48; const E_3 = 11.08;

    const omg_10 = 734.0; const l_10_01 = 0.04633; const l_10_12 = 0.03148;
    const omg_12 = 771.0; const l_12_01 = 0.03540; const l_12_12 = 0.03607;

    const omg_18 = 1193.0; const k_18_0 = -0.02203; const k_18_1 = 0.09074; const k_18_2 = 0.02748; const k_18_3 = -0.04054; const l_18_02 = -0.03538; const l_18_13 = 0.08077; const g_18_0 = 0.01938; const g_18_1 = 0.00694; const g_18_2 = -0.00294; const g_18_3 = 0.00752;
    const omg_20 = 1383.0; const k_20_0 = -0.12147; const k_20_1 = 0.05316; const k_20_2 = 0.11233; const k_20_3 =  0.00747; const l_20_02 = -0.02049; const l_20_13 = 0.00000; const g_20_0 = 0.01489; const g_20_1 = 0.00828; const g_20_2 =  0.00183; const g_20_3 = 0.00546;
    const omg_21 = 1406.0; const k_21_0 = -0.09468; const k_21_1 = 0.04454; const k_21_2 = 0.14539; const k_21_3 =  0.00050; const l_21_02 =  0.00000; const l_21_13 = 0.07284; const g_21_0 = 0.00970; const g_21_1 = 0.00096; const g_21_2 = -0.00114; const g_21_3 = 0.01108;

    const l_24_02 = 0.0000; const l_24_13 = -0.18132;
    const d_24_0 = 41.89704; const a_24_0 =  0.04719; const q_24_0 =  0.81440; const e_24_0 = -0.06431;
    const d_24_1 = 38.37122; const a_24_1 =  0.05231; const q_24_1 =  0.37488; const e_24_1 = -0.01505;
    const d_24_2 = 39.25691; const a_24_2 =  0.05286; const q_24_2 =  0.14859; const e_24_2 = -0.00244;
    const d_24_3 = 37.97847; const a_24_3 =  0.05431; const q_24_3 = -0.18152; const e_24_3 = -0.00366;

    const l_25_02 = 0.00114; const l_25_13 = 0.12606;
    const d_25_0 =  4.80270; const a_25_0 = -0.13675; const q_25_0 =  0.02883; const e_25_0 = -0.00007;
    const d_25_1 = 74.15995; const a_25_1 = -0.03064; const q_25_1 = -1.34468; const e_25_1 = -0.12082;
    const d_25_2 = 90.76928; const a_25_2 = -0.03374; const q_25_2 = -0.29923; const e_25_2 = -0.00916;
    const d_25_3 = 20.56979; const a_25_3 = -0.08044; const q_25_3 =  0.38841; const e_25_3 = -0.02071;

    const l_26_02 = 0.13035; const l_26_13 = 0.14272;
    const d_26_0 = 22.92802; const a_26_0 =  0.07438; const q_26_0 = -0.32069; const e_26_0 = -0.01274;
    const d_26_1 = 18.27440; const a_26_1 =  0.07911; const q_26_1 = -0.01711; const e_26_1 = -0.00003;
    const d_26_2 =  9.46894; const a_26_2 =  0.08653; const q_26_2 =  0.37635; const e_26_2 = -0.01037;
    const d_26_3 = 65.09678; const a_26_3 =  0.03660; const q_26_3 =  1.66312; const e_26_3 = -0.25639;

    U.ptr(0, 0).* = E_0 / au2ev + 0.5 * (omg_10 * r.at(0) * r.at(0) + omg_12 * r.at(1) * r.at(1) + omg_18 * r.at(2) * r.at(2) + omg_20 * r.at(3) * r.at(3) + omg_21 * r.at(4) * r.at(4)) / au2cm + ((d_24_0 * (std.math.exp(a_24_0 * (r.at(5) - q_24_0)) - 1) * (std.math.exp(a_24_0 * (r.at(5) - q_24_0)) - 1) + e_24_0) + (d_25_0 * (std.math.exp(a_25_0 * (r.at(6) - q_25_0)) - 1) * (std.math.exp(a_25_0 * (r.at(6) - q_25_0)) - 1) + e_25_0) + (d_26_0 * (std.math.exp(a_26_0 * (r.at(7) - q_26_0)) - 1) * (std.math.exp(a_26_0 * (r.at(7) - q_26_0)) - 1) + e_26_0)) / au2ev + (k_18_0 * r.at(2) + k_20_0 * r.at(3) + k_21_0 * r.at(4)) / au2ev + 0.5 * (g_18_0 * r.at(2) * r.at(2) + g_20_0 * r.at(3) * r.at(3) + g_21_0 * r.at(4) * r.at(4)) / au2ev;
    U.ptr(0, 1).* = (l_10_01 * r.at(0) + l_12_01 * r.at(1)) / au2ev;
    U.ptr(0, 2).* = (l_18_02 * r.at(2) + l_20_02 * r.at(3) + l_21_02 * r.at(4) + l_24_02 * r.at(5) + l_25_02 * r.at(6) + l_26_02 * r.at(7)) / au2ev;
    U.ptr(0, 3).* = 0;

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = E_1 / au2ev + 0.5 * (omg_10 * r.at(0) * r.at(0) + omg_12 * r.at(1) * r.at(1) + omg_18 * r.at(2) * r.at(2) + omg_20 * r.at(3) * r.at(3) + omg_21 * r.at(4) * r.at(4)) / au2cm + ((d_24_1 * (std.math.exp(a_24_1 * (r.at(5) - q_24_1)) - 1) * (std.math.exp(a_24_1 * (r.at(5) - q_24_1)) - 1) + e_24_1) + (d_25_1 * (std.math.exp(a_25_1 * (r.at(6) - q_25_1)) - 1) * (std.math.exp(a_25_1 * (r.at(6) - q_25_1)) - 1) + e_25_1) + (d_26_1 * (std.math.exp(a_26_1 * (r.at(7) - q_26_1)) - 1) * (std.math.exp(a_26_1 * (r.at(7) - q_26_1)) - 1) + e_26_1)) / au2ev + (k_18_1 * r.at(2) + k_20_1 * r.at(3) + k_21_1 * r.at(4)) / au2ev + 0.5 * (g_18_1 * r.at(2) * r.at(2) + g_20_1 * r.at(3) * r.at(3) + g_21_1 * r.at(4) * r.at(4)) / au2ev;
    U.ptr(1, 2).* = (l_10_12 * r.at(0) + l_12_12 * r.at(1)) / au2ev;
    U.ptr(1, 3).* = (l_18_13 * r.at(2) + l_20_13 * r.at(3) + l_21_13 * r.at(4) + l_24_13 * r.at(5) + l_25_13 * r.at(6) + l_26_13 * r.at(7)) / au2ev;

    U.ptr(2, 0).* = U.at(0, 2);
    U.ptr(2, 1).* = U.at(1, 2);
    U.ptr(2, 2).* = E_2 / au2ev + 0.5 * (omg_10 * r.at(0) * r.at(0) + omg_12 * r.at(1) * r.at(1) + omg_18 * r.at(2) * r.at(2) + omg_20 * r.at(3) * r.at(3) + omg_21 * r.at(4) * r.at(4)) / au2cm + ((d_24_2 * (std.math.exp(a_24_2 * (r.at(5) - q_24_2)) - 1) * (std.math.exp(a_24_2 * (r.at(5) - q_24_2)) - 1) + e_24_2) + (d_25_2 * (std.math.exp(a_25_2 * (r.at(6) - q_25_2)) - 1) * (std.math.exp(a_25_2 * (r.at(6) - q_25_2)) - 1) + e_25_2) + (d_26_2 * (std.math.exp(a_26_2 * (r.at(7) - q_26_2)) - 1) * (std.math.exp(a_26_2 * (r.at(7) - q_26_2)) - 1) + e_26_2)) / au2ev + (k_18_2 * r.at(2) + k_20_2 * r.at(3) + k_21_2 * r.at(4)) / au2ev + 0.5 * (g_18_2 * r.at(2) * r.at(2) + g_20_2 * r.at(3) * r.at(3) + g_21_2 * r.at(4) * r.at(4)) / au2ev;
    U.ptr(2, 3).* = 0;

    U.ptr(3, 0).* = U.at(0, 3);
    U.ptr(3, 1).* = U.at(1, 3);
    U.ptr(3, 2).* = U.at(2, 3);
    U.ptr(3, 3).* = E_3 / au2ev + 0.5 * (omg_10 * r.at(0) * r.at(0) + omg_12 * r.at(1) * r.at(1) + omg_18 * r.at(2) * r.at(2) + omg_20 * r.at(3) * r.at(3) + omg_21 * r.at(4) * r.at(4)) / au2cm + ((d_24_3 * (std.math.exp(a_24_3 * (r.at(5) - q_24_3)) - 1) * (std.math.exp(a_24_3 * (r.at(5) - q_24_3)) - 1) + e_24_3) + (d_25_3 * (std.math.exp(a_25_3 * (r.at(6) - q_25_3)) - 1) * (std.math.exp(a_25_3 * (r.at(6) - q_25_3)) - 1) + e_25_3) + (d_26_3 * (std.math.exp(a_26_3 * (r.at(7) - q_26_3)) - 1) * (std.math.exp(a_26_3 * (r.at(7) - q_26_3)) - 1) + e_26_3)) / au2ev + (k_18_3 * r.at(2) + k_20_3 * r.at(3) + k_21_3 * r.at(4)) / au2ev + 0.5 * (g_18_3 * r.at(2) * r.at(2) + g_20_3 * r.at(3) * r.at(3) + g_21_3 * r.at(4) * r.at(4)) / au2ev;
}

pub fn kgrid(comptime T: type, k: *Matrix(T), start: T, end: T, points: u32) void {
    k.fill(2 * std.math.pi / asfloat(T, points) / (end - start) * asfloat(T, points - 1));

    for (0..k.rows / points) |i| for (0..points) |j| {
        k.ptr(i * points + j, k.cols - 1).* *= asfloat(T, j) - asfloat(T, if (j < points / 2) 0 else points);
    };

    for (0..k.rows) |i| for (0..k.cols) |j| {
        k.ptr(i, j).* = k.at(i / std.math.pow(usize, points, k.cols - j - 1), k.cols - 1);
    };
}

pub fn rgrid(comptime T: type, r: *Matrix(T), start: T, end: T, points: u32) void {
    for (0..r.rows) |i| for (0..r.cols) |j| {
        r.ptr(i, j).* = start + asfloat(T, i / std.math.pow(usize, points, r.cols - j - 1) % points) * (end - start) / asfloat(T, points - 1);
    };
}

pub fn write(comptime T: type, opt: ModelPotentialOptions(T), allocator: std.mem.Allocator) !void {
    const ndim = try dims(opt.potential); const nstate = try states(opt.potential);

    var T1 = try Matrix(T).init(nstate, nstate, allocator); defer T1.deinit();
    var T2 = try Matrix(T).init(nstate, nstate, allocator); defer T2.deinit();

    var U  = try Matrix(T).init(nstate, nstate, allocator); defer  U.deinit();
    var UA = try Matrix(T).init(nstate, nstate, allocator); defer UA.deinit();
    var UC = try Matrix(T).init(nstate, nstate, allocator); defer UC.deinit();

    var R = try Matrix(T).init(std.math.pow(u32, opt.points, ndim - @as(u32, @intCast(opt.constant.len))), ndim - opt.constant.len, allocator); defer R.deinit();
    var V = try Matrix(T).init(std.math.pow(u32, opt.points, ndim - @as(u32, @intCast(opt.constant.len))), nstate * nstate,         allocator); defer V.deinit();

    rgrid(T, &R, opt.limits[0], opt.limits[1], opt.points); var r = try Vector(T).init(ndim, allocator); defer r.deinit(); 

    for (0..R.rows) |i| {

        var l: u32 = 0;

        for (0..ndim) |j| {

            var constant = false;

            for (0..opt.constant.len) |k| if (opt.constant[k].index == j) {
                r.ptr(j).* = opt.constant[k].value; constant = true;
            };

            if (!constant) {
                r.ptr(j).* = R.at(i, l); l += 1;
            }
        }

        try eval(T, &U, opt.potential, r);

        if (opt.adiabatic) {
            mat.eigh(T, &UA, &UC, U, &T1, &T2); @memcpy(U.data, UA.data);
        }

        for (U.data, 0..) |e, j| V.ptr(i, j).* = e;
    }

    var VT = try Matrix(T).init(V.rows, R.cols + V.cols, allocator); mat.hjoin(T, &VT, R, V); try VT.write(opt.output); VT.deinit();
}
