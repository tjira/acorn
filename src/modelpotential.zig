//! Module to store the potential energy surfaces of the model Hamiltonians.

const std = @import("std"); const Complex = std.math.Complex;

const cnt = @import("constant.zig");
const mat = @import("matrix.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

/// Option struct for the model potential run target.
pub fn ModelPotentialOptions(comptime T: type) type {
    return struct {
        pub const ValuePair = struct {
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

/// Return the number of dimensions of the potential energy surface.
pub fn dims(potential: []const u8) !u32 {
    if (std.mem.eql(u8, potential,    "harmonic1D_1"))  return 1;
    if (std.mem.eql(u8, potential,    "harmonic2D_1"))  return 2;
    if (std.mem.eql(u8, potential,    "harmonic3D_1"))  return 3;
    if (std.mem.eql(u8, potential,    "harmonic4D_1"))  return 4;
    if (std.mem.eql(u8, potential, "doubleState1D_1"))  return 1;
    if (std.mem.eql(u8, potential, "doubleState1D_2"))  return 1;
    if (std.mem.eql(u8, potential, "tripleState1D_1"))  return 1;
    if (std.mem.eql(u8, potential, "tripleState1D_2"))  return 1;
    if (std.mem.eql(u8, potential, "tripleState1D_3"))  return 1;
    if (std.mem.eql(u8, potential,       "tully1D_1"))  return 1;
    if (std.mem.eql(u8, potential,       "tully1D_2"))  return 1;
    if (std.mem.eql(u8, potential,       "tully1D_3"))  return 1;

    if (std.mem.eql(u8, potential,      "uracil8D_1")) return  8;
    if (std.mem.eql(u8, potential,     "uracil12D_1")) return 12;

    return error.UnknownPotential;
}

/// Evaluate the potential energy surface at a given position. The result is stored in the matrix U.
pub fn eval(comptime T: type, U: *Matrix(T), potential: []const u8, r: Vector(T)) !void {
    if (std.mem.eql(u8, potential,    "harmonic1D_1"))    return harmonic1D_1(T, U, r);
    if (std.mem.eql(u8, potential,    "harmonic2D_1"))    return harmonic2D_1(T, U, r);
    if (std.mem.eql(u8, potential,    "harmonic3D_1"))    return harmonic3D_1(T, U, r);
    if (std.mem.eql(u8, potential,    "harmonic4D_1"))    return harmonic4D_1(T, U, r);
    if (std.mem.eql(u8, potential, "doubleState1D_1")) return doubleState1D_1(T, U, r);
    if (std.mem.eql(u8, potential, "doubleState1D_2")) return doubleState1D_2(T, U, r);
    if (std.mem.eql(u8, potential, "tripleState1D_1")) return tripleState1D_1(T, U, r);
    if (std.mem.eql(u8, potential, "tripleState1D_2")) return tripleState1D_2(T, U, r);
    if (std.mem.eql(u8, potential, "tripleState1D_3")) return tripleState1D_3(T, U, r);
    if (std.mem.eql(u8, potential,       "tully1D_1"))       return tully1D_1(T, U, r);
    if (std.mem.eql(u8, potential,       "tully1D_2"))       return tully1D_2(T, U, r);
    if (std.mem.eql(u8, potential,       "tully1D_3"))       return tully1D_3(T, U, r);

    if (std.mem.eql(u8, potential,      "uracil8D_1"))     return  uracil8D_1(T, U, r);
    if (std.mem.eql(u8, potential,     "uracil12D_1"))     return uracil12D_1(T, U, r);

    return error.UnknownPotential;
}

/// Return the number of states of the potential energy surface.
pub fn states(potential: []const u8) !u32 {
    if (std.mem.eql(u8, potential,    "harmonic1D_1")) return 1;
    if (std.mem.eql(u8, potential,    "harmonic2D_1")) return 1;
    if (std.mem.eql(u8, potential,    "harmonic3D_1")) return 1;
    if (std.mem.eql(u8, potential,    "harmonic4D_1")) return 1;
    if (std.mem.eql(u8, potential, "doubleState1D_1")) return 2;
    if (std.mem.eql(u8, potential, "doubleState1D_2")) return 2;
    if (std.mem.eql(u8, potential, "tripleState1D_1")) return 3;
    if (std.mem.eql(u8, potential, "tripleState1D_2")) return 3;
    if (std.mem.eql(u8, potential, "tripleState1D_3")) return 3;
    if (std.mem.eql(u8, potential,       "tully1D_1")) return 2;
    if (std.mem.eql(u8, potential,       "tully1D_2")) return 2;
    if (std.mem.eql(u8, potential,       "tully1D_3")) return 2;

    if (std.mem.eql(u8, potential,      "uracil8D_1")) return 4;
    if (std.mem.eql(u8, potential,     "uracil12D_1")) return 4;

    return error.UnknownPotential;
}

/// One-dimensional harmonic oscillator potential energy surface.
pub fn harmonic1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.5 * r.at(0) * r.at(0);
}

/// Two-dimensional harmonic oscillator potential energy surface.
pub fn harmonic2D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.5 * (r.at(0) * r.at(0) + r.at(1) * r.at(1));
}

/// Three-dimensional harmonic oscillator potential energy surface.
pub fn harmonic3D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.5 * (r.at(0) * r.at(0) + r.at(1) * r.at(1) + r.at(2) * r.at(2));
}

/// Four-dimensional harmonic oscillator potential energy surface.
pub fn harmonic4D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.5 * (r.at(0) * r.at(0) + r.at(1) * r.at(1) + r.at(2) * r.at(2) + r.at(3) * r.at(3));
}

/// The first double-state model potential energy surface.
pub fn doubleState1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.001 * r.at(0);
    U.ptr(0, 1).* = 0.001 * std.math.exp(-0.05 * r.at(0) * r.at(0));

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -0.001 * r.at(0);
}

/// The second double-state model potential energy surface.
pub fn doubleState1D_2(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.01 * std.math.tanh(0.6 * r.at(0));
    U.ptr(0, 1).* = 0.001 * std.math.exp(-r.at(0) * r.at(0));

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -0.01 * std.math.tanh(0.6 * r.at(0));
}

/// The first triple-state model potential energy surface.
pub fn tripleState1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.001 * r.at(0);
    U.ptr(0, 1).* = 0.001 * std.math.exp(-0.01 * r.at(0) * r.at(0));
    U.ptr(0, 2).* = 0;

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = 0;
    U.ptr(1, 2).* = 0.001 * std.math.exp(-0.01 * r.at(0) * r.at(0));

    U.ptr(2, 0).* = U.at(0, 2);
    U.ptr(2, 1).* = U.at(1, 2);
    U.ptr(2, 2).* = -0.001 * r.at(0);
}

/// The second triple-state model potential energy surface.
pub fn tripleState1D_2(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0.01 * std.math.tanh(0.5 * r.at(0));
    U.ptr(0, 1).* = 0.001 * std.math.exp(-r.at(0) * r.at(0));
    U.ptr(0, 2).* = 0;

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = 0;
    U.ptr(1, 2).* = 0.001 * std.math.exp(-r.at(0) * r.at(0));

    U.ptr(2, 0).* = U.at(0, 2);
    U.ptr(2, 1).* = U.at(1, 2);
    U.ptr(2, 2).* = -0.01 * std.math.tanh(0.5 * r.at(0));
}

/// The third triple-state model potential energy surface.
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

/// The first Tully model potential energy surface.
pub fn tully1D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = if (r.at(0) > 0) 0.01 * (1 - std.math.exp(-1.6 * r.at(0))) else -0.01 * (1 - std.math.exp(1.6 * r.at(0)));
    U.ptr(0, 1).* = 0.005 * std.math.exp(-r.at(0) * r.at(0));

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -U.at(0, 0);
}

/// The second Tully model potential energy surface.
pub fn tully1D_2(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 0;
    U.ptr(0, 1).* = 0.0150 * std.math.exp(-0.06 * r.at(0) * r.at(0));

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -0.1 * std.math.exp(-0.28 * r.at(0) * r.at(0)) + 0.05;
}

/// The third Tully model potential energy surface.
pub fn tully1D_3(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.ptr(0, 0).* = 6e-4;
    U.ptr(0, 1).* = if (r.at(0) > 0) 0.1 * (2 - std.math.exp(-0.9 * r.at(0))) else 0.1 * std.math.exp(0.9 * r.at(0));

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(1, 1).* = -U.at(0, 0);
}

/// The 8-dimensional uracil model potential energy surface.
pub fn uracil8D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.fill(0);

    const omega_all = &[_]T{734, 771, 1193, 1383, 1406, 1673, 1761, 1794};

    const q = &[_]T{
        r.at(0) * std.math.sqrt(omega_all[0] / cnt.EV2RCM / cnt.AU2EV),
        r.at(1) * std.math.sqrt(omega_all[1] / cnt.EV2RCM / cnt.AU2EV),
        r.at(2) * std.math.sqrt(omega_all[2] / cnt.EV2RCM / cnt.AU2EV),
        r.at(3) * std.math.sqrt(omega_all[3] / cnt.EV2RCM / cnt.AU2EV),
        r.at(4) * std.math.sqrt(omega_all[4] / cnt.EV2RCM / cnt.AU2EV),
        r.at(5) * std.math.sqrt(omega_all[5] / cnt.EV2RCM / cnt.AU2EV),
        r.at(6) * std.math.sqrt(omega_all[6] / cnt.EV2RCM / cnt.AU2EV),
        r.at(7) * std.math.sqrt(omega_all[7] / cnt.EV2RCM / cnt.AU2EV)
    };

    const E = &[_]T{9.42, 10.11, 10.48, 11.08};

    const omega   = &[_]T{omega_all[0], omega_all[1], omega_all[2], omega_all[3], omega_all[4],  0.00000,  0.00000,  0.00000};

    const kappa_0 = &[_]T{     0.00000,      0.00000,     -0.02203,     -0.12147,     -0.09468,  0.00000,  0.00000,  0.00000};
    const kappa_1 = &[_]T{     0.00000,      0.00000,      0.09074,      0.05316,      0.04454,  0.00000,  0.00000,  0.00000};
    const kappa_2 = &[_]T{     0.00000,      0.00000,      0.02748,      0.11233,      0.14539,  0.00000,  0.00000,  0.00000};
    const kappa_3 = &[_]T{     0.00000,      0.00000,     -0.04054,      0.00747,      0.00050,  0.00000,  0.00000,  0.00000};

    const gamma_0 = &[_]T{     0.00000,      0.00000,      0.01938,      0.01489,      0.00970,  0.00000,  0.00000,  0.00000};
    const gamma_1 = &[_]T{     0.00000,      0.00000,      0.00694,      0.00828,      0.00096,  0.00000,  0.00000,  0.00000};
    const gamma_2 = &[_]T{     0.00000,      0.00000,     -0.00294,      0.00183,     -0.00114,  0.00000,  0.00000,  0.00000};
    const gamma_3 = &[_]T{     0.00000,      0.00000,      0.00752,      0.00546,      0.01108,  0.00000,  0.00000,  0.00000};

    const d0_0    = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000, 41.89704,  4.80270, 22.92802};
    const d0_1    = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000, 38.37122, 74.15995, 18.27440};
    const d0_2    = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000, 39.25691, 90.76928,  9.46894};
    const d0_3    = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000, 37.97847, 20.56079, 65.09678};

    const a_0     = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000,  0.04719, -0.13675,  0.07438};
    const a_1     = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000,  0.05231, -0.03064,  0.07911};
    const a_2     = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000,  0.05286, -0.03374,  0.08653};
    const a_3     = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000,  0.05431, -0.08044,  0.03660};

    const q0_0    = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000,  0.81440,  0.02883, -0.32069};
    const q0_1    = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000,  0.37488, -1.34468, -0.01711};
    const q0_2    = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000,  0.14859, -0.29923,  0.37635};
    const q0_3    = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000, -0.18152,  0.38841,  1.66312};

    const e0_0    = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000, -0.06431, -0.00007, -0.01274};
    const e0_1    = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000, -0.01505, -0.12082, -0.00003};
    const e0_2    = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000, -0.00244, -0.00916, -0.01037};
    const e0_3    = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000, -0.00366, -0.02071, -0.25639};

    const k_0     = &[_]T{     0.03317,      0.02979,      0.00000,      0.00000,      0.00000,  0.00000,  0.00000,  0.00000};
    const k_1     = &[_]T{     0.01157,      0.01488,      0.00000,      0.00000,      0.00000,  0.00000,  0.00000,  0.00000};
    const k_2     = &[_]T{     0.01534,      0.01671,      0.00000,      0.00000,      0.00000,  0.00000,  0.00000,  0.00000};
    const k_3     = &[_]T{     0.00000,      0.00000,      0.00000,      0.00000,      0.00000,  0.00000,  0.00000,  0.00000};

    const l_01    = &[_]T{     0.04633,      0.03540,      0.00000,      0.00000,      0.00000,  0.00000,  0.00000,  0.00000};
    const l_12    = &[_]T{     0.03148,      0.03607,      0.00000,      0.00000,      0.00000,  0.00000,  0.00000,  0.00000};
    const l_02    = &[_]T{     0.00000,      0.00000,     -0.03538,     -0.02049,      0.00000,  0.00000,  0.00114,  0.13035};
    const l_13    = &[_]T{     0.00000,      0.00000,      0.08077,      0.00000,      0.07284, -0.08132,  0.12606,  0.14272};

    U.ptr(0, 0).* += E[0]; U.ptr(1, 1).* += E[1]; U.ptr(2, 2).* += E[2]; U.ptr(3, 3).* += E[3];

    for (0..omega.len) |i| {

        U.ptr(0, 0).* += 0.5 * omega[i] * q[i] * q[i] / cnt.EV2RCM;
        U.ptr(1, 1).* += 0.5 * omega[i] * q[i] * q[i] / cnt.EV2RCM;
        U.ptr(2, 2).* += 0.5 * omega[i] * q[i] * q[i] / cnt.EV2RCM;
        U.ptr(3, 3).* += 0.5 * omega[i] * q[i] * q[i] / cnt.EV2RCM;

        U.ptr(0, 0).* += d0_0[i] * std.math.pow(T, std.math.exp(-a_0[i] * (q[i] - q0_0[i])) - 1, 2) + e0_0[i];
        U.ptr(1, 1).* += d0_1[i] * std.math.pow(T, std.math.exp(-a_1[i] * (q[i] - q0_1[i])) - 1, 2) + e0_1[i];
        U.ptr(2, 2).* += d0_2[i] * std.math.pow(T, std.math.exp(-a_2[i] * (q[i] - q0_2[i])) - 1, 2) + e0_2[i];
        U.ptr(3, 3).* += d0_3[i] * std.math.pow(T, std.math.exp(-a_3[i] * (q[i] - q0_3[i])) - 1, 2) + e0_3[i];

        U.ptr(0, 0).* += kappa_0[i] * q[i];
        U.ptr(1, 1).* += kappa_1[i] * q[i];
        U.ptr(2, 2).* += kappa_2[i] * q[i];
        U.ptr(3, 3).* += kappa_3[i] * q[i];

        U.ptr(0, 0).* += 0.5 * gamma_0[i] * q[i] * q[i];
        U.ptr(1, 1).* += 0.5 * gamma_1[i] * q[i] * q[i];
        U.ptr(2, 2).* += 0.5 * gamma_2[i] * q[i] * q[i];
        U.ptr(3, 3).* += 0.5 * gamma_3[i] * q[i] * q[i];

        U.ptr(0, 0).* += k_0[i] * q[i] * q[i] * q[i] * q[i] / 24;
        U.ptr(1, 1).* += k_1[i] * q[i] * q[i] * q[i] * q[i] / 24;
        U.ptr(2, 2).* += k_2[i] * q[i] * q[i] * q[i] * q[i] / 24;
        U.ptr(3, 3).* += k_3[i] * q[i] * q[i] * q[i] * q[i] / 24;

        U.ptr(0, 1).* += l_01[i] * q[i];
        U.ptr(1, 2).* += l_12[i] * q[i];
        U.ptr(0, 2).* += l_02[i] * q[i];
        U.ptr(1, 3).* += l_13[i] * q[i];
    }

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(2, 1).* = U.at(1, 2);
    U.ptr(2, 0).* = U.at(0, 2);
    U.ptr(3, 1).* = U.at(1, 3);

    for (0..4) |i| for (0..4) |j| {U.ptr(i, j).* /= cnt.AU2EV;};
}

/// The 12-dimensional uracil model potential energy surface.
pub fn uracil12D_1(comptime T: type, U: *Matrix(T), r: Vector(T)) !void {
    U.fill(0);

    const E = &[_]T{9.42, 10.11, 10.48, 11.08};

    const omega   = &[_]T{ 388.00000, 560.00000, 734.00000, 770.00000, 771.00000, 1193.00000, 1228.00000, 1383.00000, 1406.00000,  0.00000,  0.00000,  0.00000};

    const kappa_0 = &[_]T{   0.04139,  -0.05367,   0.00000,   0.05357,   0.00000,   -0.02203,    0.07472,   -0.12147,   -0.09468,  0.00000,  0.00000,  0.00000};
    const kappa_1 = &[_]T{  -0.02688,  -0.00503,   0.00000,  -0.00503,   0.00000,    0.09074,   -0.00582,    0.05316,    0.04454,  0.00000,  0.00000,  0.00000};
    const kappa_2 = &[_]T{  -0.05853,   0.00775,   0.00000,   0.03697,   0.00000,    0.02748,    0.00889,    0.11233,    0.14539,  0.00000,  0.00000,  0.00000};
    const kappa_3 = &[_]T{   0.00132,  -0.04581,   0.00000,  -0.02130,   0.00000,   -0.04054,   -0.02136,    0.00747,    0.00050,  0.00000,  0.00000,  0.00000};

    const gamma_0 = &[_]T{   0.00383,   0.00000,   0.00000,   0.00000,   0.00000,    0.01938,    0.01590,    0.01489,    0.00970,  0.00000,  0.00000,  0.00000};
    const gamma_1 = &[_]T{  -0.00426,   0.00000,   0.00000,   0.00000,   0.00000,    0.00694,    0.01348,    0.00828,    0.00096,  0.00000,  0.00000,  0.00000};
    const gamma_2 = &[_]T{  -0.00366,   0.00000,   0.00000,   0.00000,   0.00000,   -0.00294,    0.00901,    0.00183,   -0.00114,  0.00000,  0.00000,  0.00000};
    const gamma_3 = &[_]T{  -0.00947,   0.00000,   0.00000,   0.00000,   0.00000,    0.00752,    0.00497,    0.00546,    0.01108,  0.00000,  0.00000,  0.00000};

    const d0_0    = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000, 41.89704,  4.80270, 22.92802};
    const d0_1    = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000, 38.37122, 74.15995, 18.27440};
    const d0_2    = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000, 39.25691, 90.76928,  9.46894};
    const d0_3    = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000, 37.97847, 20.56079, 65.09678};

    const a_0     = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.04719, -0.13675,  0.07438};
    const a_1     = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.05231, -0.03064,  0.07911};
    const a_2     = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.05286, -0.03374,  0.08653};
    const a_3     = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.05431, -0.08044,  0.03660};

    const q0_0    = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.81440,  0.02883, -0.32069};
    const q0_1    = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.37488, -1.34468, -0.01711};
    const q0_2    = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.14859, -0.29923,  0.37635};
    const q0_3    = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.18152,  0.38841,  1.66312};

    const e0_0    = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.06431, -0.00007, -0.01274};
    const e0_1    = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.01505, -0.12082, -0.00003};
    const e0_2    = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.00244, -0.00916, -0.01037};
    const e0_3    = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000, -0.00366, -0.02071, -0.25639};

    const k_0     = &[_]T{   0.00000,   0.00000,   0.03317,   0.00000,   0.02979,    0.00000,    0.00000,    0.00000,    0.00000,  0.00000,  0.00000,  0.00000};
    const k_1     = &[_]T{   0.00000,   0.00000,   0.01157,   0.00000,   0.01488,    0.00000,    0.00000,    0.00000,    0.00000,  0.00000,  0.00000,  0.00000};
    const k_2     = &[_]T{   0.00000,   0.00000,   0.01534,   0.00000,   0.01671,    0.00000,    0.00000,    0.00000,    0.00000,  0.00000,  0.00000,  0.00000};
    const k_3     = &[_]T{   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  0.00000,  0.00000,  0.00000};

    const l_01    = &[_]T{   0.00000,   0.00000,   0.04633,   0.00000,   0.03540,    0.00000,    0.00000,    0.00000,    0.00000,  0.00000,  0.00000,  0.00000};
    const l_12    = &[_]T{   0.00000,   0.00000,   0.03148,   0.00000,   0.03607,    0.00000,    0.00000,    0.00000,    0.00000,  0.00000,  0.00000,  0.00000};
    const l_02    = &[_]T{   0.00000,   0.02141,   0.00000,   0.02082,   0.00000,   -0.03538,   -0.03763,   -0.02049,    0.00000,  0.00000,  0.00114,  0.13035};
    const l_13    = &[_]T{   0.00000,   0.00000,   0.00000,   0.01845,   0.00000,    0.08077,    0.05834,    0.00000,    0.07284, -0.08132,  0.12606,  0.14272};

    U.ptr(0, 0).* += E[0]; U.ptr(1, 1).* += E[1]; U.ptr(2, 2).* += E[2]; U.ptr(3, 3).* += E[3];

    for (0..12) |i| {

        U.ptr(0, 0).* += 0.5 * omega[i] * r.at(i) * r.at(i) / cnt.EV2RCM;
        U.ptr(1, 1).* += 0.5 * omega[i] * r.at(i) * r.at(i) / cnt.EV2RCM;
        U.ptr(2, 2).* += 0.5 * omega[i] * r.at(i) * r.at(i) / cnt.EV2RCM;
        U.ptr(3, 3).* += 0.5 * omega[i] * r.at(i) * r.at(i) / cnt.EV2RCM;

        U.ptr(0, 0).* += d0_0[i] * std.math.pow(T, std.math.exp(-a_0[i] * (r.at(i) - q0_0[i])) - 1, 2) + e0_0[i];
        U.ptr(1, 1).* += d0_1[i] * std.math.pow(T, std.math.exp(-a_1[i] * (r.at(i) - q0_1[i])) - 1, 2) + e0_1[i];
        U.ptr(2, 2).* += d0_2[i] * std.math.pow(T, std.math.exp(-a_2[i] * (r.at(i) - q0_2[i])) - 1, 2) + e0_2[i];
        U.ptr(3, 3).* += d0_3[i] * std.math.pow(T, std.math.exp(-a_3[i] * (r.at(i) - q0_3[i])) - 1, 2) + e0_3[i];

        U.ptr(0, 0).* += kappa_0[i] * r.at(i);
        U.ptr(1, 1).* += kappa_1[i] * r.at(i);
        U.ptr(2, 2).* += kappa_2[i] * r.at(i);
        U.ptr(3, 3).* += kappa_3[i] * r.at(i);

        U.ptr(0, 0).* += 0.5 * gamma_0[i] * r.at(i) * r.at(i);
        U.ptr(1, 1).* += 0.5 * gamma_1[i] * r.at(i) * r.at(i);
        U.ptr(2, 2).* += 0.5 * gamma_2[i] * r.at(i) * r.at(i);
        U.ptr(3, 3).* += 0.5 * gamma_3[i] * r.at(i) * r.at(i);

        U.ptr(0, 0).* += k_0[i] * r.at(i) * r.at(i) * r.at(i) * r.at(i) / 24;
        U.ptr(1, 1).* += k_1[i] * r.at(i) * r.at(i) * r.at(i) * r.at(i) / 24;
        U.ptr(2, 2).* += k_2[i] * r.at(i) * r.at(i) * r.at(i) * r.at(i) / 24;
        U.ptr(3, 3).* += k_3[i] * r.at(i) * r.at(i) * r.at(i) * r.at(i) / 24;

        U.ptr(0, 1).* += l_01[i] * r.at(i);
        U.ptr(1, 2).* += l_12[i] * r.at(i);
        U.ptr(0, 2).* += l_02[i] * r.at(i);
        U.ptr(1, 3).* += l_13[i] * r.at(i);
    }

    U.ptr(1, 0).* = U.at(0, 1);
    U.ptr(2, 1).* = U.at(1, 2);
    U.ptr(2, 0).* = U.at(0, 2);
    U.ptr(3, 1).* = U.at(1, 3);

    for (0..4) |i| for (0..4) |j| {U.ptr(i, j).* /= cnt.AU2EV;};
}

/// Generate a grid in the k-space. The result is stored in the matrix k. Both the start and end values are included.
pub fn kgrid(comptime T: type, k: *Matrix(T), start: T, end: T, points: u32) void {
    k.fill(2 * std.math.pi / asfloat(T, points) / (end - start) * asfloat(T, points - 1));

    for (0..k.rows / points) |i| for (0..points) |j| {
        k.ptr(i * points + j, k.cols - 1).* *= asfloat(T, j) - asfloat(T, if (j < points / 2) 0 else points);
    };

    for (0..k.rows) |i| for (0..k.cols) |j| {
        k.ptr(i, j).* = k.at(i / std.math.pow(usize, points, k.cols - j - 1), k.cols - 1);
    };
}

/// Generate a grid in the r-space. The result is stored in the matrix r. Both the start and end values are included.
pub fn rgrid(comptime T: type, r: *Matrix(T), start: T, end: T, points: u32) void {
    for (0..r.rows) |i| for (0..r.cols) |j| {
        r.ptr(i, j).* = start + asfloat(T, i / std.math.pow(usize, points, r.cols - j - 1) % points) * (end - start) / asfloat(T, points - 1);
    };
}

/// Write the potential energy surfaces to a file.
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
