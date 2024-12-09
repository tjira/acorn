const std = @import("std");

const hfm = @import("hartreefock.zig");
const ten = @import("tensor.zig"     );

const Matrix = @import("matrix.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

pub fn MollerPlessetOptions(comptime T: type) type {
    return struct {
        order: u32 = 2,

        hartree_fock: hfm.HartreeFockOptions(T) = .{},
    };
}

pub fn MollerPlessetOutput(comptime T: type) type {
    return struct {
        E: T,

        pub fn deinit(self: MollerPlessetOutput(T)) void {
            _ = self;
        }
    };
}

pub fn cfsMO2MS(comptime T: type, C_MS: *Matrix(T), C_MO: Matrix(T)) void {
    C_MS.fill(0);

    for (0..C_MO.rows) |i| for (0..C_MO.cols) |j| {
        C_MS.ptr(i            , 2 * j    ).* = C_MO.at(i, j);
        C_MS.ptr(i + C_MO.rows, 2 * j + 1).* = C_MO.at(i, j);
    };
}

pub fn oneAO2AS(comptime T: type, A_AS: *Matrix(T), A_AO: Matrix(T)) void {
    A_AS.fill(0);

    for (0..A_AO.rows) |i| for (0..A_AO.cols) |j| {
        A_AS.ptr(i            , j            ).* = A_AO.at(i, j);
        A_AS.ptr(i + A_AO.rows, j + A_AO.cols).* = A_AO.at(i, j);
    };
}

pub fn oneAO2MO(comptime T: type, A_MO: *Matrix(T), A_AO: Matrix(T), C_MO: Matrix(T)) void {
    A_MO.fill(0);

    for (0..A_AO.rows) |i| for (0..A_AO.cols) |j| for (0..A_MO.rows) |p| for (0..A_MO.cols) |q| {
        A_MO.ptr(p, q).* += C_MO.at(i, p) * A_AO.at(i, j) * C_MO.at(j, q);
    };
}

pub fn twoAO2AS(comptime T: type, A_AS: *Tensor(T), A_AO: Tensor(T)) void {
    A_AS.fill(0);

    for (0..A_AO.shape[0]) |i| for (0..A_AO.shape[1]) |j| for (0..A_AO.shape[2]) |k| for (0..A_AO.shape[3]) |l| {
        A_AS.ptr(&[_]usize{i,                 j,                 k                , l                }).* = A_AO.at(&[_]usize{i, j, k, l});
        A_AS.ptr(&[_]usize{i,                 j                , k + A_AO.shape[2], l + A_AO.shape[3]}).* = A_AO.at(&[_]usize{i, j, k, l});
        A_AS.ptr(&[_]usize{i + A_AO.shape[0], j + A_AO.shape[1], k,                 l                }).* = A_AO.at(&[_]usize{i, j, k, l});
        A_AS.ptr(&[_]usize{i + A_AO.shape[0], j + A_AO.shape[1], k + A_AO.shape[2], l + A_AO.shape[3]}).* = A_AO.at(&[_]usize{i, j, k, l});
    };
}

pub fn twoAO2MO(comptime T: type, A_MO: *Tensor(T), A_AO: Tensor(T), C_MO: Matrix(T)) void {
    A_MO.fill(0);

    for (0..A_AO.shape[0]) |i| for (0..A_AO.shape[1]) |j| for (0..A_AO.shape[2]) |k| for (0..A_AO.shape[3]) |l| {
        for (0..A_MO.shape[0]) |p| for (0..A_MO.shape[1]) |q| for (0..A_MO.shape[2]) |r| for (0..A_MO.shape[3]) |s| {
            A_MO.ptr(&[_]usize{p, q, r, s}).* += C_MO.at(i, p) * C_MO.at(j, q) * A_AO.at(&[_]usize{i, j, k, l}) * C_MO.at(k, r) * C_MO.at(l, s);
        };
    };
}

pub fn run(comptime T: type, opt: MollerPlessetOptions(T), print: bool, allocator: std.mem.Allocator) !MollerPlessetOutput(T) {
    if (opt.order != 2) return error.PertrubationOrderNotImplemented;

    const hf = try hfm.run(T, opt.hartree_fock, print, allocator); const nbf = 2 * hf.nbf; const nocc = 2 * hf.nocc;

    var J_MS_A = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_MS_A.deinit();
    
    {
        const J_AO = try ten.read(T, opt.hartree_fock.integral.coulomb, 4, allocator); defer J_AO.deinit();

        var J_AS = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_AS.deinit();
        var J_MS = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_MS.deinit();

        var C_MS = try Matrix(T).init(nbf, nbf, allocator); defer C_MS.deinit();

        cfsMO2MS(T, &C_MS, hf.C_MO); twoAO2AS(T, &J_AS, J_AO); twoAO2MO(T, &J_MS, J_AS, C_MS);

        for (0..J_MS.shape[0]) |i| for (0..J_MS.shape[1]) |j| for (0..J_MS.shape[2]) |k| for (0..J_MS.shape[3]) |l| {
            J_MS_A.ptr(&[_]usize{i, k, j, l}).* = J_MS.at(&[_]usize{i, j, k, l}) - J_MS.at(&[_]usize{i, l, k, j});
        };

    }

    var E: T = 0;

    for (0..nocc) |i| for (0..nocc) |j| for (nocc..nbf) |a| for (nocc..nbf) |b| {
        E += 0.25 * J_MS_A.at(&[_]usize{i, j, a, b}) * J_MS_A.at(&[_]usize{a, b, i, j}) / (hf.E_MO.at(i / 2, i / 2) + hf.E_MO.at(j / 2, j / 2) - hf.E_MO.at(a / 2, a / 2) - hf.E_MO.at(b / 2, b / 2));
    };

    if (print) try std.io.getStdOut().writer().print("\nMP{d} ENERGY: {d:.14}\n", .{opt.order, hf.E + E});

    return .{
        .E = hf.E + E,
    };
}
