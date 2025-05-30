//! Module for Configuration Interaction (CI) calculations.

const std = @import("std");

const hfm = @import("hartreefock.zig");
const inp = @import("input.zig"      );
const lpk = @import("lapack.zig"     );
const mat = @import("matrix.zig"     );
const mth = @import("math.zig"       );
const out = @import("output.zig"     );
const ten = @import("tensor.zig"     );
const tns = @import("transform.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

/// The main function to run the CI calculations.
pub fn run(comptime T: type, opt: inp.ConfigurationInteractionOptions(T), print: bool, allocator: std.mem.Allocator) !out.ConfigurationInteractionOutput(T) {
    const hf = try hfm.run(T, opt.hartree_fock, print, allocator);

    const D = try generateDeterminants(hf.S_AO.rows, hf.system.nocc, allocator); defer D.deinit();

    const nbf = 2 * hf.S_AO.rows; const nocc = 2 * hf.system.nocc;

    var J_MS_A = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_MS_A.deinit();

    var H_MS = try Matrix(T).init(nbf, nbf, allocator); defer H_MS.deinit();

    try transform(T, &H_MS, &J_MS_A, hf.T_AO, hf.V_AO, hf.J_AO, hf.C_MO, allocator);

    if (print) try std.io.getStdOut().writer().print("\nNUMBER OF CI DETERMINANTS: {d}\n", .{D.rows});

    var A = try Vector(usize).init(nocc,           allocator); defer A.deinit();
    var H = try Matrix(T    ).init(D.rows, D.rows, allocator); defer H.deinit();
    var E = try Matrix(T    ).init(D.rows, D.rows, allocator); defer E.deinit();
    var C = try Matrix(T    ).init(D.rows, D.rows, allocator); defer C.deinit();

    H.fill(0);

    for (0..D.rows) |i| for (i..D.rows) |j| {

        var so: [4]usize = undefined; var diff: u32 = 0; var k: usize = 0;

        const sign = try alignDeterminant(&A, D.row(i).vector(), D.row(j).vector());

        for (0..nocc) |l| if (A.at(l) != D.at(i, l)) {diff += 1;};

        if (diff > 2) continue;

        if (diff == 1) for (0..nocc) |l| if (A.at(l) != D.at(i, l)) {
            so[0] = A.at(l); so[1] = D.at(i, l);
        };

        if (diff == 2) for (0..nocc) |l| if (A.at(l) != D.at(i, l)) {
            so[if (k == 0) 0 else 1] = A.at(l); so[if (k == 0) 2 else 3] = D.at(i, l); k += 1;
        };

        H.ptr(i, j).* = asfloat(T, sign) * try slater(T, A, so[0..diff * 2], H_MS, J_MS_A); H.ptr(j, i).* = H.at(i, j);
    };

    lpk.dsyevd(&E, &C, H);

    if (print) try std.io.getStdOut().writer().print("\nCI ENERGY: {d:.14}\n", .{E.at(0, 0) + hf.system.nuclearRepulsion()});

    return .{
        .E = E.at(0, 0) + hf.system.nuclearRepulsion(),
    };
}

/// Aligns the vector C to the vector B. The result is stored in the vector A and the sign of the permutation is returned.
fn alignDeterminant(A: *Vector(usize), B: Vector(usize), C: Vector(usize)) !i32 {
    @memcpy(A.data, C.data); var k: i32 = 0; var sign: i32 = 1;

    while (k < A.rows) : (k += 1) {
        if (A.at(@intCast(k)) != B.at(@intCast(k))) for (@as(usize, @intCast(k)) + 1..A.rows) |l| if (A.at(@intCast(k)) == B.at(l) or A.at(l) == B.at(@intCast(k))) {
            std.mem.swap(usize, A.ptr(@intCast(k)), A.ptr(l)); sign *= -1; k -= 1; break;
        };
    }

    return sign;
}

/// Generates the determinants for the CI calculations.
fn generateDeterminants(nbf: usize, nocc: usize, allocator: std.mem.Allocator) !Matrix(usize) {
    const data_alpha = try allocator.alloc(usize, nbf); defer allocator.free(data_alpha);
    const data_beta  = try allocator.alloc(usize, nbf); defer allocator.free(data_beta );
    
    for (0..data_alpha.len) |i| data_alpha[i] = 2 * i + 0;
    for (0..data_beta.len ) |i| data_beta[i]  = 2 * i + 1;

    const dets_alpha = try mth.combinations(usize, data_alpha, nocc, allocator); defer dets_alpha.deinit();
    const dets_beta  = try mth.combinations(usize, data_beta,  nocc, allocator); defer  dets_beta.deinit();

    const D = try Matrix(usize).init(dets_alpha.items.len * dets_beta.items.len, 2 * nocc, allocator);

    for (0..dets_alpha.items.len) |i| for (0..dets_beta.items.len) |j| for (0..nocc) |k| {
        D.ptr(i * dets_beta.items.len + j, k       ).* = dets_alpha.items[i][k];
        D.ptr(i * dets_beta.items.len + j, k + nocc).* =  dets_beta.items[j][k];
    };

    return D;
}

/// Slater-Condon rules for the CI calculations.
fn slater(comptime T: type, A: Vector(usize), so: []const usize, H_MS: Matrix(T), J_MS_A: Tensor(T)) !T {
    var hij: T = 0;

    if (so.len / 2 == 0) {

        for (0..A.rows) |l| {
            hij += H_MS.at(A.at(l), A.at(l));
        }

        for (0..A.rows) |l| for (0..A.rows) |m| {
            hij += 0.5 * J_MS_A.at(&[_]usize{A.at(l), A.at(m), A.at(l), A.at(m)});
        };
    }

    else if (so.len / 2 == 1) {

        hij += H_MS.at(so[0], so[1]);

        for (0..A.rows) |m| if (A.at(m) != so[0]) {
            hij += J_MS_A.at(&[_]usize{so[0], A.at(m), so[1], A.at(m)});
        };
    }

    else if (so.len / 2 == 2) {
        hij = J_MS_A.at(&[_]usize{so[0], so[1], so[2], so[3]});
    }
    
    return hij;
}

/// Function to perform all integrals transformations used in the CI calculations.
fn transform(comptime T: type, H_MS: *Matrix(T), J_MS_A: *Tensor(T), T_AO: Matrix(T), V_AO: Matrix(T), J_AO: Tensor(T), C_MO: Matrix(T), allocator: std.mem.Allocator) !void {
    var H_AO = try Matrix(T).init(T_AO.rows, T_AO.cols, allocator); defer H_AO.deinit(); mat.add(T, &H_AO, T_AO, V_AO);

    var J_AS = try Tensor(T).init(&[_]usize{J_MS_A.shape[0], J_MS_A.shape[0], J_MS_A.shape[0], J_MS_A.shape[0]}, allocator); defer J_AS.deinit();
    var J_MS = try Tensor(T).init(&[_]usize{J_MS_A.shape[0], J_MS_A.shape[0], J_MS_A.shape[0], J_MS_A.shape[0]}, allocator); defer J_MS.deinit();

    var H_AS = try Matrix(T).init(J_MS_A.shape[0], J_MS_A.shape[0], allocator); defer H_AS.deinit();
    var C_MS = try Matrix(T).init(J_MS_A.shape[0], J_MS_A.shape[0], allocator); defer C_MS.deinit();

    tns.cfsMO2MS(T, &C_MS, C_MO); tns.oneAO2AS(T, &H_AS, H_AO); tns.oneAO2MO(T,  H_MS,  H_AS, C_MS); tns.twoAO2AS(T, &J_AS, J_AO); try tns.twoAO2MO(T, &J_MS, &J_AS, C_MS);

    for (0..J_MS.shape[0]) |i| for (0..J_MS.shape[1]) |j| for (0..J_MS.shape[2]) |k| for (0..J_MS.shape[3]) |l| {
        J_MS_A.ptr(&[_]usize{i, k, j, l}).* = J_MS.at(&[_]usize{i, j, k, l}) - J_MS.at(&[_]usize{i, l, k, j});
    };
}
