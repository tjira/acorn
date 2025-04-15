//! Module to perform the Hartree-Fock calculation.

const std = @import("std");

const A2AU  = @import("constant.zig").A2AU ;
const SM2AN = @import("constant.zig").SM2AN;

const bls = @import("blas.zig"    );
const inp = @import("input.zig"   );
const int = @import("integral.zig");
const lpk = @import("lapack.zig"  );
const lbt = @import("libint.zig"  );
const mat = @import("matrix.zig"  );
const mth = @import("math.zig"    );
const out = @import("output.zig"  );
const ten = @import("tensor.zig"  );
const vec = @import("vector.zig"  );

const ContractedGaussian = @import("contractedgaussian.zig").ContractedGaussian;
const Basis              = @import("basis.zig" ).Basis                         ;
const Matrix             = @import("matrix.zig").Matrix                        ;
const System             = @import("system.zig").System                        ;
const Tensor             = @import("tensor.zig").Tensor                        ;
const Vector             = @import("vector.zig").Vector                        ;

const asfloat = @import("helper.zig").asfloat;
const uncr    = @import("helper.zig").uncr   ;

/// Run the Hartree-Fock calculation with the given options.
pub fn run(comptime T: type, opt: inp.HartreeFockOptions(T), print: bool, allocator: std.mem.Allocator) !out.HartreeFockOutput(T) {
    if (opt.integral.basis == null and (opt.integral.overlap == null or opt.integral.kinetic == null or opt.integral.nuclear == null or opt.integral.coulomb == null)) return error.MissingIntegral;

    if (opt.system_file == null and (opt.system.coords == null or opt.system.atoms == null)) return error.SystemNotFullySpecified;

    var basis: std.ArrayList(ContractedGaussian(T)) = undefined; var system: System(T) = undefined;

    if (opt.system_file == null and opt.system.coords != null) {

        system = System(T){
            .coords = try Matrix(T).init(opt.system.coords.?.len, 3, allocator),
            .atoms  = try Vector(T).init(opt.system.atoms .?.len,    allocator),
            .nocc   = (mth.sum(u8, opt.system.atoms.?) - opt.system.charge) / 2,
        };

        for (0..opt.system.atoms.?.len) |i| {
            system.atoms.ptr(i).* = @as(T, @floatFromInt(opt.system.atoms.?[i]));
        }

        for (0..opt.system.coords.?.len) |i| for (0..3) |j| {
            system.coords.ptr(i, j).* = opt.system.coords.?[i][j] * A2AU;
        };
    }

    if (opt.system_file != null   ) system = try System(T).read(opt.system_file.?,            allocator);
    if (opt.integral.basis != null) basis  = try Basis(T).get  (system, opt.integral.basis.?, allocator);

    const nbf = basis.items.len; var npg: usize = 0; const nocc = system.nocc; const VNN = system.nuclearRepulsion();

    for (0..nbf) |i| npg += basis.items[i].c.len;

    if (print) try std.io.getStdOut().writer().print("\n# OF BASIS FUNCTIONS:      {d}\n", .{nbf});
    if (print) try std.io.getStdOut().writer().print(  "# OF PRIMITIVE GAUSSIANS:  {d}\n", .{npg});

    var timer = try std.time.Timer.start();

    const S_AO = if (opt.integral.overlap != null) try mat.read(T, opt.integral.overlap.?,    allocator) else try int.overlap(T, basis,         allocator);
    const T_AO = if (opt.integral.kinetic != null) try mat.read(T, opt.integral.kinetic.?,    allocator) else try int.kinetic(T, basis,         allocator);
    const V_AO = if (opt.integral.nuclear != null) try mat.read(T, opt.integral.nuclear.?,    allocator) else try int.nuclear(T, basis, system, allocator);
    const J_AO = if (opt.integral.coulomb != null) try ten.read(T, opt.integral.coulomb.?, 4, allocator) else try int.coulomb(T, basis,         allocator);

    if (opt.integral.libint) lbt.overlap(S_AO.data, "molecule.xyz", "sto-3g");
    if (opt.integral.libint) lbt.kinetic(T_AO.data, "molecule.xyz", "sto-3g");
    if (opt.integral.libint) lbt.nuclear(V_AO.data, "molecule.xyz", "sto-3g");
    if (opt.integral.libint) lbt.coulomb(J_AO.data, "molecule.xyz", "sto-3g");

    if (print) try std.io.getStdOut().writer().print("\nINTEGRALS OBTAINED: {}\n", .{std.fmt.fmtDuration(timer.read())});

    var T1 = try Matrix(T).init(nbf, nbf, allocator); defer T1.deinit();
    var T2 = try Matrix(T).init(nbf, nbf, allocator); defer T2.deinit();
    var T3 = try Matrix(T).init(nbf, nbf, allocator); defer T3.deinit();

    var J_AO_A = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_AO_A.deinit();

    var H_AO = try Matrix(T).init(nbf, nbf, allocator); defer H_AO.deinit();
    var ERR  = try Matrix(T).init(nbf, nbf, allocator); defer  ERR.deinit();

    var DIIS_E = std.ArrayList(Matrix(T)).init(allocator); defer DIIS_E.deinit();
    var DIIS_F = std.ArrayList(Matrix(T)).init(allocator); defer DIIS_F.deinit();

    if (opt.dsize != null) for (0..opt.dsize.?) |_| {
        try DIIS_E.append(try Matrix(T).init(H_AO.rows, H_AO.cols, allocator));
        try DIIS_F.append(try Matrix(T).init(H_AO.rows, H_AO.cols, allocator));
    };

    {
        var K_AO = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer K_AO.deinit();

        ten.transpose(T, &K_AO, J_AO, &[_]usize{0, 3, 2, 1}); ten.muls(T, &K_AO, K_AO, 0.5); ten.sub(T, &J_AO_A, J_AO, K_AO);

        mat.add(T, &H_AO, T_AO, V_AO);
    }

    var F_AO = try Matrix(T).init(nbf, nbf, allocator); F_AO.fill(0);
    var C_MO = try Matrix(T).init(nbf, nbf, allocator); C_MO.fill(0);
    var D_MO = try Matrix(T).init(nbf, nbf, allocator); D_MO.fill(0);
    var E_MO = try Matrix(T).init(nbf, nbf, allocator); E_MO.fill(0);

    var iter: u32 = 0; var E: T = 0; var EP: T = 1;

    if (print) try std.io.getStdOut().writer().print("\n{s:4} {s:20} {s:8}\n", .{"ITER", "ENERGY", "|DELTA E|"});

    while (@abs(EP - E) > opt.threshold) : (iter += 1) {

        if (iter == opt.maxiter) return error.MaxIterationsExceeded;

        @memcpy(F_AO.data, H_AO.data);

        for (0..nbf) |i| for (0..nbf) |j| for (0..nbf) |k| for (0..nbf) |l| {
            F_AO.ptr(k, l).* += D_MO.at(i, j) * J_AO_A.at(&[_]usize{i, j, k, l});
        };

        if (opt.dsize != null) {

            bls.dgemm(&T1, S_AO, false, D_MO, false); bls.dgemm(&T2, T1, false, F_AO, true); bls.dgemm(&T1, F_AO, false, D_MO, false); bls.dgemm(&T3, T1, false, S_AO, true); mat.sub(T, &ERR, T2, T3);

            @memcpy(DIIS_F.items[@intCast(@mod(@as(i32, @intCast(iter)) - 1, @as(i32, @intCast(DIIS_F.items.len))))].data, F_AO.data);
            @memcpy(DIIS_E.items[@intCast(@mod(@as(i32, @intCast(iter)) - 1, @as(i32, @intCast(DIIS_E.items.len))))].data,  ERR.data);

            try diisExtrapolate(T, &F_AO, &DIIS_F, &DIIS_E, iter, allocator);
        }

        lpk.dsygvd(&E_MO, &C_MO, F_AO, S_AO, &T1); D_MO.fill(0); EP = E; E = 0;

        for (0..nbf) |i| for (0..nocc) |j| for (0..nbf) |k| {
            D_MO.ptr(i, k).* += 2.0 * C_MO.at(i, j) * C_MO.at(k, j);
        };

        for (0..nbf) |i| for (0..nbf) |j| {
            E += 0.5 * D_MO.at(i, j) * (H_AO.at(i, j) + F_AO.at(i, j));
        };

        if (print) try std.io.getStdOut().writer().print("{d:4} {d:20.14} {e:9.3}\n", .{iter + 1, E + VNN, @abs(EP - E)});
    }

    if (print) try std.io.getStdOut().writer().print("\nHF ENERGY: {d:.14}\n", .{E + VNN});

    for (0..DIIS_E.items.len) |i| DIIS_E.items[i].deinit();
    for (0..DIIS_F.items.len) |i| DIIS_F.items[i].deinit();

    if (opt.integral.basis != null) basis.deinit();

    if (opt.write.coefficient_ao != null) try C_MO.write(opt.write.coefficient_ao.?);
    if (opt.write.coulomb_ao     != null) try J_AO.write(opt.write.coulomb_ao.?    );
    if (opt.write.density_ao     != null) try D_MO.write(opt.write.density_ao.?    );
    if (opt.write.kinetic_ao     != null) try T_AO.write(opt.write.kinetic_ao.?    );
    if (opt.write.nuclear_ao     != null) try V_AO.write(opt.write.nuclear_ao.?    );
    if (opt.write.overlap_ao     != null) try S_AO.write(opt.write.overlap_ao.?    );

    return out.HartreeFockOutput(T){
        .S_AO = S_AO, .T_AO = T_AO, .V_AO = V_AO, .J_AO = J_AO, .C_MO = C_MO, .D_MO = D_MO, .E_MO = E_MO, .F_AO = F_AO, .E = E + VNN, .system = system
    };
}

/// Extrapolate the DIIS error to obtain a new Fock matrix.
pub fn diisExtrapolate(comptime T: type, F_AO: *Matrix(T), DIIS_F: *std.ArrayList(Matrix(T)), DIIS_E: *std.ArrayList(Matrix(T)), iter: u32, allocator: std.mem.Allocator) !void {
    if (iter > 0) {

        const size = if (iter < DIIS_F.items.len) iter else DIIS_F.items.len;

        var A   = try Matrix(T  ).init(size + 1, size + 1, allocator); defer   A.deinit();
        var ALU = try Matrix(T  ).init(size + 1, size + 1, allocator); defer ALU.deinit();
        var b   = try Vector(T  ).init(size + 1,           allocator); defer   b.deinit();
        var c   = try Vector(T  ).init(size + 1,           allocator); defer   c.deinit();
        var p   = try Vector(i32).init(size + 1,           allocator); defer   p.deinit();

        A.fill(1); b.fill(0); A.ptr(A.rows - 1, A.cols - 1).* = 0; b.ptr(b.rows - 1).* = 1;

        for (0..size) |i| for (0..size) |j| {

            A.ptr(i, j).* = 0;

            const ii = @mod(@as(i32, @intCast(iter)) - @as(i32, @intCast(size)) + @as(i32, @intCast(i)), @as(i32, @intCast(DIIS_E.items.len)));
            const jj = @mod(@as(i32, @intCast(iter)) - @as(i32, @intCast(size)) + @as(i32, @intCast(j)), @as(i32, @intCast(DIIS_E.items.len)));

            for (0..DIIS_E.items[0].rows) |k| for (0..DIIS_E.items[0].cols) |l| {
                A.ptr(i, j).* += DIIS_E.items[@intCast(ii)].at(k, l) * DIIS_E.items[@intCast(jj)].at(k, l);
            };
        };

        lpk.dgesv(&c, &ALU, &p, A, b); F_AO.fill(0); 

        for (0..size) |i| {

            const ii = @mod(@as(i32, @intCast(iter)) - @as(i32, @intCast(size)) + @as(i32, @intCast(i)), @as(i32, @intCast(DIIS_F.items.len)));

            for (0..F_AO.rows) |j| for (0..F_AO.cols) |k| {
                F_AO.ptr(j, k).* += c.at(i) * DIIS_F.items[@intCast(ii)].at(j, k);
            };
        }
    }
}
