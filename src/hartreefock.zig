const std = @import("std");

const A2AU  = @import("constant.zig").A2AU ;
const SM2AN = @import("constant.zig").SM2AN;

const mat = @import("matrix.zig");
const ten = @import("tensor.zig");
const vec = @import("vector.zig");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;
const Tensor = @import("tensor.zig").Tensor;

const asfloat = @import("helper.zig").asfloat;
const uncr    = @import("helper.zig").uncr   ;

pub fn HartreeFockOptions(comptime T: type) type {
    return struct {
        const DirectInversion = struct {
            size: u32 = 5,
            start: u32 = 100
        };
        const Integral = struct {
            overlap: []const u8 = "S_AO.mat",
            kinetic: []const u8 = "T_AO.mat",
            nuclear: []const u8 = "V_AO.mat",
            coulomb: []const u8 = "J_AO.mat"
        };

        molecule: []const u8 = "molecule.xyz",
        threshold: T = 1e-12,
        maxiter: u32 = 100,

        diis: DirectInversion = .{}, integral: Integral = .{}
    };
}

pub fn HartreeFockOutput(comptime T: type) type {
    return struct {
        C_MO: Matrix(T), D_MO: Matrix(T), E_MO: Matrix(T), F_AO: Matrix(T), E: T, VNN: T, nbf: usize, nocc: usize,

        pub fn deinit(self: HartreeFockOutput(T)) void {
            self.C_MO.deinit(); self.D_MO.deinit(); self.E_MO.deinit(); self.F_AO.deinit();
        }
    };
}

pub fn run(comptime T: type, opt: HartreeFockOptions(T), print: bool, allocator: std.mem.Allocator) !HartreeFockOutput(T) {
    const S_AO = try mat.read(T, opt.integral.overlap, allocator); defer S_AO.deinit();

    const system = try parseSystem(T, opt.molecule, print, allocator);

    const nbf = S_AO.cols; const nocc = system.nocc; const VNN = system.VNN;

    var T1 = try Matrix(T).init(nbf, nbf, allocator); defer T1.deinit();
    var T2 = try Matrix(T).init(nbf, nbf, allocator); defer T2.deinit();
    var T3 = try Matrix(T).init(nbf, nbf, allocator); defer T3.deinit();
    var T4 = try Matrix(T).init(nbf, nbf, allocator); defer T4.deinit();

    var J_AO_A = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_AO_A.deinit();

    var H_AO = try Matrix(T).init(nbf, nbf, allocator); defer H_AO.deinit();
    var X    = try Matrix(T).init(nbf, nbf, allocator); defer    X.deinit();

    {
        const T_AO = try mat.read(T, opt.integral.kinetic,    allocator); defer T_AO.deinit();
        const V_AO = try mat.read(T, opt.integral.nuclear,    allocator); defer V_AO.deinit();
        const J_AO = try ten.read(T, opt.integral.coulomb, 4, allocator); defer J_AO.deinit();

        var XJ = try Matrix(T).init(nbf, nbf, allocator); defer XJ.deinit();
        var XC = try Matrix(T).init(nbf, nbf, allocator); defer XC.deinit();

        var K_AO = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer K_AO.deinit();

        ten.transpose(T, &K_AO, J_AO, &[_]usize{0, 3, 2, 1}); ten.muls(T, &K_AO, K_AO, 0.5); ten.sub(T, &J_AO_A, J_AO, K_AO);

        mat.eigh(T, &XJ, &XC, S_AO, &T1, &T2); for (0..nbf) |i| XJ.ptr(i, i).* = 1.0 / std.math.sqrt(XJ.at(i, i));

        mat.mm(T, &T1, XC, XJ); mat.transpose(T, &T2, XC); mat.mm(T, &X, T1, T2);

        mat.add(T, &H_AO, T_AO, V_AO);
    }

    var F_AO = try Matrix(T).init(nbf, nbf, allocator); F_AO.fill(0);
    var C_MO = try Matrix(T).init(nbf, nbf, allocator); C_MO.fill(0);
    var D_MO = try Matrix(T).init(nbf, nbf, allocator); D_MO.fill(0);
    var E_MO = try Matrix(T).init(nbf, nbf, allocator); E_MO.fill(0);

    var DIIS_E = std.ArrayList(Matrix(T)).init(allocator); defer DIIS_E.deinit();
    var DIIS_F = std.ArrayList(Matrix(T)).init(allocator); defer DIIS_F.deinit();

    var DIIS_B = try Matrix(T).init(opt.diis.size + 1, opt.diis.size + 1, allocator); defer DIIS_B.deinit();
    var DIIS_b = try Vector(T).init(opt.diis.size + 1,                    allocator); defer DIIS_b.deinit();
    var DIIS_c = try Vector(T).init(opt.diis.size + 1,                    allocator); defer DIIS_c.deinit();

    for (0..opt.diis.size) |_| {
        try DIIS_E.append(try Matrix(T).init(nbf, nbf, allocator));
        try DIIS_F.append(try Matrix(T).init(nbf, nbf, allocator));
    }

    var iter: u32 = 0; var E: T = 0; var EP: T = 1;

    if (print) try std.io.getStdOut().writer().print("\n{s:4} {s:20} {s:8}\n", .{"ITER", "ENERGY", "|DELTA E|"});

    while (@abs(EP - E) > opt.threshold) : (iter += 1) {

        if (iter == opt.maxiter) return error.MaxIterationsExceeded;

        @memcpy(F_AO.data, H_AO.data);

        for (0..nbf) |i| for (0..nbf) |j| for (0..nbf) |k| for (0..nbf) |l| {
            F_AO.ptr(k, l).* += D_MO.at(i, j) * J_AO_A.at(&[_]usize{i, j, k, l});
        };

        if (iter + 1 >= opt.diis.start) {

            mat.mm(T, &T1, S_AO, D_MO); mat.mm(T, &T2, T1, F_AO); mat.mm(T, &T1, F_AO, D_MO); mat.mm(T, &T3, T1, S_AO); mat.sub(T, &T1, T2, T3);

            @memcpy(DIIS_F.items[iter % opt.diis.size].data, F_AO.data); @memcpy(DIIS_E.items[iter % opt.diis.size].data, T1.data);

            DIIS_B.fill(1); DIIS_b.fill(0); DIIS_B.ptr(opt.diis.size, opt.diis.size).* = 0; DIIS_b.ptr(opt.diis.size).* = 1;

            for (0..opt.diis.size) |i| for (0..opt.diis.size) |j| for (0.. DIIS_E.items[i].data.len) |k| {
                DIIS_B.ptr(i, j).* += DIIS_E.items[(iter - opt.diis.size + i + 1) % opt.diis.size].data[k] * DIIS_E.items[(iter - opt.diis.size + j + 1) % opt.diis.size].data[k];
            };

            mat.linsolve(T, &DIIS_c, DIIS_B, DIIS_b);
        }

        mat.mm(T, &T2, X, F_AO); mat.mm(T, &T1, T2, X); mat.eigh(T, &E_MO, &T2, T1, &T3, &T4); mat.mm(T, &C_MO, X, T2);

        D_MO.fill(0); EP = E; E = 0;

        for (0..nbf) |i| for (0..nocc) |j| for (0..nbf) |k| {
            D_MO.ptr(i, k).* += 2.0 * C_MO.at(i, j) * C_MO.at(k, j);
        };

        for (0..nbf) |i| for (0..nbf) |j| {
            E += 0.5 * D_MO.at(i, j) * (H_AO.at(i, j) + F_AO.at(i, j));
        };

        if (print) try std.io.getStdOut().writer().print("{d:4} {d:20.14} {e:9.3}\n", .{iter + 1, E + VNN, @abs(EP - E)});
    }

    for (0..opt.diis.size) |i| {
        DIIS_E.items[i].deinit(); DIIS_F.items[i].deinit();
    }

    if (print) try std.io.getStdOut().writer().print("\nHF ENERGY: {d:.14}\n", .{E + VNN});

    return HartreeFockOutput(T){
        .C_MO = C_MO, .D_MO = D_MO, .E_MO = E_MO, .F_AO = F_AO, .E = E + VNN, .VNN = VNN, .nbf = nbf, .nocc = nocc
    };
}

pub fn parseSystem(comptime T: type, path: []const u8, print: bool, allocator: std.mem.Allocator) !struct {natom: u32, nocc: u32, VNN: T} {
    const file = try std.fs.cwd().openFile(path, .{}); defer file.close(); var buffer: [64]u8 = undefined;

    var buffered = std.io.bufferedReader(file.reader()); var reader = buffered.reader();
    var stream   = std.io.fixedBufferStream(&buffer);  const writer =   stream.writer();

    stream.reset(); try reader.streamUntilDelimiter(writer, '\n', 64);

    const natom = try std.fmt.parseInt(u32, uncr(stream.getWritten()), 10);

    stream.reset(); try reader.streamUntilDelimiter(writer, '\n', 64);

    var coords = try Matrix(T).init(natom, 3, allocator); defer coords.deinit(); coords.fill(0);
    var atoms  = try Vector(T).init(natom,    allocator); defer  atoms.deinit();  atoms.fill(0);

    for (0..natom) |i| {

        stream.reset(); try reader.streamUntilDelimiter(writer, '\n', 64);

        var it = std.mem.splitScalar(u8, uncr(stream.getWritten()), ' '); 

        while (it.next()) |token| if (token.len > 0 and atoms.at(i) == 0) {
            atoms.ptr(i).* = asfloat(T, SM2AN.get(token).?); break;
        };

        var j: i32 = 0;

        while (it.next()) |token| : (j += if (token.len > 0) 1 else 0) if (token.len > 0) {
            coords.ptr(i, @as(usize, @intCast(j))).* = try std.fmt.parseFloat(T, token);
        };
    }

    var VNN: T = 0; var nocc: u32 = 0; for (0..natom) |i| nocc +=  @intFromFloat(atoms.at(i));

    for (0..natom) |i| for (0..natom) |j| if (i != j) {

        const x1 = coords.at(i, 0); const y1 = coords.at(i, 1); const z1 = coords.at(i, 2);
        const x2 = coords.at(j, 0); const y2 = coords.at(j, 1); const z2 = coords.at(j, 2);

        const r = std.math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)) * A2AU;

        VNN += 0.5 * atoms.at(i) * atoms.at(j) / r;
    };

    if (print) try std.io.getStdOut().writer().print("\nMOLECULE: {s}\n", .{path});

    if (print) for (0..natom) |i| {
        try std.io.getStdOut().writer().print("{d:2} {d:14.8} {d:14.8} {d:14.8}\n", .{@as(u32, @intFromFloat(atoms.at(i))), coords.at(i, 0), coords.at(i, 1), coords.at(i, 2)});
    };

    return .{.natom = natom, .nocc = nocc / 2, .VNN = VNN};
}
