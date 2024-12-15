const std = @import("std");

const hfm = @import("hartreefock.zig");
const mat = @import("matrix.zig"     );
const ten = @import("tensor.zig"     );
const tns = @import("transform.zig"  );

const Matrix = @import("matrix.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

pub fn ConfigurationInteractionOptions(comptime T: type) type {
    return struct {
        excitation: ?[]const u32 = null,

        hartree_fock: hfm.HartreeFockOptions(T) = .{},
    };
}

pub fn ConfigurationInteractionOutput(comptime T: type) type {
    return struct {
        E: T,

        pub fn deinit(self: ConfigurationInteractionOutput(T)) void {
            _ = self;
        }
    };
}

pub fn run(comptime T: type, opt: ConfigurationInteractionOptions(T), print: bool, allocator: std.mem.Allocator) !ConfigurationInteractionOutput(T) {
    const hf = try hfm.run(T, opt.hartree_fock, print, allocator); const nbf = 2 * hf.nbf; const nocc = 2 * hf.nocc; var ndet: usize = 1;

    var J_MS_A = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_MS_A.deinit();

    var H_MS = try Matrix(T).init(nbf, nbf, allocator); defer H_MS.deinit();

    {
        const T_AO = try mat.read(T, opt.hartree_fock.integral.kinetic,    allocator); defer T_AO.deinit();
        const V_AO = try mat.read(T, opt.hartree_fock.integral.nuclear,    allocator); defer V_AO.deinit();
        const J_AO = try ten.read(T, opt.hartree_fock.integral.coulomb, 4, allocator); defer J_AO.deinit();

        var H_AO = try Matrix(T).init(T_AO.rows, T_AO.cols, allocator); defer H_AO.deinit(); mat.add(T, &H_AO, T_AO, V_AO);

        var J_AS = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_AS.deinit();
        var J_MS = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_MS.deinit();

        var H_AS = try Matrix(T).init(nbf, nbf, allocator); defer H_AS.deinit();
        var C_MS = try Matrix(T).init(nbf, nbf, allocator); defer C_MS.deinit();

        tns.cfsMO2MS(T, &C_MS, hf.C_MO);

        tns.oneAO2AS(T, &H_AS, H_AO); tns.oneAO2MO(T, &H_MS,  H_AS, C_MS);
        tns.twoAO2AS(T, &J_AS, J_AO); tns.twoAO2MO(T, &J_MS, &J_AS, C_MS);

        for (0..J_MS.shape[0]) |i| for (0..J_MS.shape[1]) |j| for (0..J_MS.shape[2]) |k| for (0..J_MS.shape[3]) |l| {
            J_MS_A.ptr(&[_]usize{i, k, j, l}).* = J_MS.at(&[_]usize{i, j, k, l}) - J_MS.at(&[_]usize{i, l, k, j});
        };

    }

    for (nocc + 1..nbf + 1) |i| ndet *= i;
    for (2..nbf - nocc + 1) |i| ndet /= i;

    try std.io.getStdOut().writer().print("\nNUMBER OF CI DETERMINANTS: {}\n", .{ndet});

    var T1 = try Matrix(T).init(ndet, ndet, allocator); defer T1.deinit();
    var T2 = try Matrix(T).init(ndet, ndet, allocator); defer T1.deinit();

    var A = try Vector(usize).init(nocc,       allocator); defer A.deinit();
    var D = try Matrix(usize).init(ndet, nocc, allocator); defer D.deinit();
    var H = try Matrix(T    ).init(ndet, ndet, allocator); defer H.deinit();
    var E = try Matrix(T    ).init(ndet, ndet, allocator); defer E.deinit();
    var C = try Matrix(T    ).init(ndet, ndet, allocator); defer C.deinit();

    for (0..D.cols) |i| D.ptr(0, i).* = i;

    for (1..D.rows) |i| {

        @memcpy(D.row(i).data, D.row(i - 1).data); var index: usize = undefined;

        for (0..D.cols) |j| if (D.at(i, D.cols - j - 1) != nbf - j - 1) {
            D.ptr(i, D.cols - j - 1).* += 1; index = D.cols - j - 1; break;
        };

        for (index + 1..D.cols) |j| D.ptr(i, j).* = D.at(i, j - 1) + 1;
    }

    H.fill(0);

    for (0..ndet) |i| for (i..ndet) |j| {

        @memcpy(A.data, D.row(j).data); var k: usize = 0; var so: [4]usize = undefined; var diff: u32 = 0; var sign: i32 = 1;

        while (k < nocc) : (k += 1) if (A.at(k) != D.at(i, k)) for (k + 1..nocc) |l| if (A.at(k) == D.at(i, l)) {

            for (k..l) |m| {
                std.mem.swap(usize, A.ptr(m), A.ptr(m + 1)); sign *= -1;
            }

            k -= 1; break;
        };

        for (0..nocc) |l| if (A.at(l) != D.at(i, l)) {diff += 1;}; var up: usize = 0;

        for (0..nocc) |l| if (A.at(l) != D.at(i, l)) {

            if (diff == 1) {
                so[0] = A.at(l); so[1] = D.at(i, l);
            } 

            if (diff == 2) {
                so[if (up == 0) 0 else 1] = A.at(l); so[if (up == 0) 2 else 3] = D.at(i, l); up += 1;
            } 
        };

        if (diff == 0) {

            for (0..nocc) |l| {
                H.ptr(i, j).* += H_MS.at(A.at(l), A.at(l));
            }

            for (0..nocc) |l| for (0..nocc) |m| {
                H.ptr(i, j).* += 0.5 * J_MS_A.at(&[_]usize{A.at(l), A.at(m), A.at(l), A.at(m)});
            };
        }

        if (diff == 1) {

            H.ptr(i, j).* += H_MS.at(so[0], so[1]);

            for (0..nocc) |m| if (A.at(m) != so[0] and D.at(i, m) != so[1]) {
                H.ptr(i, j).* += J_MS_A.at(&[_]usize{so[0], A.at(m), so[1], A.at(m)});
            };
        }

        if (diff == 2) {
            H.ptr(i, j).* = J_MS_A.at(&[_]usize{so[0], so[1], so[2], so[3]});
        }

        // if (i == 0 and j == 35) {
        //     std.debug.print("{}/{} - {d} {d} - {d:20.14}\n{any}\n{any}\n{any}\n{any}\n\n", .{
        //         i, j, diff, sign, asfloat(T, sign) * H.at(i, j), D.row(i).data, D.row(j).data, A.data, so
        //     });
        // }
        // if (i == 2 and j == 18) {
        //     std.debug.print("{}/{} - {d} {d} - {d:20.14}\n{any}\n{any}\n{any}\n{any}\n\n", .{
        //         i, j, diff, sign, asfloat(T, sign) * H.at(i, j), D.row(i).data, D.row(j).data, A.data, so
        //     });
        // }

        H.ptr(i, j).* = asfloat(T, sign) * H.at(i, j); H.ptr(j, i).* = H.at(i, j);
    };

    // try H.write("HCI_Z.mat");

    mat.eigh(T, &E, &C, H, &T1, &T2);

    std.debug.print("\nCI ENERGY: {d:.14}\n", .{E.at(0, 0) + hf.VNN});

    return .{
        .E = E.at(0, 0) + hf.VNN,
    };
}
