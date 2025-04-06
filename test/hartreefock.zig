const std = @import("std");

const hartree_fock = @import("acorn").hartree_fock;
const input        = @import("acorn").input;

const expect = @import("main.zig").expect;
const log    = @import("main.zig").log   ;

const allocator = std.testing.allocator;

test "hartree_fock_hydrogen" {
    const opt_sto3g = input.HartreeFockOptions(f64){
        .integral = .{
            .basis = "sto-3g"
        },
        .system = .{
            .atoms = &[_]u8{1, 1},
            .coords = &[_][3]f64{
                .{0.14, 0, 0},
                .{0.86, 0, 0},
            },
        },
    };

    var opt_631g   = opt_sto3g;   opt_631g.integral.basis = "6-31g"  ;
    var opt_ccpvdz = opt_sto3g; opt_ccpvdz.integral.basis = "cc-pvdz";
    var opt_ccpvtz = opt_sto3g; opt_ccpvtz.integral.basis = "cc-pvtz";

    const output_sto3g  = try hartree_fock.run(f64, opt_sto3g,  false, allocator); defer  output_sto3g.deinit();
    const output_631g   = try hartree_fock.run(f64, opt_631g,   false, allocator); defer   output_631g.deinit();
    const output_ccpvdz = try hartree_fock.run(f64, opt_ccpvdz, false, allocator); defer output_ccpvdz.deinit();
    const output_ccpvtz = try hartree_fock.run(f64, opt_ccpvtz, false, allocator); defer output_ccpvtz.deinit();

    try expect(output_sto3g.E,  -1.11744507742767);
    try expect(output_631g.E,   -1.12675331625770);
    try expect(output_ccpvdz.E, -1.12815326004132);
    try expect(output_ccpvtz.E, -1.13285933203976);
}
