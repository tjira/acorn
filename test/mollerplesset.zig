const std = @import("std");

const input          = @import("acorn").input;
const moller_plesset = @import("acorn").moller_plesset;

const expect = @import("main.zig").expect;
const log    = @import("main.zig").log   ;

const allocator = std.testing.allocator;

test "moller_plesset_hydrogen_2" {
    const opt_sto3g = input.MollerPlessetOptions(f64){
        .order = 2,

        .hartree_fock = input.HartreeFockOptions(f64){
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
        }
    };

    var opt_631g   = opt_sto3g;   opt_631g.hartree_fock.integral.basis = "6-31g"  ;
    var opt_ccpvdz = opt_sto3g; opt_ccpvdz.hartree_fock.integral.basis = "cc-pvdz";
    var opt_ccpvtz = opt_sto3g; opt_ccpvtz.hartree_fock.integral.basis = "cc-pvtz";

    const output_sto3g  = try moller_plesset.run(f64, opt_sto3g,  false, allocator); defer  output_sto3g.deinit();
    const output_631g   = try moller_plesset.run(f64, opt_631g,   false, allocator); defer   output_631g.deinit();
    const output_ccpvdz = try moller_plesset.run(f64, opt_ccpvdz, false, allocator); defer output_ccpvdz.deinit();
    const output_ccpvtz = try moller_plesset.run(f64, opt_ccpvtz, false, allocator); defer output_ccpvtz.deinit();

    try expect(output_sto3g.E,  -1.13012364534534);
    try expect(output_631g.E,   -1.14391755479002);
    try expect(output_ccpvdz.E, -1.15434251439068);
    try expect(output_ccpvtz.E, -1.16458063441218);
}
