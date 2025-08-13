const std = @import("std");

const configuration_interaction = @import("acorn").configuration_interaction;
const input                     = @import("acorn").input;

const expect = @import("acorn.zig").expect;
const log    = @import("acorn.zig").log   ;

const allocator = std.testing.allocator;

test "configuration_interaction_hydrogen_full" {
    const opt_sto3g = input.ConfigurationInteractionOptions(f64){
        .active_space = null,

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

    var opt_631g = opt_sto3g; opt_631g.hartree_fock.integral.basis = "6-31g";

    const output_sto3g = try configuration_interaction.run(f64, opt_sto3g, false, allocator); defer output_sto3g.deinit();
    const output_631g  = try configuration_interaction.run(f64, opt_631g,  false, allocator); defer  output_631g.deinit();

    try expect(output_sto3g.E, -1.13711171524389);
    try expect(output_631g.E,  -1.15122475612240);
}
