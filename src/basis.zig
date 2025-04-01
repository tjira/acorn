//! Basis library.

const std = @import("std");

const ContractedGaussian = @import("contractedgaussian.zig").ContractedGaussian;
const System             = @import("system.zig"            ).System            ;

/// Basis struct.
pub fn Basis(comptime T: type) type {
    return struct {

        /// Get the basis set for the given system and name.
        pub fn get(system: System(T), name: []const u8, allocator: std.mem.Allocator) !std.ArrayList(ContractedGaussian(T)) {

            if (std.mem.eql(u8, name, "STO-3G")) return STO_3G(system, allocator);

            return error.InvalidBasisName;
        }

        /// Get the STO-3G basis set for the given system.
        pub fn STO_3G(system: System(T), allocator: std.mem.Allocator) !std.ArrayList(ContractedGaussian(T)) {
            var basis = std.ArrayList(ContractedGaussian(T)).init(allocator);

            for (0..system.atoms.rows) |i| {

                var a     = std.ArrayList([3]      T).init(allocator); defer     a.deinit();
                var c     = std.ArrayList([ ]const T).init(allocator); defer     c.deinit();
                var alpha = std.ArrayList([ ]const T).init(allocator); defer alpha.deinit();

                switch (@as(u32, @intFromFloat(system.atoms.at(i)))) {
                    1  => {
                        try     a.append( [3]T{0.0,          0.0,          0.0         });
                        try     c.append(&[_]T{0.1543289673, 0.5353281423, 0.4446345422});
                        try alpha.append(&[_]T{3.4252509140, 0.6239137298, 0.1688554040});
                    },
                    else => {return error.InvalidAtomicNumber;}
                }

                for (0..a.items.len) |j| try basis.append(try ContractedGaussian(T).init(system.getCoords(i), a.items[j], c.items[j], alpha.items[j], allocator));
            }

            return basis;
        }
    };
}
