//! Basis library.

const std = @import("std");

const ContractedGaussian = @import("contractedgaussian.zig").ContractedGaussian;
const System             = @import("system.zig"            ).System            ;

/// Basis struct.
pub fn Basis(comptime T: type) type {
    return struct {

        /// Get the basis set for the given system and name.
        pub fn get(system: System(T), name: []const u8, allocator: std.mem.Allocator) !std.ArrayList(ContractedGaussian(T)) {

            const lower = try allocator.alloc(u8, name.len); defer allocator.free(lower);

            _ = std.ascii.lowerString(lower, name);

            if (std.mem.eql(u8, lower, "sto-3g")) return STO_3G(system, allocator);

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
                        try     a.append( [3]T{ 0.0,              0.0,              0.0            });
                        try     c.append(&[_]T{ 0.1543289673e+0,  0.5353281423e+0,  0.4446345422e+0});
                        try alpha.append(&[_]T{ 0.3425250914e+1,  0.6239137298e+0,  0.1688554040e+0});
                    },
                    6  => {
                        try     a.append( [3]T{ 0.0,              0.0,              0.0            });
                        try     c.append(&[_]T{ 0.1543289673e+0,  0.5353281423e+0,  0.4446345422e+0});
                        try alpha.append(&[_]T{ 0.7161683735e+2,  0.1304509632e+2,  0.3530512160e+1});

                        try     a.append( [3]T{ 0.0,              0.0,              0.0            });
                        try     c.append(&[_]T{-0.9996722919e-1,  0.3995128261e+0,  0.7001154689e+0});
                        try alpha.append(&[_]T{ 0.2941249355e+1,  0.6834830964e+0,  0.2222899159e+0});

                        try     a.append( [3]T{ 1.0,              0.0,              0.0            });
                        try     c.append(&[_]T{ 0.1559162750e+0,  0.6076837186e+0,  0.3919573931e+0});
                        try alpha.append(&[_]T{ 0.2941249355e+1,  0.6834830964e+0,  0.2222899159e+0});

                        try     a.append( [3]T{ 0.0,              1.0,              0.0            });
                        try     c.append(&[_]T{ 0.1559162750e+0,  0.6076837186e+0,  0.3919573931e+0});
                        try alpha.append(&[_]T{ 0.2941249355e+1,  0.6834830964e+0,  0.2222899159e+0});

                        try     a.append( [3]T{ 0.0,              0.0,              1.0            });
                        try     c.append(&[_]T{ 0.1559162750e+0,  0.6076837186e+0,  0.3919573931e+0});
                        try alpha.append(&[_]T{ 0.2941249355e+1,  0.6834830964e+0,  0.2222899159e+0});
                    },
                    8  => {
                        try     a.append( [3]T{ 0.0,              0.0,              0.0            });
                        try     c.append(&[_]T{ 0.1543289673e+0,  0.5353281423e+0,  0.4446345422e+0});
                        try alpha.append(&[_]T{ 0.1307093214e+3,  0.2380886605e+2,  0.6443608313e+1});

                        try     a.append( [3]T{ 0.0,              0.0,              0.0            });
                        try     c.append(&[_]T{-0.9996722919e-1,  0.3995128261e+0,  0.7001154689e+0});
                        try alpha.append(&[_]T{ 0.5033151319e+1,  0.1169596125e+1,  0.3803889600e+0});

                        try     a.append( [3]T{ 1.0,              0.0,              0.0            });
                        try     c.append(&[_]T{ 0.1559162750e+0,  0.6076837186e+0,  0.3919573931e+0});
                        try alpha.append(&[_]T{ 0.5033151319e+1,  0.1169596125e+1,  0.3803889600e+0});

                        try     a.append( [3]T{ 0.0,              1.0,              0.0            });
                        try     c.append(&[_]T{ 0.1559162750e+0,  0.6076837186e+0,  0.3919573931e+0});
                        try alpha.append(&[_]T{ 0.5033151319e+1,  0.1169596125e+1,  0.3803889600e+0});

                        try     a.append( [3]T{ 0.0,              0.0,              1.0            });
                        try     c.append(&[_]T{ 0.1559162750e+0,  0.6076837186e+0,  0.3919573931e+0});
                        try alpha.append(&[_]T{ 0.5033151319e+1,  0.1169596125e+1,  0.3803889600e+0});
                    },
                    else => {return error.InvalidAtomicNumber;}
                }

                for (0..a.items.len) |j| try basis.append(try ContractedGaussian(T).init(system.getCoords(i), a.items[j], c.items[j], alpha.items[j], allocator));
            }

            return basis;
        }
    };
}
