//! Basis library.

const std = @import("std");

const BSE   = @import("bse.zig"     ).BSE  ;
const SM2AN = @import("constant.zig").SM2AN;

const ContractedGaussian = @import("contractedgaussian.zig").ContractedGaussian;
const System             = @import("system.zig"            ).System            ;

/// Basis struct.
pub fn Basis(comptime T: type) type {
    return struct {

        /// Get the basis set for the given system and name.
        pub fn get(system: System(T), name: []const u8, allocator: std.mem.Allocator) !std.ArrayList(ContractedGaussian(T)) {
            var basis = std.ArrayList(ContractedGaussian(T)).init(allocator);

            const lower = try allocator.alloc(u8, name.len); defer allocator.free(lower); _ = std.ascii.lowerString(lower, name);

            if (BSE.get(lower) == null) return error.BasisNotFound;

            for (0..system.atoms.rows) |i| {

                var a     = std.ArrayList([3]      T).init(allocator); defer     a.deinit();
                var c     = std.ArrayList([ ]const T).init(allocator); defer     c.deinit();
                var alpha = std.ArrayList([ ]const T).init(allocator); defer alpha.deinit();

                for (SM2AN.keys()) |symbol| if (SM2AN.get(symbol).? == @as(u32, @intFromFloat(system.atoms.at(i)))) {

                    if (BSE.get(lower).?.get(symbol) == null) return error.AtomNotFoundInBasis;

                    const array = BSE.get(lower).?.get(symbol).?; const ncgs = @as(usize, @intFromFloat(array[0]));

                    for (0..ncgs) |j| {

                        const npgs = @as(usize, @intFromFloat(array[1 + j])); var li = 1 + ncgs;

                        for (0..j) |k| li += 2 * @as(usize, @intFromFloat(array[1 + k])) + 1;

                        const l = @as(usize, @intFromFloat(array[li]));

                        for (0..2 * l + 1) |k| {

                            var ak: [3]T = .{0, 0, 0};

                            if (l == 1) {ak[k] = 1;} else if (l > 1) return error.InvalidAngularMomentum;

                            try a.append(ak); try c.append(array[li + 1..li + npgs + 1]); try alpha.append(array[li + npgs + 1..li + 2 * npgs + 1]);
                        }
                    }
                };

                for (0..a.items.len) |j| try basis.append(try ContractedGaussian(T).init(system.getCoords(i), a.items[j], c.items[j], alpha.items[j], allocator));
            }

            return basis;
        }
    };
}
