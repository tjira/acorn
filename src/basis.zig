//! Basis library.

const std = @import("std");

const ContractedGaussian = @import("contractedgaussian.zig").ContractedGaussian;
const System             = @import("system.zig"            ).System            ;

/// Basis struct.
pub fn Basis(comptime T: type) type {
    return struct {

        /// Get the basis set for the given system and name.
        pub fn get(system: System(T), name: []const u8, allocator: std.mem.Allocator) !std.ArrayList(ContractedGaussian(T)) {
            var basis = std.ArrayList(ContractedGaussian(T)).init(allocator); var parsed: std.json.Parsed(std.json.Value) = undefined; defer parsed.deinit();

            const lower = try allocator.alloc(u8, name.len); defer allocator.free(lower); _ = std.ascii.lowerString(lower, name);

            if      (std.mem.eql(u8, lower, "sto-3g")) {parsed = try std.json.parseFromSlice(std.json.Value, allocator, @embedFile("basis/sto-3g.json"), .{});}
            else if (std.mem.eql(u8, lower, "6-31g" )) {parsed = try std.json.parseFromSlice(std.json.Value, allocator, @embedFile("basis/6-31g.json" ), .{});}

            else return error.BasisNameNotFound;

            for (0..system.atoms.rows) |i| {

                var an: [3]u8 = undefined; const ans = std.fmt.formatIntBuf(&an, @as(u8, @intFromFloat(system.atoms.at(i))), 10, std.fmt.Case.lower, .{});

                if (parsed.value.object.get("elements").?.object.get(an[0..ans]) == null) return error.AtomNotFoundInBasis;

                const shells = parsed.value.object.get("elements").?.object.get(an[0..ans]).?.object.get("electron_shells").?.array.items;

                for (shells) |shell| {

                    const am        = shell.object.get("angular_momentum").?.array.items[0].integer    ;
                    const exponents = shell.object.get("exponents"       ).?.array.items               ;
                    const coefs     = shell.object.get("coefficients"    ).?.array.items[0].array.items;

                    var a     = try allocator.alloc(T, 3            ); defer     allocator.free(a);
                    var c     = try allocator.alloc(T, coefs.len    ); defer     allocator.free(c);
                    var alpha = try allocator.alloc(T, exponents.len); defer allocator.free(alpha);

                    for (0..@as(usize, @intCast(2 * am + 1))) |j| {

                        a[0] = 0; a[1] = 0; a[2] = 0; if (am == 1) {a[j] = @as(T, @floatFromInt(am));} else if (am > 1) return error.InvalidAngularMomentum;

                        for (0..c.len    ) |k| c[k]     = try std.fmt.parseFloat(T,     coefs[k].string);
                        for (0..alpha.len) |k| alpha[k] = try std.fmt.parseFloat(T, exponents[k].string);

                        try basis.append(try ContractedGaussian(T).init(system.getCoords(i), .{a[0], a[1], a[2]}, c, alpha, allocator));
                    }
                }
            }

            return basis;
        }
    };
}
