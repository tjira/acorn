//! Basis library.

const std = @import("std");

const System = @import("system.zig").System;

/// Basis struct.
pub fn Basis(comptime T: type) type {
    return struct {

        /// Return the basis set as an array of consecutive numbers.
        pub fn array(system: System(T), name: []const u8, allocator: std.mem.Allocator) ![]const T {
            var basis = std.ArrayList(T).init(allocator); const parsed = try std.json.parseFromSlice(std.json.Value, allocator, try embedded(name, "json", allocator), .{});

            for (0..system.atoms.rows) |i| {

                var an: [3]u8 = undefined; const ans = std.fmt.formatIntBuf(&an, @as(u8, @intFromFloat(system.atoms.at(i))), 10, std.fmt.Case.lower, .{});

                if (parsed.value.object.get("elements").?.object.get(an[0..ans]) == null) return error.AtomNotFoundInBasis;

                const shells = parsed.value.object.get("elements").?.object.get(an[0..ans]).?.object.get("electron_shells").?.array.items;

                for (shells) |shell| {

                    const am        = shell.object.get("angular_momentum").?.array.items[0].integer;
                    const exponents = shell.object.get("exponents"       ).?.array.items           ;
                    const coefs     = shell.object.get("coefficients"    ).?.array.items           ;

                    var c     = try allocator.alloc(T, exponents.len); defer     allocator.free(c);
                    var alpha = try allocator.alloc(T, exponents.len); defer allocator.free(alpha);

                    for (0..alpha.len) |j| alpha[j] = try std.fmt.parseFloat(T, exponents[j].string);

                    for (coefs) |coef| {

                        try basis.append(@as(T, @floatFromInt(c.len)));
                        try basis.append(@as(T, @floatFromInt(am   )));

                        for (0..3) |j| try basis.append(system.getCoords(i)[j]);

                        for (0..c.len) |j| c[j] = try std.fmt.parseFloat(T, coef.array.items[j].string);

                        for (0..alpha.len) |j| try basis.append(alpha[j]);
                        for (0..c.len    ) |j| try basis.append(    c[j]);
                    }
                }
            }

            parsed.deinit(); return basis.items;
        }

        /// Returns the basis set for the given name and format.
        pub fn embedded(basis: []const u8, format: []const u8, allocator: std.mem.Allocator) ![]const u8 {
            const lower = try allocator.alloc(u8, basis.len); defer allocator.free(lower); _ = std.ascii.lowerString(lower, basis);

            if      (std.mem.eql(u8, lower, "sto-3g" ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/sto-3g.json" );}
            else if (std.mem.eql(u8, lower, "6-31g"  ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-31g.json"  );}
            else if (std.mem.eql(u8, lower, "cc-pvdz") and std.mem.eql(u8, format, "json")) {return @embedFile("basis/cc-pvdz.json");}
            else if (std.mem.eql(u8, lower, "cc-pvtz") and std.mem.eql(u8, format, "json")) {return @embedFile("basis/cc-pvtz.json");}

            else return error.BasisNameNotFound;
        }
    };
}
