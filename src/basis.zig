//! Basis library.

const std = @import("std");

const ContractedGaussian = @import("contractedgaussian.zig").ContractedGaussian;
const System             = @import("system.zig"            ).System            ;

/// Basis struct.
pub fn Basis(comptime T: type) type {
    return struct {

        pub fn embedded(basis: []const u8, format: []const u8, allocator: std.mem.Allocator) ![]const u8 {
            const lower = try allocator.alloc(u8, basis.len); defer allocator.free(lower); _ = std.ascii.lowerString(lower, basis);

            if      (std.mem.eql(u8, lower, "sto-3g" ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/sto-3g.json" );}
            else if (std.mem.eql(u8, lower, "sto-3g" ) and std.mem.eql(u8, format, "g94" )) {return @embedFile("basis/sto-3g.g94"  );}
            else if (std.mem.eql(u8, lower, "6-31g"  ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-31g.json"  );}
            else if (std.mem.eql(u8, lower, "6-31g"  ) and std.mem.eql(u8, format, "g94" )) {return @embedFile("basis/6-31g.g94"   );}
            else if (std.mem.eql(u8, lower, "cc-pvdz") and std.mem.eql(u8, format, "json")) {return @embedFile("basis/cc-pvdz.json");}
            else if (std.mem.eql(u8, lower, "cc-pvdz") and std.mem.eql(u8, format, "g94" )) {return @embedFile("basis/cc-pvdz.g94" );}
            else if (std.mem.eql(u8, lower, "cc-pvtz") and std.mem.eql(u8, format, "json")) {return @embedFile("basis/cc-pvtz.json");}
            else if (std.mem.eql(u8, lower, "cc-pvtz") and std.mem.eql(u8, format, "g94" )) {return @embedFile("basis/cc-pvtz.g94" );}

            else return error.BasisNameNotFound;
        }

        /// Get the basis set for the given system and name.
        pub fn get(system: System(T), name: []const u8, allocator: std.mem.Allocator) !std.ArrayList(ContractedGaussian(T)) {
            var basis = std.ArrayList(ContractedGaussian(T)).init(allocator); const parsed = try std.json.parseFromSlice(std.json.Value, allocator, try embedded(name, "json", allocator), .{});

            for (0..system.atoms.rows) |i| {

                var an: [3]u8 = undefined; const ans = std.fmt.formatIntBuf(&an, @as(u8, @intFromFloat(system.atoms.at(i))), 10, std.fmt.Case.lower, .{});

                if (parsed.value.object.get("elements").?.object.get(an[0..ans]) == null) return error.AtomNotFoundInBasis;

                const shells = parsed.value.object.get("elements").?.object.get(an[0..ans]).?.object.get("electron_shells").?.array.items;

                for (shells) |shell| {

                    const am        = shell.object.get("angular_momentum").?.array.items[0].integer;
                    const exponents = shell.object.get("exponents"       ).?.array.items           ;
                    const coefs     = shell.object.get("coefficients"    ).?.array.items           ;

                    var a     = try allocator.alloc(T, 3            ); defer     allocator.free(a);
                    var c     = try allocator.alloc(T, exponents.len); defer     allocator.free(c);
                    var alpha = try allocator.alloc(T, exponents.len); defer allocator.free(alpha);

                    for (coefs) |coef| {

                        for (0..alpha.len) |k| c[k] = try std.fmt.parseFloat(T, coef.array.items[k].string);

                        for (0..@as(usize, @intCast(am + 1))) |lx| {
                            for (0..@as(usize, @intCast(am + 1)) - lx) |ly| {

                                a[0] = @as(T, @floatFromInt(lx)); a[1] = @as(T, @floatFromInt(ly)); a[2] = @as(T, @floatFromInt(am)) - @as(T, @floatFromInt(lx + ly));

                                for (0..alpha.len) |k| alpha[k] = try std.fmt.parseFloat(T, exponents[k].string);

                                try basis.append(try ContractedGaussian(T).init(system.getCoords(i), .{a[0], a[1], a[2]}, c, alpha, allocator));
                            }
                        }
                    }
                }
            }

            parsed.deinit(); return basis;
        }
    };
}
