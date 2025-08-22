//! Basis library.

const std = @import("std");

const System = @import("system.zig").System;
const Vector = @import("vector.zig").Vector;

/// Basis struct.
pub fn Basis(comptime T: type) type {
    return struct {

        /// Return the basis set as an array of consecutive numbers.
        pub fn array(system: System(T), name: []const u8, allocator: std.mem.Allocator) !Vector(T){
            var basis = try Vector(T).init(0, allocator); const parsed = try std.json.parseFromSlice(std.json.Value, allocator, try embedded(name, "json", allocator), .{});

            for (0..system.atoms.rows) |i| {

                var an: [3]u8 = undefined; const ans = (try std.fmt.bufPrint(&an, "{d}", .{@as(u8, @intFromFloat(system.atoms.at(i)))})).len;

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

            parsed.deinit(); return basis;
        }

        /// Returns the basis set for the given name and format.
        pub fn embedded(basis: []const u8, format: []const u8, allocator: std.mem.Allocator) ![]const u8 {
            const lower = try allocator.alloc(u8, basis.len); defer allocator.free(lower); _ = std.ascii.lowerString(lower, basis);

            if (std.mem.eql(u8, lower, "sto-2g"     ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/sto-2g.json"     );}
            if (std.mem.eql(u8, lower, "sto-3g"     ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/sto-3g.json"     );}
            if (std.mem.eql(u8, lower, "sto-4g"     ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/sto-4g.json"     );}
            if (std.mem.eql(u8, lower, "sto-5g"     ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/sto-5g.json"     );}
            if (std.mem.eql(u8, lower, "sto-6g"     ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/sto-6g.json"     );}
            if (std.mem.eql(u8, lower, "3-21g"      ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/3-21g.json"      );}
            if (std.mem.eql(u8, lower, "6-31g"      ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-31g.json"      );}
            if (std.mem.eql(u8, lower, "6-31g*"     ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-31gs.json"     );}
            if (std.mem.eql(u8, lower, "6-31g**"    ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-31gss.json"    );}
            if (std.mem.eql(u8, lower, "6-31+g"     ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-31pg.json"     );}
            if (std.mem.eql(u8, lower, "6-31+g*"    ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-31pgs.json"    );}
            if (std.mem.eql(u8, lower, "6-31+g**"   ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-31pgss.json"   );}
            if (std.mem.eql(u8, lower, "6-31++g"    ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-31ppg.json"    );}
            if (std.mem.eql(u8, lower, "6-31++g*"   ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-31ppgs.json"   );}
            if (std.mem.eql(u8, lower, "6-31++g**"  ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-31ppgss.json"  );}
            if (std.mem.eql(u8, lower, "6-311g"     ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-311g.json"     );}
            if (std.mem.eql(u8, lower, "6-311g*"    ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-311gs.json"    );}
            if (std.mem.eql(u8, lower, "6-311g**"   ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-311gss.json"   );}
            if (std.mem.eql(u8, lower, "6-311+g"    ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-311pg.json"    );}
            if (std.mem.eql(u8, lower, "6-311+g*"   ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-311pgs.json"   );}
            if (std.mem.eql(u8, lower, "6-311+g**"  ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-311pgss.json"  );}
            if (std.mem.eql(u8, lower, "6-311++g"   ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-311ppg.json"   );}
            if (std.mem.eql(u8, lower, "6-311++g*"  ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-311ppgs.json"  );}
            if (std.mem.eql(u8, lower, "6-311++g**" ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/6-311ppgss.json" );}
            if (std.mem.eql(u8, lower, "def2-svp"   ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/def2-svp.json"   );}
            if (std.mem.eql(u8, lower, "def2-svpd"  ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/def2-svpd.json"  );}
            if (std.mem.eql(u8, lower, "def2-tzvp"  ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/def2-tzvp.json"  );}
            if (std.mem.eql(u8, lower, "def2-tzvpd" ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/def2-tzvpd.json" );}
            if (std.mem.eql(u8, lower, "def2-tzvpp" ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/def2-tzvpp.json" );}
            if (std.mem.eql(u8, lower, "def2-tzvppd") and std.mem.eql(u8, format, "json")) {return @embedFile("basis/def2-tzvppd.json");}
            if (std.mem.eql(u8, lower, "def2-qzvp"  ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/def2-qzvp.json"  );}
            if (std.mem.eql(u8, lower, "def2-qzvpd" ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/def2-qzvpd.json" );}
            if (std.mem.eql(u8, lower, "def2-qzvpp" ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/def2-qzvpp.json" );}
            if (std.mem.eql(u8, lower, "def2-qzvppd") and std.mem.eql(u8, format, "json")) {return @embedFile("basis/def2-qzvppd.json");}
            if (std.mem.eql(u8, lower, "cc-pvdz"    ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/cc-pvdz.json"    );}
            if (std.mem.eql(u8, lower, "cc-pvtz"    ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/cc-pvtz.json"    );}
            if (std.mem.eql(u8, lower, "cc-pvqz"    ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/cc-pvqz.json"    );}
            if (std.mem.eql(u8, lower, "cc-pv5z"    ) and std.mem.eql(u8, format, "json")) {return @embedFile("basis/cc-pv5z.json"    );}
            if (std.mem.eql(u8, lower, "aug-cc-pvdz") and std.mem.eql(u8, format, "json")) {return @embedFile("basis/aug-cc-pvdz.json");}
            if (std.mem.eql(u8, lower, "aug-cc-pvtz") and std.mem.eql(u8, format, "json")) {return @embedFile("basis/aug-cc-pvtz.json");}
            if (std.mem.eql(u8, lower, "aug-cc-pvqz") and std.mem.eql(u8, format, "json")) {return @embedFile("basis/aug-cc-pvqz.json");}
            if (std.mem.eql(u8, lower, "aug-cc-pv5z") and std.mem.eql(u8, format, "json")) {return @embedFile("basis/aug-cc-pv5z.json");}

            return error.BasisNameNotFound;
        }
    };
}
