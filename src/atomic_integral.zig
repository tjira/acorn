//! Module for calculating integrals over atomic basis functions without doing any other calculation.

const std = @import("std"); const cwp = @import("cwrapper.zig");

const hlp = @import("helper.zig");
const inp = @import("input.zig" );
const sys = @import("system.zig");

const Basis  = @import("basis.zig" ).Basis ;
const Vector = @import("vector.zig").Vector;
const Matrix = @import("matrix.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

/// The main function for generating the integrals.
pub fn run(comptime T: type, opt: inp.AtomicIntegralOptions(T), print: bool, allocator: std.mem.Allocator) !void {
    var system = try sys.load(T, opt.system, opt.system_file, allocator); defer system.deinit();

    var basis = try Basis(T).array(system, opt.basis, allocator); defer basis.deinit();

    var nbf: usize = 0; var npgs: usize = 0;

    {
        var i: usize = 0; while (i < basis.rows) : (i += 2 * @as(usize, @intFromFloat(basis.at(i))) + 5) {
            const cgs: usize = @as(usize, @intFromFloat((basis.at(i + 1) + 1) * (basis.at(i + 1) + 2))) / 2; nbf += cgs; npgs += @as(usize, @intFromFloat(basis.at(i))) * cgs;
        }
    }

    if (print) try hlp.print(std.fs.File.stdout(), "\n# OF CONTRACTED GAUSSIAN SHELLS: {d}\n", .{nbf });
    if (print) try hlp.print(std.fs.File.stdout(),   "# OF PRIMITIVE  GAUSSIAN SHELLS: {d}\n", .{npgs});

    var timer = try std.time.Timer.start(); if (print) try hlp.print(std.fs.File.stdout(), "\n", .{});

    {

        const mem = try hlp.memFormat(T, nbf * nbf, allocator); defer allocator.free(mem);

        if (print) try hlp.print(std.fs.File.stdout(), "# OVERLAP INTEGRALS ({s}): ", .{mem});

        var S_AO = try Matrix(T).init(nbf, nbf, allocator); defer S_AO.deinit(); cwp.Libint(T).overlap(&S_AO, system, basis);

        if (print) try hlp.print(std.fs.File.stdout(), "{D}\n", .{timer.read()});

        if (print and opt.print.overlap) try S_AO.print(std.fs.File.stdout());

        if (opt.write.overlap != null) try S_AO.write(opt.write.overlap.?);
    }

    timer = try std.time.Timer.start();

    {

        const mem = try hlp.memFormat(T, nbf * nbf, allocator); defer allocator.free(mem);

        if (print) try hlp.print(std.fs.File.stdout(), "# KINETIC INTEGRALS ({s}): ", .{mem});

        var T_AO = try Matrix(T).init(nbf, nbf, allocator); defer T_AO.deinit(); cwp.Libint(T).kinetic(&T_AO, system, basis);

        if (print) try hlp.print(std.fs.File.stdout(), "{D}\n", .{timer.read()});

        if (print and opt.print.kinetic) try T_AO.print(std.fs.File.stdout());

        if (opt.write.kinetic != null) try T_AO.write(opt.write.kinetic.?);
    }

    timer = try std.time.Timer.start();

    {

        const mem = try hlp.memFormat(T, nbf * nbf, allocator); defer allocator.free(mem);

        if (print) try hlp.print(std.fs.File.stdout(), "# NUCLEAR INTEGRALS ({s}): ", .{mem});

        var V_AO = try Matrix(T).init(nbf, nbf, allocator); defer V_AO.deinit(); cwp.Libint(T).nuclear(&V_AO, system, basis);

        if (print) try hlp.print(std.fs.File.stdout(), "{D}\n", .{timer.read()});

        if (print and opt.print.nuclear) try V_AO.print(std.fs.File.stdout());

        if (opt.write.nuclear != null) try V_AO.write(opt.write.nuclear.?);
    }

    timer = try std.time.Timer.start();

    {

        const mem = try hlp.memFormat(T, nbf * nbf * nbf * nbf, allocator); defer allocator.free(mem);

        if (print) try hlp.print(std.fs.File.stdout(), "# COULOMB INTEGRALS ({s}): ", .{mem});

        var J_AO = try Tensor(T).init(&[_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_AO.deinit(); cwp.Libint(T).coulomb(&J_AO, system, basis);

        if (print) try hlp.print(std.fs.File.stdout(), "{D}\n", .{timer.read()});

        if (print and opt.print.coulomb) try J_AO.print(std.fs.File.stdout());

        if (opt.write.coulomb != null) try J_AO.write(opt.write.coulomb.?);
    }
}
