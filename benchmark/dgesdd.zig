const std = @import("std"); const acorn = @import("acorn");

pub const Args = struct {
    A: acorn.Matrix(f64), AU: acorn.Matrix(f64), AS: acorn.Matrix(f64), AVT: acorn.Matrix(f64),

    pub fn init(allocator: std.mem.Allocator, D: usize) !Args {
        return .{
            .A   = try acorn.Matrix(f64).init(D, D, allocator),
            .AU  = try acorn.Matrix(f64).init(D, D, allocator),
            .AS  = try acorn.Matrix(f64).init(D, D, allocator),
            .AVT = try acorn.Matrix(f64).init(D, D, allocator),
        };
    }

    pub fn deinit(self: *Args) void {
        self.A.deinit(); self.AU.deinit(); self.AS.deinit(); self.AVT.deinit();
    }

    pub fn function(self: *Args) !void {
        try acorn.cwrapper.Lapack(f64).dgesdd(&self.AU, &self.AS, &self.AVT, self.A);
    }

    pub fn print(self: Args) !void {
        try std.io.getStdOut().writer().print("MATRIX A:\n",     .{}); try   self.A.print(std.io.getStdOut().writer());
        try std.io.getStdOut().writer().print("MATRIX U:\n",     .{}); try  self.AU.print(std.io.getStdOut().writer());
        try std.io.getStdOut().writer().print("MATRIX SIGMA:\n", .{}); try  self.AS.print(std.io.getStdOut().writer());
        try std.io.getStdOut().writer().print("MATRIX V^T:\n",   .{}); try self.AVT.print(std.io.getStdOut().writer());
    }

    pub fn randomize(self: *Args, seed: usize) !void {
        self.A.randn(0, 1, seed);
    }
};
