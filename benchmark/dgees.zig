const std = @import("std"); const acorn = @import("acorn");

pub const Args = struct {
    A: acorn.Matrix(f64), AT: acorn.Matrix(f64), AQ: acorn.Matrix(f64),

    pub fn init(allocator: std.mem.Allocator, D: usize) !Args {
        return .{
            .A   = try acorn.Matrix(f64).init(D, D, allocator),
            .AT  = try acorn.Matrix(f64).init(D, D, allocator),
            .AQ  = try acorn.Matrix(f64).init(D, D, allocator),
        };
    }

    pub fn deinit(self: *Args) void {
        self.A.deinit(); self.AT.deinit(); self.AQ.deinit();
    }

    pub fn function(self: *Args) !void {
        try acorn.cwrapper.Lapack(f64).dgees(&self.AT, &self.AQ, self.A);
    }

    pub fn print(self: Args) !void {
        try std.io.getStdOut().writer().print("MATRIX A:\n",        .{}); try self.A.print( std.io.getStdOut().writer());
        try std.io.getStdOut().writer().print("SCHUR FORM:\n",      .{}); try self.AT.print(std.io.getStdOut().writer());
        try std.io.getStdOut().writer().print("ORTHOGONAL FORM:\n", .{}); try self.AQ.print(std.io.getStdOut().writer());
    }

    pub fn randomize(self: *Args, seed: usize) !void {
        self.A.randn(0, 1, seed);
    }
};
