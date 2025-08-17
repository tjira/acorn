const std = @import("std"); const acorn = @import("acorn");

pub const Args = struct {
    A: acorn.Matrix(f64), AJ: acorn.Matrix(f64), AC: acorn.Matrix(f64),

    pub fn init(allocator: std.mem.Allocator, D: usize) !Args {
        return .{
            .A  = try acorn.Matrix(f64).init(D, D, allocator),
            .AJ = try acorn.Matrix(f64).init(D, D, allocator),
            .AC = try acorn.Matrix(f64).init(D, D, allocator),
        };
    }

    pub fn deinit(self: *Args) void {
        self.A.deinit(); self.AJ.deinit(); self.AC.deinit();
    }

    pub fn function(self: *Args) !void {
        try acorn.cwrapper.Lapack(f64).dsyevd(&self.AJ, &self.AC, self.A);
    }

    pub fn print(self: Args) !void {
        try std.io.getStdOut().writer().print("MATRIX A:\n",     .{}); try self.A.print( std.io.getStdOut().writer());
        try std.io.getStdOut().writer().print("EIGENVALUES:\n",  .{}); try self.AJ.print(std.io.getStdOut().writer());
        try std.io.getStdOut().writer().print("EIGENVECTORS:\n", .{}); try self.AC.print(std.io.getStdOut().writer());
    }

    pub fn randomize(self: *Args, seed: usize) !void {
        self.A.randn(0, 1, seed); try self.A.symmetrize();
    }
};
