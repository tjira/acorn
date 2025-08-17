const std = @import("std"); const acorn = @import("acorn");

pub const Args = struct {
    A: acorn.Tensor(f64), B: acorn.Tensor(f64), C: acorn.Tensor(f64), allocator: std.mem.Allocator,

    pub fn init(allocator: std.mem.Allocator, D: usize) !Args {
        return .{
            .A = try acorn.Tensor(f64).init(&[_]usize{D, D, D, D}, allocator),
            .B = try acorn.Tensor(f64).init(&[_]usize{D, D      }, allocator),
            .C = try acorn.Tensor(f64).init(&[_]usize{D, D      }, allocator),
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Args) void {
        self.A.deinit(); self.B.deinit(); self.C.deinit();
    }

    pub fn function(self: *Args) !void {
        try acorn.tensor.contract(f64, &self.C, self.A, &[_]usize{1, 2}, self.B, &[_]usize{0, 1}, self.allocator);
    }

    pub fn print(self: Args) !void {
        try std.io.getStdOut().writer().print("MATRIX C:\n", .{}); try self.C.matrix(1).print(std.io.getStdOut().writer());
    }

    pub fn randomize(self: *Args, seed: usize) !void {
        self.A.randn(0, 1, 2 * seed + 0);
        self.B.randn(0, 1, 2 * seed + 1);
    }
};
