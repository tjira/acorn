const std = @import("std");

pub fn Matrix(comptime T: type) type {
    return struct {
        data: []T, rows: usize, cols: usize, allocator: std.mem.Allocator,

        pub fn init(rows: usize, cols: usize, allocator: std.mem.Allocator) !Matrix(T) {
            return Matrix(T){.data = try allocator.alloc(T, rows * cols), .rows = rows, .cols = cols, .allocator = allocator};
        }
        pub fn clone(self: Matrix(T)) !Matrix(T) {
            const other = try Matrix(T).init(self.rows, self.cols, self.allocator); @memcpy(other.data, self.data); return other;
        }
        pub fn deinit(self: Matrix(T)) void {
            self.allocator.free(self.data);
        }

        pub fn at(self: Matrix(T), i: usize, j: usize) T {
            return self.data[i * self.cols + j];
        }
        pub fn ptr(self: Matrix(T), i: usize, j: usize) *T {
            return &self.data[i * self.cols + j];
        }

        pub fn fill(self: Matrix(T), value: T) void {
            for (0..self.data.len) |i| self.data[i] = value;
        }
    };
}

pub fn write(comptime T: type, path: []const u8, matrix: Matrix(T)) !void {
    const file = try std.fs.cwd().createFile(path, .{}); defer file.close(); const writer = file.writer();
    try writer.print("{d} {d}\n", .{matrix.rows, matrix.cols});
    for (0..matrix.rows) |i| {
        for (0..matrix.cols) |j| {
            try writer.print("{d:20.14}{s}", .{matrix.at(i, j), if(j == matrix.cols - 1) "\n" else " "});
        }
    }
}
