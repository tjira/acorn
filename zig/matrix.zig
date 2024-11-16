const std = @import("std");

const swap = @import("helper.zig").swap;

pub fn Matrix(comptime T: type) type {
    return struct {
        data: []T, rows: usize, cols: usize, allocator: std.mem.Allocator,

        pub fn init(rows: usize, cols: usize, allocator: std.mem.Allocator) !Matrix(T) {
            return Matrix(T){.data = try allocator.alloc(T, rows * cols), .rows = rows, .cols = cols, .allocator = allocator};
        }
        pub fn clone(self: Matrix(T)) !Matrix(T) {
            const other = try Matrix(T).init(self.rows, self.cols, self.allocator); other.set(self.data); return other;
        }
        pub fn deinit(self: Matrix(T)) void {
            self.allocator.free(self.data);
        }

        pub fn eq(self: Matrix(T), other: Matrix(T), tol: T) bool {
            for (0..self.data.len) |i| {if (@abs(other.data[i] - self.data[i]) > tol) return false;} return true;
        }
        pub fn isdiag(self: Matrix(T), tol: T) bool {
            for (0..self.rows) |i| {for (i + 1..self.cols) |j| if (@abs(self.at(i, j)) > tol or @abs(self.at(j, i)) > tol) return false;} return true;
        }

        pub fn at(self: Matrix(T), i: usize, j: usize) T {
            return self.data[i * self.cols + j];
        }
        pub fn ptr(self: Matrix(T), i: usize, j: usize) *T {
            return &self.data[i * self.cols + j];
        }
        pub fn set(self: Matrix(T), data: []const T) void {
            @memcpy(self.data[0..data.len], data);
        }

        pub fn fill(self: Matrix(T), value: T) void {
            for (0..self.data.len) |i| self.data[i] = value;
        }
        pub fn identity(self: Matrix(T)) void {
            self.fill(0); for (0..self.rows) |i| self.ptr(i, i).* = 1;
        }
        pub fn linspace(self: Matrix(T), start: T, end: T) void {
            for (0..self.data.len) |i| self.data[i] = start + @as(T, @floatFromInt(i)) * (end - start) / @as(T, @floatFromInt(self.rows * self.cols - 1));
        }

        pub fn hjoin(self: Matrix(T), output: *Matrix(T), other: Matrix(T)) void {
            for (0..self.rows) |i| {for (0..self.cols) |j| output.ptr(i, j).* = self.at(i, j); for (0..other.cols) |j| output.ptr(i, self.cols + j).* = other.at(i, j);}
        }

        pub fn print(self: Matrix(T), device: anytype) !void {
            try device.print("{d} {d}\n", .{self.rows, self.cols}); for (self.data, 1..) |e, i| try device.print("{d:20.14}{s}", .{e, if(i % self.cols == 0) "\n" else " "});
        }
        pub fn write(self: Matrix(T), path: []const u8) !void {
            const file = try std.fs.cwd().createFile(path, .{}); defer file.close(); try self.print(file.writer());
        }
    };
}

pub fn mm(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    C.fill(0); for (0..A.rows) |i| for (0..B.cols) |j| for (0..A.cols) |k| {C.ptr(i, j).* += A.at(i, k) * B.at(k, j);};
}

pub fn eigh(comptime T: type, J: *Matrix(T), C: *Matrix(T), A: Matrix(T), tol: T, T1: *Matrix(T), T2: *Matrix(T)) void {
    @memcpy(J.data, A.data); var maxi: usize = undefined; var maxj: usize = undefined; var maxv: T = undefined; var phi: T = undefined; 

    while (!J.isdiag(tol)) : ({maxi = 0; maxj = 1; maxv = J.at(maxi, maxj);}) {

        for (0..A.rows) |i| for (i + 1..A.cols) |j| if (J.at(i, j) > maxv) {maxi = i; maxj = j; maxv = J.at(i, j);};

        phi = 0.5 * std.math.atan(2 * maxv / (J.at(maxi, maxi) - J.at(maxj, maxj))); C.identity();

        C.ptr(maxi, maxi).* = std.math.cos(phi); C.ptr(maxj, maxj).* =  C.at(maxi, maxi);
        C.ptr(maxj, maxi).* = std.math.sin(phi); C.ptr(maxi, maxj).* = -C.at(maxj, maxi);

        mm(T, T1, J.*, C.*); transpose(T, T2, C.*); mm(T, J, T2.*, T1.*);
    }

    for (0..J.rows) |i| for (i..J.cols) |j| if (J.at(i, i) > J.at(j, j)) {
        swap(J.ptr(j, j), J.ptr(i, i)); for (0..C.rows) |k| swap(C.ptr(k, i), C.ptr(k, j));
    };
}

pub fn exph(comptime T: type, E: *Matrix(T), J: Matrix(T), C: Matrix(T), T1: *Matrix(T), T2: *Matrix(T)) void {
    @memcpy(T1.data, J.data); for (0..J.rows) |i| T1.ptr(i, i).* = std.math.exp(J.at(i, i));

    mm(T, T2, C, T1.*); transpose(T, T1, C); mm(T, E, T2.*, T1.*);
}

pub fn logh(comptime T: type, L: *Matrix(T), J: Matrix(T), C: Matrix(T), T1: *Matrix(T), T2: *Matrix(T)) void {
    @memcpy(T1.data, J.data); for (0..J.rows) |i| T1.ptr(i, i).* = std.math.log10(J.at(i, i)) / std.math.log10e;

    mm(T, T2, C, T1.*); transpose(T, T1, C); mm(T, L, T2.*, T1.*);
}

pub fn transpose(comptime T: type, B: *Matrix(T), A: Matrix(T)) void {
    for (0..A.rows) |i| for (0..A.cols) |j| {B.ptr(i, j).* = A.at(j, i);};
}
