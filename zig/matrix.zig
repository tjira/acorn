const std = @import("std");

const gsl = @cImport(@cInclude("gsl/gsl_eigen.h"));

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
            var sumsq: T = 0; for (0..self.rows) |i| {for (i + 1..self.cols) |j| sumsq += self.at(i, j) * self.at(i, j) + self.at(j, i) * self.at(j, i);} return sumsq < tol;
        }

        pub fn at(self: Matrix(T), i: usize, j: usize) T {
            return self.data[i * self.cols + j];
        }
        pub fn col(self: Matrix(T), output: *Matrix(T), j: usize) void {
            for (0..self.rows) |i| output.data[i] = self.at(i, j);
        }
        pub fn ptr(self: Matrix(T), i: usize, j: usize) *T {
            return &self.data[i * self.cols + j];
        }
        pub fn row(self: Matrix(T), output: *Matrix(T), i: usize) void {
            @memcpy(output.data, self.data[i * self.cols..(i + 1) * self.cols]);
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
        pub fn setcol(self: Matrix(T), j: usize, other: Matrix(T)) void {
            for (0..self.rows) |i| self.ptr(i, j).* = other.at(i, 0);
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

pub fn eigh(comptime T: type, J: *Matrix(T), C: *Matrix(T), A: Matrix(T)) void {
    const m    = gsl.gsl_matrix_alloc(A.rows, A.cols);
    const eval = gsl.gsl_vector_alloc(A.rows        );
    const evec = gsl.gsl_matrix_alloc(A.rows, A.cols);

    for (0..A.rows) |i| for (0..A.cols) |j| gsl.gsl_matrix_set(m, i, j, A.at(i, j));

    const w = gsl.gsl_eigen_symmv_alloc(A.rows);

    _ = gsl.gsl_eigen_symmv(m, eval, evec, w);

    _ = gsl.gsl_eigen_symmv_sort(eval, evec, gsl.GSL_EIGEN_SORT_VAL_DESC);

    for (0..A.rows) |i| for (0..A.cols) |j| {C.ptr(i, j).* = gsl.gsl_matrix_get(evec, i, j);};
    for (0..A.rows) |i| J.ptr(i, i).* = gsl.gsl_vector_get(eval, i);

    gsl.gsl_eigen_symmv_free(w);
    gsl.gsl_matrix_free(m);
}

pub fn transpose(comptime T: type, B: *Matrix(T), A: Matrix(T)) void {
    for (0..A.rows) |i| for (0..A.cols) |j| {B.ptr(i, j).* = A.at(j, i);};
}
