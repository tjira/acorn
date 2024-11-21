const std = @import("std"); const Complex = std.math.Complex; const gsl = @cImport(@cInclude("gsl/gsl_eigen.h"));

const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

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
        pub fn complex(self: Matrix(T)) !Matrix(Complex(T)) {
            var other = try Matrix(Complex(T)).init(self.rows, self.cols, self.allocator); for (0..self.data.len) |i| other.data[i] = Complex(T).init(self.data[i], 0); return other;
        }

        pub fn eq(self: Matrix(T), other: Matrix(T), tol: T) bool {
            for (0..self.data.len) |i| {if (@abs(other.data[i] - self.data[i]) > tol) return false;} return true;
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
        pub fn rowptr(self: Matrix(T), i: usize) Matrix(T) {
            return Matrix(T){.data = self.data[i * self.cols..(i + 1) * self.cols], .rows = 1, .cols = self.cols, .allocator = self.allocator};
        }
        pub fn vectorptr(self: Matrix(T)) Vector(T) {
            return Vector(T){.data = self.data[0..], .rows = self.rows * self.cols, .allocator = self.allocator};
        }

        pub fn fill(self: Matrix(T), value: T) void {
            for (0..self.data.len) |i| self.data[i] = value;
        }
        pub fn identity(self: Matrix(T)) void {
            self.fill(0); for (0..self.rows) |i| self.ptr(i, i).* = 1;
        }
        pub fn linspace(self: Matrix(T), start: T, end: T) void {
            for (0..self.data.len) |i| self.data[i] = start + asfloat(T, i) * (end - start) / asfloat(T, self.rows * self.cols - 1);
        }
        pub fn setcol(self: Matrix(T), j: usize, other: Matrix(T)) void {
            for (0..self.rows) |i| self.ptr(i, j).* = other.at(i, 0);
        }
        pub fn transpose(self: Matrix(T)) void {
            for (0..self.rows) |i| for (i + 1..self.cols) |j| std.mem.swap(T, self.ptr(i, j), self.ptr(j, i));
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

pub fn eig(comptime T: type, J: *Matrix(T), C: *Matrix(T), A: Matrix(T), GSLW: *gsl.gsl_eigen_nonsymmv_workspace) void {
    var GSLA = gsl.gsl_matrix_view_array(A.data.ptr, A.rows, A.cols);
    var GSLC = gsl.gsl_matrix_view_array(C.data.ptr, A.rows, A.cols);
    var GSLJ = gsl.gsl_vector_view_array(J.data.ptr, A.rows        );

    J.fill(0); _ = gsl.gsl_eigen_nonsymmv(&GSLA.matrix, &GSLJ.vector, &GSLC.matrix, GSLW);

    _ = gsl.gsl_eigen_nonsymmv_sort(&GSLJ.vector, &GSLC.matrix, gsl.GSL_EIGEN_SORT_VAL_ASC);

    for (1..A.rows) |i| {J.ptr(i, i).* = J.at(0, i); J.ptr(0, i).* = 0;}
}

pub fn eigh(comptime T: type, J: *Matrix(T), C: *Matrix(T), A: Matrix(T), GSLW: *gsl.gsl_eigen_symmv_workspace) void {
    var GSLA = gsl.gsl_matrix_view_array(A.data.ptr, A.rows, A.cols);
    var GSLC = gsl.gsl_matrix_view_array(C.data.ptr, A.rows, A.cols);
    var GSLJ = gsl.gsl_vector_view_array(J.data.ptr, A.rows        );

    J.fill(0); _ = gsl.gsl_eigen_symmv(&GSLA.matrix, &GSLJ.vector, &GSLC.matrix, GSLW);

    _ = gsl.gsl_eigen_symmv_sort(&GSLJ.vector, &GSLC.matrix, gsl.GSL_EIGEN_SORT_VAL_ASC);

    for (1..A.rows) |i| {J.ptr(i, i).* = J.at(0, i); J.ptr(0, i).* = 0;}
}

pub fn hjoin(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    for (0..A.rows) |i| {for (0..A.cols) |j| C.ptr(i, j).* = A.at(i, j); for (0..B.cols) |j| C.ptr(i, A.cols + j).* = B.at(i, j);}
}

pub fn transpose(comptime T: type, B: *Matrix(T), A: Matrix(T)) void {
    for (0..A.rows) |i| for (0..A.cols) |j| {B.ptr(j, i).* = A.at(i, j);};
}
