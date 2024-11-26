const std = @import("std"); const Complex = std.math.Complex; const gsl_eigen = @cImport(@cInclude("gsl/gsl_eigen.h"));

const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;

pub fn Matrix(comptime T: type) type {
    return struct {
        data: []T, rows: usize, cols: usize, allocator: std.mem.Allocator,

        pub fn init(rows: usize, cols: usize, allocator: std.mem.Allocator) !Matrix(T) {
            return Matrix(T){.data = try allocator.alloc(T, rows * cols), .rows = rows, .cols = cols, .allocator = allocator};
        }
        pub fn deinit(self: Matrix(T)) void {
            self.allocator.free(self.data);
        }

        pub fn clone(self: Matrix(T)) !Matrix(T) {
            const other = try Matrix(T).init(self.rows, self.cols, self.allocator); @memcpy(other.data, self.data); return other;
        }
        pub fn complex(self: Matrix(T)) !Matrix(Complex(T)) {
            var other = try Matrix(Complex(T)).init(self.rows, self.cols, self.allocator); for (0..self.data.len) |i| other.data[i] = Complex(T).init(self.data[i], 0); return other;
        }

        pub fn at(self: Matrix(T), i: usize, j: usize) T {
            return self.data[i * self.cols + j];
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

        pub fn print(self: Matrix(T), device: anytype) !void {
            try device.print("{d} {d}\n", .{self.rows, self.cols}); for (self.data, 1..) |e, i| try device.print("{d:20.14}{s}", .{e, if(i % self.cols == 0) "\n" else " "});
        }
        pub fn write(self: Matrix(T), path: []const u8) !void {
            const file = try std.fs.cwd().createFile(path, .{}); defer file.close(); try self.print(file.writer());
        }
    };
}

pub fn mm(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if ( @hasField(T, "re")) {C.fill(T.init(0, 0)); for (0..A.rows) |i| for (0..B.cols) |j| for (0..A.cols) |k| {C.ptr(i, j).* = C.at(i, j).add(A.at(i, k).mul(B.at(k, j)));};}
    if (!@hasField(T, "re")) {C.fill(0)           ; for (0..A.rows) |i| for (0..B.cols) |j| for (0..A.cols) |k| {C.ptr(i, j).* = C.at(i, j) +   A.at(i, k) *   B.at(k, j)  ;};}
}

pub fn mam(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if ( @hasField(T, "re")) {C.fill(T.init(0, 0)); for (0..A.cols) |i| for (0..B.cols) |j| for (0..A.cols) |k| {C.ptr(i, j).* = C.at(i, j).add(A.at(k, i).conjugate().mul(B.at(k, j)));};}
    if (!@hasField(T, "re")) {C.fill(0)           ; for (0..A.cols) |i| for (0..B.cols) |j| for (0..A.cols) |k| {C.ptr(i, j).* = C.at(i, j) +   A.at(k, i) *               B.at(k, j)  ;};}
}

pub fn mma(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if ( @hasField(T, "re")) {C.fill(T.init(0, 0)); for (0..A.rows) |i| for (0..B.rows) |j| for (0..A.cols) |k| {C.ptr(i, j).* = C.at(i, j).add(A.at(i, k).mul(B.at(j, k).conjugate()));};}
    if (!@hasField(T, "re")) {C.fill(0)           ; for (0..A.rows) |i| for (0..B.rows) |j| for (0..A.cols) |k| {C.ptr(i, j).* = C.at(i, j) +   A.at(i, k) *   B.at(j, k)              ;};}
}

pub fn mtma(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if ( @hasField(T, "re")) {C.fill(T.init(0, 0)); for (0..A.cols) |i| for (0..B.rows) |j| for (0..A.cols) |k| {C.ptr(j, i).* = C.at(j, i).add(A.at(i, k).mul(B.at(j, k).conjugate()));};}
    if (!@hasField(T, "re")) {C.fill(0)           ; for (0..A.cols) |i| for (0..B.rows) |j| for (0..A.cols) |k| {C.ptr(j, i).* = C.at(j, i) +   A.at(i, k) *   B.at(j, k)              ;};}
}

pub fn mamt(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if ( @hasField(T, "re")) {C.fill(T.init(0, 0)); for (0..A.cols) |i| for (0..B.rows) |j| for (0..A.cols) |k| {C.ptr(j, i).* = C.at(j, i).add(A.at(k, i).conjugate().mul(B.at(j, k)));};}
    if (!@hasField(T, "re")) {C.fill(0)           ; for (0..A.cols) |i| for (0..B.rows) |j| for (0..A.cols) |k| {C.ptr(j, i).* = C.at(j, i) +   A.at(k, i) *               B.at(j, k)  ;};}
}

pub fn eigh(comptime T: type, J: *Matrix(T), C: *Matrix(T), A: Matrix(T), AT: *Matrix(T), GSLW: *gsl_eigen.gsl_eigen_symmv_workspace) void {
    var sumsq: T = 0; for (0..A.rows) |i| for (i + 1..A.cols) |j| {sumsq += 2 * A.at(i, j) * A.at(i, j); sumsq = std.math.sqrt(sumsq);};

    J.fill(0); C.fill(0); @memcpy(AT.data, A.data);

    var GSLA = gsl_eigen.gsl_matrix_view_array(AT.data.ptr, A.rows, A.cols);
    var GSLC = gsl_eigen.gsl_matrix_view_array( C.data.ptr, A.rows, A.cols);
    var GSLJ = gsl_eigen.gsl_vector_view_array( J.data.ptr, A.rows        );

    if (sumsq >  1e-14) _ = gsl_eigen.gsl_eigen_symmv(&GSLA.matrix, &GSLJ.vector, &GSLC.matrix, GSLW);
    if (sumsq <= 1e-14) {for (0..A.rows) |i| J.ptr(0, i).* = A.at(i, i); C.identity();}

    _ = gsl_eigen.gsl_eigen_symmv_sort(&GSLJ.vector, &GSLC.matrix, gsl_eigen.GSL_EIGEN_SORT_VAL_ASC);

    for (1..A.rows) |i| {J.ptr(i, i).* = J.at(0, i); J.ptr(0, i).* = 0;}
}

pub fn hjoin(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    for (0..A.rows) |i| {for (0..A.cols) |j| C.ptr(i, j).* = A.at(i, j); for (0..B.cols) |j| C.ptr(i, A.cols + j).* = B.at(i, j);}
}
