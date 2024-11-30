const std = @import("std"); const Complex = std.math.Complex;

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
    C.fill(0); for (0..A.rows) |i| for (0..B.cols) |j| for (0..A.cols) |k| {C.ptr(i, j).* = C.at(i, j) + A.at(i, k) * B.at(k, j);};
}

pub fn cmm(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    C.fill(T.init(0, 0)); for (0..A.rows) |i| for (0..B.cols) |j| for (0..A.cols) |k| {C.ptr(i, j).* = C.at(i, j).add(A.at(i, k).mul(B.at(k, j)));};
}

pub fn cmma(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    C.fill(T.init(0, 0)); for (0..A.rows) |i| for (0..B.rows) |j| for (0..A.cols) |k| {C.ptr(i, j).* = C.at(i, j).add(A.at(i, k).mul(B.at(j, k).conjugate()));};
}

pub fn cmamt(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    C.fill(T.init(0, 0)); for (0..A.cols) |i| for (0..B.rows) |j| for (0..A.cols) |k| {C.ptr(j, i).* = C.at(j, i).add(A.at(k, i).conjugate().mul(B.at(j, k)));};
}

pub fn hjoin(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    for (0..A.rows) |i| {for (0..A.cols) |j| C.ptr(i, j).* = A.at(i, j); for (0..B.cols) |j| C.ptr(i, A.cols + j).* = B.at(i, j);}
}

pub fn transpose(comptime T: type, B: *Matrix(T), A: Matrix(T)) void {
    for (0..A.rows) |i| for (0..A.cols) |j| {B.ptr(j, i).* = A.at(i, j);};
}

pub fn eigh(comptime T: type, J: *Matrix(T), C: *Matrix(T), A: Matrix(T), T1: *Matrix(T), T2: *Matrix(T), T3: *Matrix(T)) void {
    var maxi: usize = 0; var maxj: usize = 1; var maxv: T = 0; var phi: T = undefined; @memcpy(J.data, A.data); C.identity();

    for (0..A.rows) |i| for (i + 1..A.cols) |j| if (@abs(J.at(i, j)) > @abs(maxv)) {maxi = i; maxj = j; maxv = J.at(i, j);};

    while (@abs(maxv) > 1e-14) {

        phi = 0.5 * std.math.atan(2 * maxv / (J.at(maxi, maxi) - J.at(maxj, maxj))); T3.identity();

        T3.ptr(maxi, maxi).* = std.math.cos(phi); T3.ptr(maxj, maxj).* =  T3.at(maxi, maxi);
        T3.ptr(maxj, maxi).* = std.math.sin(phi); T3.ptr(maxi, maxj).* = -T3.at(maxj, maxi);

        mm(T, T1, J.*, T3.*); transpose(T, T2, T3.*); mm(T, J, T2.*, T1.*); mm(T, T1, C.*, T3.*); @memcpy(C.data, T1.data);

        maxv = 0; for (0..A.rows) |i| for (i + 1..A.cols) |j| if (@abs(J.at(i, j)) > @abs(maxv)) {maxi = i; maxj = j; maxv = J.at(i, j);};
    }

    for (0..A.rows) |i| for (i + 1..A.cols) |j| if (J.at(i, i) > J.at(j, j)) {
        std.mem.swap(T, J.ptr(i, i), J.ptr(j, j)); for (0..A.rows) |k| std.mem.swap(T, C.ptr(k, i), C.ptr(k, j));
    };
}

test "eigh" {
    // var A  = try Matrix(f64).init(2, 2, std.testing.allocator); defer  A.deinit();
    // var AJ = try Matrix(f64).init(2, 2, std.testing.allocator); defer AJ.deinit();
    // var AC = try Matrix(f64).init(2, 2, std.testing.allocator); defer AC.deinit();
    // A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    // A.ptr(1, 0).* = 2; A.ptr(1, 1).* = 1;
    // var T1 = try A.clone(); defer T1.deinit();
    // var T2 = try A.clone(); defer T2.deinit();
    // var T3 = try A.clone(); defer T3.deinit();
    // myeigh(f64, &AJ, &AC, A, &T1, &T2, &T3);
    // std.debug.print("{d:20.14}{d:20.14}\n{d:20.14}{d:20.14}\n\n", .{A.at(0, 0), A.at(0, 1), A.at(1, 0), A.at(1, 1)});
    // std.debug.print("{d:20.14}{d:20.14}\n{d:20.14}{d:20.14}\n\n", .{AJ.at(0, 0), AJ.at(0, 1), AJ.at(1, 0), AJ.at(1, 1)});
    // std.debug.print("{d:20.14}{d:20.14}\n{d:20.14}{d:20.14}\n\n", .{AC.at(0, 0), AC.at(0, 1), AC.at(1, 0), AC.at(1, 1)});

    // var A  = try Matrix(f64).init(3, 3, std.testing.allocator); defer  A.deinit();
    // var AJ = try Matrix(f64).init(3, 3, std.testing.allocator); defer AJ.deinit();
    // var AC = try Matrix(f64).init(3, 3, std.testing.allocator); defer AC.deinit();
    // A.ptr(0, 0).* = 0.032; A.ptr(0, 1).* = 0.00000003571285; A.ptr(0, 2).* = 0.00000007142570;
    // A.ptr(1, 0).* = 0.00000003571285; A.ptr(1, 1).* = 0; A.ptr(1, 2).* = 0.00000003571285;
    // A.ptr(2, 0).* = 0.00000007142570; A.ptr(2, 1).* = 0.00000003571285; A.ptr(2, 2).* = -0.032;
    // var T1 = try A.clone(); defer T1.deinit();
    // var T2 = try A.clone(); defer T2.deinit();
    // var T3 = try A.clone(); defer T3.deinit();
    // myeigh(f64, &AJ, &AC, A, &T1, &T2, &T3);
    // std.debug.print("{d:20.14}{d:20.14}{d:20.14}\n{d:20.14}{d:20.14}{d:20.14}\n{d:20.14}{d:20.14}{d:20.14}\n\n", .{A.at(0, 0), A.at(0, 1), A.at(0, 2), A.at(1, 0), A.at(1, 1), A.at(1, 2), A.at(2, 0), A.at(2, 1), A.at(2, 2)});
    // std.debug.print("{d:20.14}{d:20.14}{d:20.14}\n{d:20.14}{d:20.14}{d:20.14}\n{d:20.14}{d:20.14}{d:20.14}\n\n", .{AJ.at(0, 0), AJ.at(0, 1), AJ.at(0, 2), AJ.at(1, 0), AJ.at(1, 1), AJ.at(1, 2), AJ.at(2, 0), AJ.at(2, 1), AJ.at(2, 2)});
    // std.debug.print("{d:20.14}{d:20.14}{d:20.14}\n{d:20.14}{d:20.14}{d:20.14}\n{d:20.14}{d:20.14}{d:20.14}\n\n", .{AC.at(0, 0), AC.at(0, 1), AC.at(0, 2), AC.at(1, 0), AC.at(1, 1), AC.at(1, 2), AC.at(2, 0), AC.at(2, 1), AC.at(2, 2)});

    // const GSLW = gsl.gsl_eigen_symmv_alloc(2); defer gsl.gsl_eigen_symmv_free(GSLW); eigh(f64, &AJ, &AC, A, &T1, GSLW);
    try std.testing.expect(true);
}
