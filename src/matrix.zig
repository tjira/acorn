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

pub fn eigh(comptime T: type, J: *Matrix(T), C: *Matrix(T), A: Matrix(T), T1: *Matrix(T), T2: *Matrix(T)) void {
    var maxi: usize = 0; var maxj: usize = 1; var maxv: T = 0; var phi: T = undefined; @memcpy(J.data, A.data); C.identity();

    for (0..A.rows) |i| for (i + 1..A.cols) |j| if (@abs(J.at(i, j)) > @abs(maxv)) {maxi = i; maxj = j; maxv = J.at(i, j);};

    while (@abs(maxv) > 1e-14) {

        phi = 0.5 * std.math.atan(2 * maxv / (J.at(maxi, maxi) - J.at(maxj, maxj))); T1.identity();

        T1.ptr(maxi, maxi).* = std.math.cos(phi); T1.ptr(maxj, maxj).* =  T1.at(maxi, maxi);
        T1.ptr(maxj, maxi).* = std.math.sin(phi); T1.ptr(maxi, maxj).* = -T1.at(maxj, maxi);

        mm(T, T2, J.*, T1.*); mam(T, J, T1.*, T2.*); mm(T, T2, C.*, T1.*); @memcpy(C.data, T2.data);

        maxv = 0; for (0..A.rows) |i| for (i + 1..A.cols) |j| if (@abs(J.at(i, j)) > @abs(maxv)) {maxi = i; maxj = j; maxv = J.at(i, j);};
    }

    for (0..A.rows) |i| for (i + 1..A.cols) |j| if (J.at(i, i) > J.at(j, j)) {
        std.mem.swap(T, J.ptr(i, i), J.ptr(j, j)); for (0..A.rows) |k| std.mem.swap(T, C.ptr(k, i), C.ptr(k, j));
    };
}

pub fn eq(comptime T: type, A: Matrix(T), B: Matrix(T), epsilon: T) bool {
    if (A.rows != B.rows or A.cols != B.cols) return false; for (0..A.data.len) |i| if (@abs(A.data[i] - B.data[i]) > epsilon) return false; return true;
}

pub fn ceq(comptime T: type, A: Matrix(Complex(T)), B: Matrix(Complex(T)), epsilon: T) bool {
    if (A.rows != B.rows or A.cols != B.cols) return false; for (0..A.data.len) |i| if (@abs(A.data[i].re - B.data[i].re) > epsilon or @abs(A.data[i].im - B.data[i].im) > epsilon) return false; return true;
}

pub fn mm(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    C.fill(0); for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {C.ptr(i, j).* += A.at(i, k) * B.at(k, j);};
}

pub fn cmm(comptime T: type, C: *Matrix(Complex(T)), A: Matrix(Complex(T)), B: Matrix(Complex(T))) void {
    C.fill(Complex(T).init(0, 0)); for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {C.ptr(i, j).* = C.at(i, j).add(A.at(i, k).mul(B.at(k, j)));};
}

pub fn mam(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    C.fill(0); for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.rows) |k| {C.ptr(i, j).* += A.at(k, i) * B.at(k, j);};
}

pub fn cmam(comptime T: type, C: *Matrix(Complex(T)), A: Matrix(Complex(T)), B: Matrix(Complex(T))) void {
    C.fill(Complex(T).init(0, 0)); for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.rows) |k| {C.ptr(i, j).* = C.at(i, j).add(A.at(k, i).conjugate().mul(B.at(k, j)));};
}

pub fn mma(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    C.fill(0); for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {C.ptr(i, j).* += A.at(i, k) * B.at(j, k);};
}

pub fn cmma(comptime T: type, C: *Matrix(Complex(T)), A: Matrix(Complex(T)), B: Matrix(Complex(T))) void {
    C.fill(Complex(T).init(0, 0)); for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {C.ptr(i, j).* = C.at(i, j).add(A.at(i, k).mul(B.at(j, k).conjugate()));};
}

pub fn tmamt(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    C.fill(0); for (0..C.cols) |i| for (0..C.rows) |j| for (0..A.rows) |k| {C.ptr(j, i).* += A.at(k, i) * B.at(j, k);};
}

pub fn ctmamt(comptime T: type, C: *Matrix(Complex(T)), A: Matrix(Complex(T)), B: Matrix(Complex(T))) void {
    C.fill(Complex(T).init(0, 0)); for (0..C.cols) |i| for (0..C.rows) |j| for (0..A.rows) |k| {C.ptr(j, i).* = C.at(j, i).add(A.at(k, i).conjugate().mul(B.at(j, k)));};
}

pub fn hjoin(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    for (0..A.rows) |i| {for (0..A.cols) |j| C.ptr(i, j).* = A.at(i, j); for (0..B.cols) |j| C.ptr(i, A.cols + j).* = B.at(i, j);}
}

// TESTS ===============================================================================================================================================================================================

test "matrix_init_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();

    try std.testing.expect(A.rows == 1 and A.cols == 1 and A.data.len == 1);
}

test "matrix_init_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();

    try std.testing.expect(A.rows == 2 and A.cols == 2 and A.data.len == 4);
}

test "matrix_clone_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();

    var B = try A.clone(); defer B.deinit();

    try std.testing.expect(eq(f64, A, B, 1e-12));
}

test "matrix_clone_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();

    var B = try A.clone(); defer B.deinit();

    try std.testing.expect(eq(f64, A, B, 1e-12));
}

test "matrix_complex_1x1" {
    var A = try Matrix(f64         ).init(1, 1, std.testing.allocator); defer A.deinit();
    var R = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1;

    R.ptr(0, 0).* = Complex(f64).init(1, 0);

    var B = try A.complex(); defer B.deinit();

    try std.testing.expect(ceq(f64, B, R, 1e-12));
}

test "matrix_complex_2x2" {
    var A = try Matrix(f64         ).init(2, 2, std.testing.allocator); defer A.deinit();
    var R = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    R.ptr(0, 0).* = Complex(f64).init(1, 0); R.ptr(0, 1).* = Complex(f64).init(2, 0);
    R.ptr(1, 0).* = Complex(f64).init(3, 0); R.ptr(1, 1).* = Complex(f64).init(4, 0);

    var B = try A.complex(); defer B.deinit();

    try std.testing.expect(ceq(f64, B, R, 1e-12));
}

test "matrix_at_ptr_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();

    A.ptr(0, 0).* = 1;

    try std.testing.expect(A.at(0, 0) == 1);
}

test "matrix_at_ptr_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    try std.testing.expect(A.at(0, 0) == 1 and A.at(0, 1) == 2 and A.at(1, 0) == 3 and A.at(1, 1) == 4);
}

test "matrix_rowptr_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1;

    R.ptr(0, 0).* = 1;

    const B = A.rowptr(0);

    try std.testing.expect(eq(f64, B, R, 1e-12));
}

test "matrix_rowptr_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(1, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    R.ptr(0, 0).* = 1; R.ptr(0, 1).* = 2;

    const B = A.rowptr(0);

    try std.testing.expect(eq(f64, B, R, 1e-12));
}

test "matrix_vectorptr_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1;

    R.ptr(0, 0).* = 1;

    const B = A.vectorptr().matrixptr();

    try std.testing.expect(eq(f64, B, R, 1e-12));
}

test "matrix_vectorptr_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(4, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    R.ptr(0, 0).* = 1; R.ptr(1, 0).* = 2; R.ptr(2, 0).* = 3; R.ptr(3, 0).* = 4;

    const B = A.vectorptr().matrixptr();

    try std.testing.expect(eq(f64, B, R, 1e-12));
}

test "matrix_fill_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.fill(1);

    R.ptr(0, 0).* = 1;

    try std.testing.expect(eq(f64, A, R, 1e-12));
}

test "matrix_fill_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.fill(1);

    R.ptr(0, 0).* = 1; R.ptr(0, 1).* = 1;
    R.ptr(1, 0).* = 1; R.ptr(1, 1).* = 1;

    try std.testing.expect(eq(f64, A, R, 1e-12));
}

test "matrix_identity_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.identity();

    R.ptr(0, 0).* = 1;

    try std.testing.expect(eq(f64, A, R, 1e-12));
}

test "matrix_identity_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.identity();

    R.ptr(0, 0).* = 1; R.ptr(0, 1).* = 0;
    R.ptr(1, 0).* = 0; R.ptr(1, 1).* = 1;

    try std.testing.expect(eq(f64, A, R, 1e-12));
}

test "matrix_linspace_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.linspace(1, 1);

    R.ptr(0, 0).* = 1;

    try std.testing.expect(eq(f64, A, R, 1e-12));
}

test "matrix_linspace_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.linspace(1, 4);

    R.ptr(0, 0).* = 1; R.ptr(0, 1).* = 2;
    R.ptr(1, 0).* = 3; R.ptr(1, 1).* = 4;

    try std.testing.expect(eq(f64, A, R, 1e-12));
}

test "eigh_1x1" {
    var T1 = try Matrix(f64).init(1, 1, std.testing.allocator); defer T1.deinit();
    var T2 = try Matrix(f64).init(1, 1, std.testing.allocator); defer T2.deinit();

    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var J = try Matrix(f64).init(1, 1, std.testing.allocator); defer J.deinit();
    var C = try Matrix(f64).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();
    var S = try Matrix(f64).init(1, 1, std.testing.allocator); defer S.deinit();

    A.ptr(0, 0).* = 1;

    R.ptr(0, 0).* = 1;

    S.ptr(0, 0).* = 1;

    eigh(f64, &J, &C, A, &T1, &T2);

    try std.testing.expect(eq(f64, J, R, 1e-12) and eq(f64, C, S, 1e-12));
}

test "eigh_2x2" {
    var T1 = try Matrix(f64).init(2, 2, std.testing.allocator); defer T1.deinit();
    var T2 = try Matrix(f64).init(2, 2, std.testing.allocator); defer T2.deinit();

    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var J = try Matrix(f64).init(2, 2, std.testing.allocator); defer J.deinit();
    var C = try Matrix(f64).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();
    var S = try Matrix(f64).init(2, 2, std.testing.allocator); defer S.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 2; A.ptr(1, 1).* = 1;

    R.ptr(0, 0).* = -1; R.ptr(1, 1).* = 3;

    S.ptr(0, 0).* = -std.math.sqrt1_2; S.ptr(0, 1).* = std.math.sqrt1_2;
    S.ptr(1, 0).* =  std.math.sqrt1_2; S.ptr(1, 1).* = std.math.sqrt1_2;

    eigh(f64, &J, &C, A, &T1, &T2);

    try std.testing.expect(eq(f64, J, R, 1e-12) and eq(f64, C, S, 1e-12));
}

test "eq_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(1, 1, std.testing.allocator); defer B.deinit();

    A.ptr(0, 0).* = 1;

    B.ptr(0, 0).* = 1;

    try std.testing.expect(eq(f64, A, B, 1e-12));
}

test "eq_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    B.ptr(0, 0).* = 1; B.ptr(0, 1).* = 2;
    B.ptr(1, 0).* = 3; B.ptr(1, 1).* = 4;

    try std.testing.expect(eq(f64, A, B, 1e-12));
}

test "ceq_1x1" {
    var A = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer B.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2);

    B.ptr(0, 0).* = Complex(f64).init(1, 2);

    try std.testing.expect(ceq(f64, A, B, 1e-12));
}

test "ceq_2x2" {
    var A = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer B.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2); A.ptr(0, 1).* = Complex(f64).init(3, 4);
    A.ptr(1, 0).* = Complex(f64).init(5, 6); A.ptr(1, 1).* = Complex(f64).init(7, 8);

    B.ptr(0, 0).* = Complex(f64).init(1, 2); B.ptr(0, 1).* = Complex(f64).init(3, 4);
    B.ptr(1, 0).* = Complex(f64).init(5, 6); B.ptr(1, 1).* = Complex(f64).init(7, 8);

    try std.testing.expect(ceq(f64, A, B, 1e-12));
}

test "hjoin_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(1, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(1, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1;

    B.ptr(0, 0).* = 2;

    R.ptr(0, 0).* = 1; R.ptr(0, 1).* = 2;

    hjoin(f64, &C, A, B);

    try std.testing.expect(eq(f64, C, R, 1e-12));
}

test "hjoin_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(2, 4, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(2, 4, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    B.ptr(0, 0).* = 5; B.ptr(0, 1).* = 6;
    B.ptr(1, 0).* = 7; B.ptr(1, 1).* = 8;

    R.ptr(0, 0).* = 1; R.ptr(0, 1).* = 2; R.ptr(0, 2).* = 5; R.ptr(0, 3).* = 6;
    R.ptr(1, 0).* = 3; R.ptr(1, 1).* = 4; R.ptr(1, 2).* = 7; R.ptr(1, 3).* = 8;

    hjoin(f64, &C, A, B);

    try std.testing.expect(eq(f64, C, R, 1e-12));
}

test "mm_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 2;

    B.ptr(0, 0).* = 2;

    R.ptr(0, 0).* = 4;

    mm(f64, &C, A, B);

    try std.testing.expect(eq(f64, C, R, 1e-12));
}

test "mm_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    B.ptr(0, 0).* = 5; B.ptr(0, 1).* = 6;
    B.ptr(1, 0).* = 7; B.ptr(1, 1).* = 8;

    R.ptr(0, 0).* = 19; R.ptr(0, 1).* = 22;
    R.ptr(1, 0).* = 43; R.ptr(1, 1).* = 50;

    mm(f64, &C, A, B);

    try std.testing.expect(eq(f64, C, R, 1e-12));
}

test "cmm_1x1" {
    var A = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2);

    B.ptr(0, 0).* = Complex(f64).init(3, 4);

    R.ptr(0, 0).* = Complex(f64).init(-5, 10);

    cmm(f64, &C, A, B);

    try std.testing.expect(ceq(f64, C, R, 1e-12));
}

test "cmm_2x2" {
    var A = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2); A.ptr(0, 1).* = Complex(f64).init(3, 4);
    A.ptr(1, 0).* = Complex(f64).init(5, 6); A.ptr(1, 1).* = Complex(f64).init(7, 8);

    B.ptr(0, 0).* = Complex(f64).init( 9, 10); B.ptr(0, 1).* = Complex(f64).init(11, 12);
    B.ptr(1, 0).* = Complex(f64).init(13, 14); B.ptr(1, 1).* = Complex(f64).init(15, 16);

    R.ptr(0, 0).* = Complex(f64).init(-28, 122); R.ptr(0, 1).* = Complex(f64).init(-32, 142);
    R.ptr(1, 0).* = Complex(f64).init(-36, 306); R.ptr(1, 1).* = Complex(f64).init(-40, 358);

    cmm(f64, &C, A, B);

    try std.testing.expect(ceq(f64, C, R, 1e-12));
}

test "mam_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 2;

    B.ptr(0, 0).* = 2;

    R.ptr(0, 0).* = 4;

    mam(f64, &C, A, B);

    try std.testing.expect(eq(f64, C, R, 1e-12));
}

test "mam_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    B.ptr(0, 0).* = 5; B.ptr(0, 1).* = 6;
    B.ptr(1, 0).* = 7; B.ptr(1, 1).* = 8;

    R.ptr(0, 0).* = 26; R.ptr(0, 1).* = 30;
    R.ptr(1, 0).* = 38; R.ptr(1, 1).* = 44;

    mam(f64, &C, A, B);

    try std.testing.expect(eq(f64, C, R, 1e-12));
}

test "cmam_1x1" {
    var A = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2);

    B.ptr(0, 0).* = Complex(f64).init(3, 4);

    R.ptr(0, 0).* = Complex(f64).init(11, -2);

    cmam(f64, &C, A, B);

    try std.testing.expect(ceq(f64, C, R, 1e-12));
}

test "cmam_2x2" {
    var A = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2); A.ptr(0, 1).* = Complex(f64).init(3, 4);
    A.ptr(1, 0).* = Complex(f64).init(5, 6); A.ptr(1, 1).* = Complex(f64).init(7, 8);

    B.ptr(0, 0).* = Complex(f64).init( 9, 10); B.ptr(0, 1).* = Complex(f64).init(11, 12);
    B.ptr(1, 0).* = Complex(f64).init(13, 14); B.ptr(1, 1).* = Complex(f64).init(15, 16);

    R.ptr(0, 0).* = Complex(f64).init(178, -16); R.ptr(0, 1).* = Complex(f64).init(206, -20);
    R.ptr(1, 0).* = Complex(f64).init(270, -12); R.ptr(1, 1).* = Complex(f64).init(314, -16);

    cmam(f64, &C, A, B);

    try std.testing.expect(ceq(f64, C, R, 1e-12));
}

test "mma_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 2;

    B.ptr(0, 0).* = 2;

    R.ptr(0, 0).* = 4;

    mma(f64, &C, A, B);

    try std.testing.expect(eq(f64, C, R, 1e-12));
}

test "mma_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    B.ptr(0, 0).* = 5; B.ptr(0, 1).* = 6;
    B.ptr(1, 0).* = 7; B.ptr(1, 1).* = 8;

    R.ptr(0, 0).* = 17; R.ptr(0, 1).* = 23;
    R.ptr(1, 0).* = 39; R.ptr(1, 1).* = 53;

    mma(f64, &C, A, B);

    try std.testing.expect(eq(f64, C, R, 1e-12));
}

test "cmma_1x1" {
    var A = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2);

    B.ptr(0, 0).* = Complex(f64).init(3, 4);

    R.ptr(0, 0).* = Complex(f64).init(11, 2);

    cmma(f64, &C, A, B);

    try std.testing.expect(ceq(f64, C, R, 1e-12));
}

test "cmma_2x2" {
    var A = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2); A.ptr(0, 1).* = Complex(f64).init(3, 4);
    A.ptr(1, 0).* = Complex(f64).init(5, 6); A.ptr(1, 1).* = Complex(f64).init(7, 8);

    B.ptr(0, 0).* = Complex(f64).init( 9, 10); B.ptr(0, 1).* = Complex(f64).init(11, 12);
    B.ptr(1, 0).* = Complex(f64).init(13, 14); B.ptr(1, 1).* = Complex(f64).init(15, 16);

    R.ptr(0, 0).* = Complex(f64).init(110, 16); R.ptr(0, 1).* = Complex(f64).init(150, 24);
    R.ptr(1, 0).* = Complex(f64).init(278,  8); R.ptr(1, 1).* = Complex(f64).init(382, 16);

    cmma(f64, &C, A, B);

    try std.testing.expect(ceq(f64, C, R, 1e-12));
}

test "tmamt_1x1" {
    var A = try Matrix(f64).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 2;

    B.ptr(0, 0).* = 2;

    R.ptr(0, 0).* = 4;

    tmamt(f64, &C, A, B);

    try std.testing.expect(eq(f64, C, R, 1e-12));
}

test "tmamt_2x2" {
    var A = try Matrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(f64).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(f64).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2;
    A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;

    B.ptr(0, 0).* = 5; B.ptr(0, 1).* = 6;
    B.ptr(1, 0).* = 7; B.ptr(1, 1).* = 8;

    R.ptr(0, 0).* = 23; R.ptr(0, 1).* = 34;
    R.ptr(1, 0).* = 31; R.ptr(1, 1).* = 46;

    tmamt(f64, &C, A, B);

    try std.testing.expect(eq(f64, C, R, 1e-12));
}

test "ctmamt_1x1" {
    var A = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(1, 1, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2);

    B.ptr(0, 0).* = Complex(f64).init(3, 4);

    R.ptr(0, 0).* = Complex(f64).init(11, -2);

    ctmamt(f64, &C, A, B);

    try std.testing.expect(ceq(f64, C, R, 1e-12));
}

test "ctmamt_2x2" {
    var A = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer C.deinit();
    var R = try Matrix(Complex(f64)).init(2, 2, std.testing.allocator); defer R.deinit();

    A.ptr(0, 0).* = Complex(f64).init(1, 2); A.ptr(0, 1).* = Complex(f64).init(3, 4);
    A.ptr(1, 0).* = Complex(f64).init(5, 6); A.ptr(1, 1).* = Complex(f64).init(7, 8);

    B.ptr(0, 0).* = Complex(f64).init( 9, 10); B.ptr(0, 1).* = Complex(f64).init(11, 12);
    B.ptr(1, 0).* = Complex(f64).init(13, 14); B.ptr(1, 1).* = Complex(f64).init(15, 16);

    R.ptr(0, 0).* = Complex(f64).init(156, -14); R.ptr(0, 1).* = Complex(f64).init(240, -10);
    R.ptr(1, 0).* = Complex(f64).init(212, -22); R.ptr(1, 1).* = Complex(f64).init(328, -18);

    ctmamt(f64, &C, A, B);

    try std.testing.expect(ceq(f64, C, R, 1e-12));
}
