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
        pub fn reshaped(self: Matrix(T), rows: usize, cols: usize) Matrix(T) {
            return Matrix(T){.data = self.data, .rows = rows, .cols = cols, .allocator = self.allocator};
        }
        pub fn row(self: Matrix(T), i: usize) Matrix(T) {
            return Matrix(T){.data = self.data[i * self.cols..(i + 1) * self.cols], .rows = 1, .cols = self.cols, .allocator = self.allocator};
        }
        pub fn vector(self: Matrix(T)) Vector(T) {
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
    if (A.rows != B.rows or A.cols != B.cols) return false;

    for (0..A.data.len) |i| if (@abs(A.data[i] - B.data[i]) > epsilon) return false;

    return true;
}

pub fn ceq(comptime T: type, A: Matrix(Complex(T)), B: Matrix(Complex(T)), epsilon: T) bool {
    if (A.rows != B.rows or A.cols != B.cols) return false;

    for (0..A.data.len) |i| if (@abs(A.data[i].re - B.data[i].re) > epsilon or @abs(A.data[i].im - B.data[i].im) > epsilon) return false;

    return true;
}

pub fn mm(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    C.fill(0);

    for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {
        C.ptr(i, j).* += A.at(i, k) * B.at(k, j);
    };
}

pub fn cmm(comptime T: type, C: *Matrix(Complex(T)), A: Matrix(Complex(T)), B: Matrix(Complex(T))) void {
    C.fill(Complex(T).init(0, 0));

    for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {
        C.ptr(i, j).* = C.at(i, j).add(A.at(i, k).mul(B.at(k, j)));
    };
}

pub fn mam(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    C.fill(0);

    for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.rows) |k| {
        C.ptr(i, j).* += A.at(k, i) * B.at(k, j);
    };
}

pub fn cmam(comptime T: type, C: *Matrix(Complex(T)), A: Matrix(Complex(T)), B: Matrix(Complex(T))) void {
    C.fill(Complex(T).init(0, 0));

    for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.rows) |k| {
        C.ptr(i, j).* = C.at(i, j).add(A.at(k, i).conjugate().mul(B.at(k, j)));
    };
}

pub fn mma(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    C.fill(0);

    for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {
        C.ptr(i, j).* += A.at(i, k) * B.at(j, k);
    };
}

pub fn cmma(comptime T: type, C: *Matrix(Complex(T)), A: Matrix(Complex(T)), B: Matrix(Complex(T))) void {
    C.fill(Complex(T).init(0, 0));

    for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {
        C.ptr(i, j).* = C.at(i, j).add(A.at(i, k).mul(B.at(j, k).conjugate()));
    };
}

pub fn hjoin(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    for (0..A.rows) |i| {
        for (0..A.cols) |j| C.ptr(i, j).*          = A.at(i, j);
        for (0..B.cols) |j| C.ptr(i, A.cols + j).* = B.at(i, j);
    }
}
