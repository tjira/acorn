const std = @import("std"); const Complex = std.math.Complex;

const Matrix = @import("matrix.zig").Matrix;

pub fn Vector(comptime T: type) type {
    return struct {
        data: []T, rows: usize, allocator: std.mem.Allocator,

        pub fn init(rows: usize, allocator: std.mem.Allocator) !Vector(T) {
            return Vector(T){.data = try allocator.alloc(T, rows), .rows = rows, .allocator = allocator};
        }
        pub fn deinit(self: Vector(T)) void {
            self.allocator.free(self.data);
        }

        pub fn clone(self: Vector(T)) !Vector(T) {
            const other = try Vector(T).init(self.rows, self.allocator); @memcpy(other.data, self.data); return other;
        }
        pub fn complex(self: Vector(T)) !Vector(Complex(T)) {
            var other = try Vector(Complex(T)).init(self.rows, self.allocator); for (0..self.data.len) |i| other.data[i] = Complex(T).init(self.data[i], 0); return other;
        }

        pub fn at(self: Vector(T), i: usize) T {
            return self.data[i];
        }
        pub fn ptr(self: Vector(T), i: usize) *T {
            return &self.data[i];
        }
        pub fn matrixptr(self: Vector(T)) Matrix(T) {
            return Matrix(T){.data = self.data[0..], .rows = self.rows, .cols = 1, .allocator = self.allocator};
        }

        pub fn fill(self: Vector(T), value: T) void {
            for (0..self.data.len) |i| self.data[i] = value;
        }
    };
}

pub fn eq(comptime T: type, u: Vector(T), v: Vector(T), epsilon: T) bool {
    if (u.rows != v.rows) return false; for (0..u.data.len) |i| if (@abs(u.data[i] - v.data[i]) > epsilon) return false; return true;
}

pub fn ceq(comptime T: type, u: Vector(Complex(T)), v: Vector(Complex(T)), epsilon: T) bool {
    if (u.rows != v.rows) return false; for (0..u.data.len) |i| if (@abs(u.data[i].re - v.data[i].re) > epsilon or @abs(u.data[i].im - v.data[i].im) > epsilon) return false; return true;
}

// TESTS ===============================================================================================================================================================================================

test "vector_init_1" {
    var u = try Vector(f64).init(1, std.testing.allocator); defer u.deinit();

    try std.testing.expect(u.data.len == 1 and u.rows == 1);
}

test "vector_init_2" {
    var u = try Vector(f64).init(2, std.testing.allocator); defer u.deinit();

    try std.testing.expect(u.data.len == 2 and u.rows == 2);
}

test "vector_clone_1" {
    var u = try Vector(f64).init(1, std.testing.allocator); defer u.deinit();

    u.ptr(0).* = 3;

    var v = try u.clone(); defer v.deinit();

    try std.testing.expect(eq(f64, u, v, 1e-12));
}

test "vector_clone_2" {
    var u = try Vector(f64).init(2, std.testing.allocator); defer u.deinit();

    u.ptr(0).* = 3;
    u.ptr(1).* = 4;

    var v = try u.clone(); defer v.deinit();

    try std.testing.expect(eq(f64, u, v, 1e-12));
}

test "vector_complex_1" {
    var u = try Vector(f64         ).init(1, std.testing.allocator); defer u.deinit();
    var r = try Vector(Complex(f64)).init(1, std.testing.allocator); defer r.deinit();

    u.ptr(0).* = 3;

    r.ptr(0).* = Complex(f64).init(3, 0);

    var v = try u.complex(); defer v.deinit();

    try std.testing.expect(ceq(f64, v, r, 1e-12));
}

test "vector_complex_2" {
    var u = try Vector(f64         ).init(2, std.testing.allocator); defer u.deinit();
    var r = try Vector(Complex(f64)).init(2, std.testing.allocator); defer r.deinit();

    u.ptr(0).* = 3;
    u.ptr(1).* = 4;

    r.ptr(0).* = Complex(f64).init(3, 0);
    r.ptr(1).* = Complex(f64).init(4, 0);

    var v = try u.complex(); defer v.deinit();

    try std.testing.expect(ceq(f64, v, r, 1e-12));
}

test "vector_at_ptr_1" {
    var u = try Vector(f64).init(1, std.testing.allocator); defer u.deinit();

    u.ptr(0).* = 3;

    try std.testing.expect(u.at(0) == 3);
}

test "vector_at_ptr_2" {
    var u = try Vector(f64).init(2, std.testing.allocator); defer u.deinit();

    u.ptr(0).* = 3;
    u.ptr(1).* = 4;

    try std.testing.expect(u.at(0) == 3 and u.at(1) == 4);
}

test "vector_matrixptr_1" {
    var u = try Vector(f64).init(1, std.testing.allocator); defer u.deinit();

    u.ptr(0).* = 3;

    const v = u.matrixptr().vectorptr();

    try std.testing.expect(eq(f64, u, v, 1e-12));
}

test "vector_matrixptr_2" {
    var u = try Vector(f64).init(2, std.testing.allocator); defer u.deinit();

    u.ptr(0).* = 3;
    u.ptr(1).* = 4;

    const v = u.matrixptr().vectorptr();

    try std.testing.expect(eq(f64, u, v, 1e-12));
}

test "vector_fill_1" {
    var u = try Vector(f64).init(1, std.testing.allocator); defer u.deinit();
    var r = try Vector(f64).init(1, std.testing.allocator); defer r.deinit();

    r.ptr(0).* = 3;

    u.fill(3);

    try std.testing.expect(eq(f64, u, r, 1e-12));
}

test "vector_fill_2" {
    var u = try Vector(f64).init(2, std.testing.allocator); defer u.deinit();
    var r = try Vector(f64).init(2, std.testing.allocator); defer r.deinit();

    r.ptr(0).* = 3;
    r.ptr(1).* = 3;

    u.fill(3);

    try std.testing.expect(eq(f64, u, r, 1e-12));
}

test "eq_1" {
    var u = try Vector(f64).init(1, std.testing.allocator); defer u.deinit();
    var v = try Vector(f64).init(1, std.testing.allocator); defer v.deinit();

    u.ptr(0).* = 3;

    v.ptr(0).* = 3;

    try std.testing.expect(eq(f64, u, v, 1e-12));
}

test "eq_2" {
    var u = try Vector(f64).init(2, std.testing.allocator); defer u.deinit();
    var v = try Vector(f64).init(2, std.testing.allocator); defer v.deinit();

    u.ptr(0).* = 3;
    u.ptr(1).* = 4;

    v.ptr(0).* = 3;
    v.ptr(1).* = 4;

    try std.testing.expect(eq(f64, u, v, 1e-12));
}

test "ceq_1" {
    var u = try Vector(Complex(f64)).init(1, std.testing.allocator); defer u.deinit();
    var v = try Vector(Complex(f64)).init(1, std.testing.allocator); defer v.deinit();

    u.ptr(0).* = Complex(f64).init(3, 1);

    v.ptr(0).* = Complex(f64).init(3, 1);

    try std.testing.expect(ceq(f64, u, v, 1e-12));
}

test "ceq_2" {
    var u = try Vector(Complex(f64)).init(2, std.testing.allocator); defer u.deinit();
    var v = try Vector(Complex(f64)).init(2, std.testing.allocator); defer v.deinit();

    u.ptr(0).* = Complex(f64).init(3, 1);
    u.ptr(1).* = Complex(f64).init(4, 2);

    v.ptr(0).* = Complex(f64).init(3, 1);
    v.ptr(1).* = Complex(f64).init(4, 2);

    try std.testing.expect(ceq(f64, u, v, 1e-12));
}
