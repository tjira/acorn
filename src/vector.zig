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
            const other = try Vector(T).init(self.rows, self.allocator);

            @memcpy(other.data, self.data);

            return other;
        }
        pub fn complex(self: Vector(T)) !Vector(Complex(T)) {
            var other = try Vector(Complex(T)).init(self.rows, self.allocator);

            for (0..self.data.len) |i| other.data[i] = Complex(T).init(self.data[i], 0);

            return other;
        }

        pub fn at(self: Vector(T), i: usize) T {
            return self.ptr(i).*;
        }
        pub fn ptr(self: Vector(T), i: usize) *T {
            return &self.data[i];
        }
        pub fn matrix(self: Vector(T)) Matrix(T) {
            return Matrix(T){.data = self.data, .rows = self.rows, .cols = 1, .allocator = self.allocator};
        }

        pub fn fill(self: Vector(T), value: T) void {
            for (0..self.data.len) |i| self.data[i] = value;
        }
    };
}

pub fn eq(comptime T: type, u: Vector(T), v: Vector(T), epsilon: f64) bool {
    if (u.rows != v.rows) return false;

    if (@typeInfo(T) != .Struct) for (0..u.data.len) |i| if (@abs(u.data[i] - v.data[i]) > epsilon or std.math.isNan(@abs(u.data[i] - v.data[i]))) {
        return false;
    };

    if (@typeInfo(T) == .Struct) for (0..u.data.len) |i| if (@abs(u.data[i].re - v.data[i].re) > epsilon or @abs(u.data[i].im - v.data[i].im) > epsilon) {
        return false;
    };

    return true;
}
