const std = @import("std"); const Complex = std.math.Complex;

pub fn Vector(comptime T: type) type {
    return struct {
        data: []T, rows: usize, allocator: std.mem.Allocator,

        pub fn init(rows: usize, allocator: std.mem.Allocator) !Vector(T) {
            return Vector(T){.data = try allocator.alloc(T, rows), .rows = rows, .allocator = allocator};
        }
        pub fn deinit(self: Vector(T)) void {
            self.allocator.free(self.data);
        }

        pub fn eq(self: Vector(T), other: Vector(T), epsilon: f64) bool {
            if (self.rows != other.rows) return false; for (0..self.data.len) |i| if (@abs(self.data[i] - other.data[i]) > epsilon) return false; return true;
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

        pub fn fill(self: Vector(T), value: T) void {
            for (0..self.data.len) |i| self.data[i] = value;
        }
    };
}
