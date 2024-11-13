const std = @import("std");

pub fn Vector(comptime T: type) type {
    return struct {
        data: std.ArrayList(T),
        rows: usize,

        pub fn init(m: usize, data: []const T, allocator: std.mem.Allocator) !Vector(T) {
            var vector = Vector(T){.data = std.ArrayList(T).init(allocator), .rows = m};
            try vector.data.appendSlice(data);
            return vector;
        }
        pub fn zero(m: usize, allocator: std.mem.Allocator) !Vector(T) {
            var vector = Vector(T){.data = std.ArrayList(T).init(allocator), .rows = m};
            try vector.data.appendNTimes(0, m);
            return vector;
        }

        pub fn at(self: Vector(T), i: usize) !T {
            if (i >= self.rows) return error.IndexOutOfRange;
            return self.data.items[i];
        }
        pub fn clone(self: Vector(T)) !Vector(T) {
            return Vector(T){.data = try self.data.clone(), .rows = self.rows};
        }

        pub fn add(self: Vector(T), other: Vector(T)) !Vector(T) {
            if (self.rows != other.rows) return error.IncompatibleVectors;
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] += other.data.items[i];
            return result;
        }
        pub fn addScalar(self: Vector(T), scalar: T) !Vector(T) {
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] += scalar;
            return result;
        }
        pub fn divScalar(self: Vector(T), scalar: T) !Vector(T) {
            return self.mulScalar(1 / scalar);
        }
        pub fn mulScalar(self: Vector(T), scalar: T) !Vector(T) {
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] *= scalar;
            return result;
        }
        pub fn subScalar(self: Vector(T), scalar: T) !Vector(T) {
            return self.addScalar(-scalar);
        }

        pub fn set(self: Vector(T), i: usize, value: T) !void {
            if (i >= self.rows) return error.IndexOutOfRange;
            self.data.items[i] = value;
        }
    };
}
