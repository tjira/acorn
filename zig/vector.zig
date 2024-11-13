const std = @import("std");

pub fn Vector(comptime T: type) type {
    return struct {
        data: std.ArrayList(T),
        rows: usize,

        // =============================================================================================================================================================================================

        pub fn clone(self: Vector(T)) !Vector(T) {
            return Vector(T){.data = try self.data.clone(), .rows = self.rows};
        }
        pub fn constant(m: usize, value: T, allocator: std.mem.Allocator) !Vector(T) {
            var vector = Vector(T){.data = std.ArrayList(T).init(allocator), .rows = m};
            try vector.data.appendNTimes(value, m);
            return vector;
        }
        pub fn init(m: usize, data: []const T, allocator: std.mem.Allocator) !Vector(T) {
            var vector = Vector(T){.data = std.ArrayList(T).init(allocator), .rows = m};
            try vector.data.appendSlice(data);
            return vector;
        }
        pub fn zero(m: usize, allocator: std.mem.Allocator) !Vector(T) {
            return Vector(T).constant(m, 0, allocator);
        }

        // =============================================================================================================================================================================================

        pub fn at(self: Vector(T), i: usize) T {
            return self.data.items[i];
        }

        // =============================================================================================================================================================================================

        pub fn eqApprox(self: Vector(T), other: Vector(T), tolerance: f64) bool {
            for (0..self.rows) |i| if (self.at(i) - other.at(i) > tolerance) return false;
            return true;
        }

        // =============================================================================================================================================================================================

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
        pub fn div(self: Vector(T), other: Vector(T)) !Vector(T) {
            if (self.rows != other.rows) return error.IncompatibleVectors;
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] /= other.data.items[i];
            return result;
        }
        pub fn divScalar(self: Vector(T), scalar: T) !Vector(T) {
            return self.mulScalar(1 / scalar);
        }
        pub fn mul(self: Vector(T), other: Vector(T)) !Vector(T) {
            if (self.rows != other.rows) return error.IncompatibleVectors;
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] *= other.data.items[i];
            return result;
        }
        pub fn mulScalar(self: Vector(T), scalar: T) !Vector(T) {
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] *= scalar;
            return result;
        }
        pub fn sub(self: Vector(T), other: Vector(T)) !Vector(T) {
            if (self.rows != other.rows) return error.IncompatibleVectors;
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] -= other.data.items[i];
            return result;
        }
        pub fn subScalar(self: Vector(T), scalar: T) !Vector(T) {
            return self.addScalar(-scalar);
        }

        // =============================================================================================================================================================================================


        pub fn sum(self: Vector(T)) T {
            var result: T = 0;
            for (self.data.items) |e| result += e;
            return result;
        }

        // =============================================================================================================================================================================================

        pub fn set(self: Vector(T), i: usize, value: T) void {
            self.data.items[i] = value;
        }

        // =============================================================================================================================================================================================

        pub fn deinit(self: Vector(T)) void {
            self.data.deinit();
        }
    };
}
