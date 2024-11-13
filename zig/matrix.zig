const std = @import("std");

const Vector = @import("vector.zig").Vector;

pub fn Matrix(comptime T: type) type {
    return struct {
        data: std.ArrayList(T),
        rows: usize,
        cols: usize,

        // =============================================================================================================================================================================================

        pub fn clone(self: Matrix(T)) !Matrix(T) {
            return Matrix(T){.data = try self.data.clone(), .rows = self.rows, .cols = self.cols};
        }
        pub fn constant(m: usize, n: usize, value: T, allocator: std.mem.Allocator) !Matrix(T) {
            var matrix = Matrix(T){.data = std.ArrayList(T).init(allocator), .rows = m, .cols = n};
            try matrix.data.appendNTimes(value, m * n);
            return matrix;
        }
        pub fn init(m: usize, n: usize, data: []const T, allocator: std.mem.Allocator) !Matrix(T) {
            var matrix = Matrix(T){.data = std.ArrayList(T).init(allocator), .rows = m, .cols = n};
            try matrix.data.appendSlice(data);
            return matrix;
        }
        pub fn randNorm(m: usize, n: usize, mu: T, sigma: T, seed: usize, allocator: std.mem.Allocator) !Matrix(T) {
            var matrix = Matrix(T){.data = std.ArrayList(T).init(allocator), .rows = m, .cols = n};
            var prng = std.Random.DefaultPrng.init(seed); const rand = prng.random();
            for (0..m * n) |_| try matrix.data.append(mu + sigma * rand.floatNorm(T));
            return matrix;
        }
        pub fn zero(m: usize, n: usize, allocator: std.mem.Allocator) !Matrix(T) {
            return Matrix(T).constant(m, n, 0, allocator);
        }

        // =============================================================================================================================================================================================

        pub fn at(self: Matrix(T), i: usize, j: usize) T {
            return self.data.items[i * self.cols + j];
        }
        pub fn row(self: Matrix(T), i: usize) !Vector(T) {
            var vector = try Vector(T).zero(self.cols, self.data.allocator);
            for (0..self.cols) |j| vector.set(j, self.at(i, j));
            return vector;
        }

        // =============================================================================================================================================================================================

        pub fn eqApprox(self: Matrix(T), other: Matrix(T), tolerance: f64) bool {
            if (self.rows != other.rows or self.cols != other.cols) return false;
            for (0..self.rows) |i| for (0..self.cols) |j| if (self.at(i, j) - other.at(i, j) > tolerance) return false;
            return true;
        }

        // =============================================================================================================================================================================================

        pub fn add(self: Matrix(T), other: Matrix(T)) !Matrix(T) {
            if (self.rows != other.rows or self.cols != other.cols) return error.IncompatibleMatrices;
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] += other.data.items[i];
            return result;
        }
        pub fn addScalar(self: Matrix(T), scalar: T) !Matrix(T) {
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] += scalar;
            return result;
        }
        pub fn div(self: Matrix(T), other: Matrix(T)) !Matrix(T) {
            if (self.rows != other.rows or self.cols != other.cols) return error.IncompatibleMatrices;
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] /= other.data.items[i];
            return result;
        }
        pub fn divScalar(self: Matrix(T), scalar: T) !Matrix(T) {
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] /= scalar;
            return result;
        }
        pub fn mul(self: Matrix(T), other: Matrix(T)) !Matrix(T) {
            if (self.rows != other.rows or self.cols != other.cols) return error.IncompatibleMatrices;
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] *= other.data.items[i];
            return result;
        }
        pub fn mulScalar(self: Matrix(T), scalar: T) !Matrix(T) {
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] *= scalar;
            return result;
        }
        pub fn sub(self: Matrix(T), other: Matrix(T)) !Matrix(T) {
            if (self.rows != other.rows or self.cols != other.cols) return error.IncompatibleMatrices;
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] -= other.data.items[i];
            return result;
        }
        pub fn subScalar(self: Matrix(T), scalar: T) !Matrix(T) {
            var result = try self.clone();
            for (0..self.data.items.len) |i| result.data.items[i] -= scalar;
            return result;
        }

        // =============================================================================================================================================================================================

        pub fn set(self: Matrix(T), i: usize, j: usize, value: T) void {
            self.data.items[i * self.cols + j] = value;
        }

        // =============================================================================================================================================================================================

        pub fn deinit(self: Matrix(T)) void {
            self.data.deinit();
        }
        pub fn write(self: Matrix(T), path: []const u8) !void {
            const file = try std.fs.cwd().createFile(path, .{}); defer file.close();
            _ = try file.writeAll(try std.fmt.allocPrint(self.data.allocator, "{} {}\n", .{self.rows, self.cols}));
            for (0..self.rows) |i| {
                for (0..self.cols) |j| {
                    _ = try file.writeAll(try std.fmt.allocPrint(self.data.allocator, "{d:24.16}{s}", .{self.at(i, j), if (j + 1 < self.cols) " " else "\n"}));
                }
            }
        }
    };
}
