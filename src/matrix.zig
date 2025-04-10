//! This module provides a matrix class and functions to manipulate matrices.

const std = @import("std"); const Complex = std.math.Complex;

const mth = @import("math.zig"  );
const vec = @import("vector.zig");

const StridedArray = @import("stridedarray.zig").StridedArray;
const Vector       = @import("vector.zig"      ).Vector      ;

const asfloat = @import("helper.zig").asfloat;
const istruct = @import("helper.zig").istruct;
const uncr    = @import("helper.zig").uncr   ;

/// Matrix class.
pub fn Matrix(comptime T: type) type {
    return struct {
        data: []T, rows: usize, cols: usize, allocator: std.mem.Allocator,

        /// Initialize a matrix with a given number of rows and columns and specify an allocator. The function returns an error if the allocation fails.
        pub fn init(rows: usize, cols: usize, allocator: std.mem.Allocator) !Matrix(T) {
            return Matrix(T){
                .data = try allocator.alloc(T, rows * cols),
                .rows = rows,
                .cols = cols,
                .allocator = allocator
            };
        }

        /// Free the memory allocated for the matrix.
        pub fn deinit(self: Matrix(T)) void {
            self.allocator.free(self.data);
        }

        /// Access the element at the given row and column as a value.
        pub fn at(self: Matrix(T), i: usize, j: usize) T {
            return self.ptr(i, j).*;
        }

        /// Clone the matrix. The function returns an error if the allocation of the new matrix fails.
        pub fn clone(self: Matrix(T)) !Matrix(T) {
            const other = try Matrix(T).init(self.rows, self.cols, self.allocator);

            @memcpy(other.data, self.data);

            return other;
        }

        /// Return a column of a matrix as a strided array. No memory is allocated.
        pub fn column(self: Matrix(T), j: usize) StridedArray(T) {
            return StridedArray(T){
                .data = self.data,
                .len = self.rows,
                .stride = self.cols,
                .zero = j
            };
        }

        /// Convert the matrix to a complex matrix. The function returns an error if the allocation of the new matrix fails.
        pub fn complex(self: Matrix(T)) !Matrix(Complex(T)) {
            var other = try Matrix(Complex(T)).init(self.rows, self.cols, self.allocator);

            for (0..self.rows) |i| for (0..self.cols) |j| {
                other.ptr(i, j).* = Complex(T).init(self.at(i, j), 0);
            };

            return other;
        }

        /// Fill the matrix with a given value.
        pub fn fill(self: Matrix(T), value: T) void {
            for (0..self.data.len) |i| self.data[i] = value;
        }

        /// Fill the diagonal of the matrix with ones.
        pub fn identity(self: Matrix(T)) void {
            self.fill(0);

            for (0..mth.min(self.rows, self.cols)) |i| {
                self.ptr(i, i).* = 1;
            }
        }

        /// Fill the data of the matrix with a linearly spaced vector from start to end. Both start and end are included. If the matrix has only one element, the start value is assigned to it.
        pub fn linspace(self: Matrix(T), start: T, end: T) void {
            self.data[0] = start;

            for (1..self.data.len) |i| {
                self.data[i] = start + asfloat(T, i) * (end - start) / asfloat(T, self.rows * self.cols - 1);
            }
        }

        /// Print the matrix to the given device.
        pub fn print(self: Matrix(T), device: anytype) !void {
            var buffered = std.io.bufferedWriter(device); var writer = buffered.writer();

            try writer.print("{d} {d}\n", .{self.rows, self.cols});

            for (self.data, 1..) |e, i| {
                try writer.print("{d:20.14}{s}", .{e, if(i % self.cols == 0) "\n" else " "});
            }

            try buffered.flush();
        }

        /// Access the element at the given row and column as a mutable reference.
        pub fn ptr(self: Matrix(T), i: usize, j: usize) *T {
            return &self.data[i * self.cols + j];
        }

        /// Fill the matrix with normally distributed random numbers.
        pub fn randn(self: Matrix(T), mean: T, stdev: T, seed: usize) void {
            var rng = std.rand.DefaultPrng.init(seed); const random = rng.random();

            for (0..self.data.len) |i| {
                self.data[i] = mean + stdev * random.floatNorm(T);
            }
        }

        /// Returns a reference to the row at the given index. The row is returned as a new matrix. No memory is allocated.
        pub fn row(self: Matrix(T), i: usize) Matrix(T) {
            return Matrix(T){
                .data = self.data[i * self.cols..(i + 1) * self.cols],
                .rows = 1,
                .cols = self.cols,
                .allocator = self.allocator
            };
        }

        /// Returns the matrix in the form of a vector. No memory is allocated.
        pub fn vector(self: Matrix(T)) Vector(T) {
            return Vector(T){
                .data = self.data[0..],
                .rows = self.rows * self.cols,
                .allocator = self.allocator
            };
        }

        /// Write the matrix to a file.
        pub fn write(self: Matrix(T), path: []const u8) !void {
            const file = try std.fs.cwd().createFile(path, .{}); defer file.close();

            try self.print(file.writer());
        }
    };
}

/// Add two matrices element-wise. The output matrix is stored in the matrix C.
pub fn add(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (comptime !istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] + B.data[i];
    };

    if (comptime istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].add(B.data[i]);
    };
}

/// Add a scalar to a matrix element-wise. The output matrix is stored in the matrix C.
pub fn adds(comptime T: type, C: *Matrix(T), A: Matrix(T), s: T) void {
    if (comptime !istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] + s;
    };

    if (comptime istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].add(s);
    };
}

/// Extracts the diagonal of a matrix A. The output vector is stored in the vector D.
pub fn diag(comptime T: type, D: *Vector(T), A: Matrix(T)) void {
    for (0..A.rows) |i| D.ptr(i).* = A.at(i, i);
}

/// Compares the two matrices with tolerance. The function returns true if the matrices are equal, false otherwise.
pub fn eq(comptime T: type, A: Matrix(T), B: Matrix(T), tol: T) bool {
    if (A.rows != B.rows or A.cols != B.cols) return false;

    if (comptime !istruct(T)) for (0..A.data.len) |i| if (@abs(A.data[i] - B.data[i]) > tol) return false;

    if (comptime istruct(T)) for (0..A.data.len) |i| if (A.data[i].sub(B.data[i]).abs() > tol) return false;

    return true;
}

/// Horizontally concatenate two matrices. The output matrix is stored in the matrix C.
pub fn hjoin(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    for (0..A.rows) |i| {

        for (0..A.cols) |j| {
            C.ptr(i, j).* = A.at(i, j);
        }

        for (0..B.cols) |j| {
            C.ptr(i, A.cols + j).* = B.at(i, j);
        }
    }
}

/// Matrix multiplication of the adjoint of matrix A and matrix B. The output matrix is stored in the matrix C.
pub fn mam(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (comptime istruct(T)) C.fill(T.init(0, 0)) else C.fill(0);

    if (comptime !istruct(T)) for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.rows) |k| {
        C.ptr(i, j).* += A.at(k, i) * B.at(k, j);
    };

    if (comptime istruct(T)) for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.rows) |k| {
        C.ptr(i, j).* = C.at(i, j).add(A.at(k, i).conjugate().mul(B.at(k, j)));
    };
}

/// Matrix multiplication of matrix A and matrix B. The output matrix is stored in the matrix C.
pub fn mm(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (comptime istruct(T)) C.fill(T.init(0, 0)) else C.fill(0);

    if (comptime !istruct(T)) for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {
        C.ptr(i, j).* += A.at(i, k) * B.at(k, j);
    };

    if (comptime istruct(T)) for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {
        C.ptr(i, j).* = C.at(i, j).add(A.at(i, k).mul(B.at(k, j)));
    };
}

/// Matrix multiplication of matrix A and the adjoint of matrix B. The output matrix is stored in the matrix C.
pub fn mma(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (comptime istruct(T)) C.fill(T.init(0, 0)) else C.fill(0);

    if (comptime !istruct(T)) for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {
        C.ptr(i, j).* += A.at(i, k) * B.at(j, k);
    };

    if (comptime istruct(T)) for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {
        C.ptr(i, j).* = C.at(i, j).add(A.at(i, k).mul(B.at(j, k)));
    };
}

/// Multiply two matrices element-wise. The output matrix is stored in the matrix C.
pub fn mul(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (comptime !istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] * B.data[i];
    };

    if (comptime istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].mul(B.data[i]);
    };
}

/// Multiply a matrix by a scalar. The output matrix is stored in the matrix C.
pub fn muls(comptime T: type, C: *Matrix(T), A: Matrix(T), s: T) void {
    if (comptime !istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] * s;
    };

    if (comptime istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].mul(s);
    };
}

/// Read a matrix from a file. The function returns an error if the file cannot be opened or if the matrix cannot be read.
pub fn read(comptime T: type, path: []const u8, allocator: std.mem.Allocator) !Matrix(T) {
    const file = try std.fs.cwd().openFile(path, .{}); defer file.close(); var buffer: [16]u8 = undefined;

    var buffered = std.io.bufferedReader(file.reader()); var reader = buffered.reader();
    var stream   = std.io.fixedBufferStream(&buffer);  const writer =   stream.writer();

    stream.reset(); try reader.streamUntilDelimiter(writer, ' ',  10); const rows = try std.fmt.parseInt(usize, uncr(stream.getWritten()), 10);
    stream.reset(); try reader.streamUntilDelimiter(writer, '\n', 10); const cols = try std.fmt.parseInt(usize, uncr(stream.getWritten()), 10);

    const mat = try Matrix(T).init(rows, cols, allocator); var bytes: [100]u8 = undefined;

    for (0..rows * cols) |i| {

        var slice: []u8 = bytes[0..0];

        while (slice.len == 0) slice = try reader.readUntilDelimiter(bytes[0..], if ((i + 1) % cols == 0) '\n' else ' ');

        var j: u32 = 0; while (slice[j] == ' ') : (j += 1) {}

        mat.data[i] = try std.fmt.parseFloat(T, uncr(slice[j..]));
    }

    return mat;
}

/// Subtract two matrices element-wise. The output matrix is stored in the matrix C.
pub fn sub(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (comptime !istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] - B.data[i];
    };

    if (comptime istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].sub(B.data[i]);
    };
}

/// Transpose the matrix A. The output matrix is stored in the matrix C.
pub fn transpose(comptime T: type, C: *Matrix(T), A: Matrix(T)) void {
    for (0..A.rows) |i| for (0..A.cols) |j| {
        C.ptr(j, i).* = A.at(i, j);
    };
}
