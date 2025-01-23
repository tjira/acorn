//! This module provides a matrix class and functions to manipulate matrices.

const std = @import("std"); const Complex = std.math.Complex;

const StridedArray = @import("stridedarray.zig").StridedArray;
const Vector       = @import("vector.zig"      ).Vector      ;

const asfloat = @import("helper.zig").asfloat;
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

        /// Clone the matrix. The function returns an error if the allocation of the new matrix fails.
        pub fn clone(self: Matrix(T)) !Matrix(T) {
            const other = try Matrix(T).init(self.rows, self.cols, self.allocator);

            @memcpy(other.data, self.data);

            return other;
        }

        /// Convert the matrix to a complex matrix. The function returns an error if the allocation of the new matrix fails.
        pub fn complex(self: Matrix(T)) !Matrix(Complex(T)) {
            var other = try Matrix(Complex(T)).init(self.rows, self.cols, self.allocator);

            for (0..self.data.len) |i| {
                other.data[i] = Complex(T).init(self.data[i], 0);
            }

            return other;
        }

        /// Access the element at the given row and column as a value.
        pub fn at(self: Matrix(T), i: usize, j: usize) T {
            return self.ptr(i, j).*;
        }

        /// Access the element at the given row and column as a mutable reference.
        pub fn ptr(self: Matrix(T), i: usize, j: usize) *T {
            return &self.data[i * self.cols + j];
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

        /// Returns the matrix in the form of a strided array. No memory is allocated.
        pub fn sa(self: Matrix(T)) StridedArray(T) {
            return StridedArray(T){.data = self.data, .len = self.rows * self.cols, .stride = 1, .zero = 0};
        }

        /// Returns the matrix in the form of a vector. No memory is allocated.
        pub fn vector(self: Matrix(T)) Vector(T) {
            return Vector(T){
                .data = self.data[0..],
                .rows = self.rows * self.cols,
                .allocator = self.allocator
            };
        }

        /// Fill the matrix with a given value.
        pub fn fill(self: Matrix(T), value: T) void {
            for (0..self.data.len) |i| self.data[i] = value;
        }

        /// Fill the diagonal of the matrix with ones.
        pub fn identity(self: Matrix(T)) void {
            self.fill(0);

            for (0..self.rows) |i| {
                self.ptr(i, i).* = 1;
            }
        }

        /// Fill the data of the matrix with a linearly spaced vector from start to end. Both start and end are included.
        pub fn linspace(self: Matrix(T), start: T, end: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = start + asfloat(T, i) * (end - start) / asfloat(T, self.rows * self.cols - 1);
            }
        }

        /// Print the matrix to the given device.
        pub fn print(self: Matrix(T), device: anytype) !void {
            try device.print("{d} {d}\n", .{self.rows, self.cols});

            for (self.data, 1..) |e, i| {
                try device.print("{d:20.14}{s}", .{e, if(i % self.cols == 0) "\n" else " "});
            }
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
    if (@typeInfo(T) != .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] + B.data[i];
    };

    if (@typeInfo(T) == .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].add(B.data[i]);
    };
}

/// Add a scalar to a matrix element-wise. The output matrix is stored in the matrix C.
pub fn adds(comptime T: type, C: *Matrix(T), A: Matrix(T), s: T) void {
    if (@typeInfo(T) != .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] + s;
    };

    if (@typeInfo(T) == .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].add(s);
    };
}

/// Find the eigenvalues and eigenvectors of a real symmetric matrix A. The eigenvalues are stored in the diagonal of the matrix J, and the eigenvectors are stored in the columns of the matrix C. The matrices T1 and T2 are temporary matrices used in the computation.
pub fn eigh(comptime T: type, J: *Matrix(T), C: *Matrix(T), A: Matrix(T), T1: *Matrix(T), T2: *Matrix(T)) void {
    var maxi: usize = 0; var maxj: usize = 1; var maxv: T = 0; var phi: T = undefined; @memcpy(J.data, A.data); C.identity();

    for (0..A.rows) |i| for (i + 1..A.cols) |j| if (@abs(J.at(i, j)) > @abs(maxv)) {
        maxi = i; maxj = j; maxv = J.at(i, j);
    };

    while (@abs(maxv) > 1e-14) {

        phi = 0.5 * std.math.atan(2 * maxv / (J.at(maxi, maxi) - J.at(maxj, maxj))); T1.identity();

        T1.ptr(maxi, maxi).* = std.math.cos(phi); T1.ptr(maxj, maxj).* =  T1.at(maxi, maxi);
        T1.ptr(maxj, maxi).* = std.math.sin(phi); T1.ptr(maxi, maxj).* = -T1.at(maxj, maxi);

        mm(T, T2, J.*, T1.*); mam(T, J, T1.*, T2.*); mm(T, T2, C.*, T1.*); @memcpy(C.data, T2.data);

        maxv = 0;

        for (0..A.rows) |i| for (i + 1..A.cols) |j| if (@abs(J.at(i, j)) > @abs(maxv)) {
            maxi = i; maxj = j; maxv = J.at(i, j);
        };
    }

    for (0..A.rows) |i| for (i + 1..A.cols) |j| if (J.at(i, i) > J.at(j, j)) {

        std.mem.swap(T, J.ptr(i, i), J.ptr(j, j));

        for (0..A.rows) |k| {
            std.mem.swap(T, C.ptr(k, i), C.ptr(k, j));
        }
    };
}

/// Check if two matrices are equal within a given tolerance. The function returns true if the matrices are equal and false otherwise.
pub fn eq(comptime T: type, A: Matrix(T), B: Matrix(T), epsilon: T) bool {
    if (A.rows != B.rows or A.cols != B.cols) return false;

    if (@typeInfo(T) != .Struct) for (0..A.data.len) |i| if (@abs(A.data[i] - B.data[i]) > epsilon or std.math.isNan(@abs(A.data[i] - B.data[i]))) {
        return false;
    };

    if (@typeInfo(T) == .Struct) for (0..A.data.len) |i| if (@abs(A.data[i].re - B.data[i].re) > epsilon or @abs(A.data[i].im - B.data[i].im) > epsilon) {
        return false;
    };

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

/// Solve the linear system Ax = b. The output vector is stored in the vector x. The matrix A is converted to a row-echelon form. Vector b is also modified.
pub fn linsolve(comptime T: type, x: *Vector(T), A: *Matrix(T), b: *Vector(T)) void {
    for (0..A.cols - 1) |j| for (j + 1..A.rows) |i| {

        const factor = A.at(i, j) / A.at(j, j);

        for (j..A.cols) |k| {
            A.ptr(i, k).* -= factor * A.at(j, k);
        }

        b.ptr(i).* -= factor * b.at(j);
    };

    for (0..x.rows) |i| {

        for (0..i) |j| {
            b.ptr(b.rows - i - 1).* -= A.at(A.rows - i - 1, A.cols - j - 1) * x.at(x.rows - j - 1);
        }

        x.ptr(x.rows - i - 1).* = b.at(x.rows - i - 1) / A.at(x.rows - i - 1, x.rows - i - 1);
    }
}

/// Matrix multiplication of the adjoint of matrix A and matrix B. The output matrix is stored in the matrix C.
pub fn mam(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (@typeInfo(T) == .Struct) C.fill(T.init(0, 0)) else C.fill(0);

    if (@typeInfo(T) != .Struct) for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.rows) |k| {
        C.ptr(i, j).* += A.at(k, i) * B.at(k, j);
    };

    if (@typeInfo(T) == .Struct) for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.rows) |k| {
        C.ptr(i, j).* = C.at(i, j).add(A.at(k, i).conjugate().mul(B.at(k, j)));
    };
}

/// Matrix multiplication of matrix A and matrix B. The output matrix is stored in the matrix C.
pub fn mm(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (@typeInfo(T) == .Struct) C.fill(T.init(0, 0)) else C.fill(0);

    if (@typeInfo(T) != .Struct) for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {
        C.ptr(i, j).* += A.at(i, k) * B.at(k, j);
    };

    if (@typeInfo(T) == .Struct) for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {
        C.ptr(i, j).* = C.at(i, j).add(A.at(i, k).mul(B.at(k, j)));
    };
}

/// Matrix multiplication of matrix A and the adjoint of matrix B. The output matrix is stored in the matrix C.
pub fn mma(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (@typeInfo(T) == .Struct) C.fill(T.init(0, 0)) else C.fill(0);

    if (@typeInfo(T) != .Struct) for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {
        C.ptr(i, j).* += A.at(i, k) * B.at(j, k);
    };

    if (@typeInfo(T) == .Struct) for (0..C.rows) |i| for (0..C.cols) |j| for (0..A.cols) |k| {
        C.ptr(i, j).* = C.at(i, j).add(A.at(i, k).mul(B.at(j, k)));
    };
}

/// Multiply two matrices element-wise. The output matrix is stored in the matrix C.
pub fn mul(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (@typeInfo(T) != .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] * B.data[i];
    };

    if (@typeInfo(T) == .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].mul(B.data[i]);
    };
}

/// Multiply a matrix by a scalar element-wise. The output matrix is stored in the matrix C.
pub fn muls(comptime T: type, C: *Matrix(T), A: Matrix(T), s: T) void {
    if (@typeInfo(T) != .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] * s;
    };

    if (@typeInfo(T) == .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].mul(s);
    };
}

/// Read a matrix from a file. The function returns an error if the file cannot be opened or if the matrix cannot be read.
pub fn read(comptime T: type, path: []const u8, allocator: std.mem.Allocator) !Matrix(T) {
    const file = try std.fs.cwd().openFile(path, .{}); defer file.close(); var buffer: [16]u8 = undefined;

    var buffered = std.io.bufferedReader(file.reader()); var reader = buffered.reader();
    var stream   = std.io.fixedBufferStream(&buffer);  const writer =   stream.writer();

    stream.reset(); try reader.streamUntilDelimiter(writer, ' ',  4); const rows = try std.fmt.parseInt(usize, uncr(stream.getWritten()), 10);
    stream.reset(); try reader.streamUntilDelimiter(writer, '\n', 4); const cols = try std.fmt.parseInt(usize, uncr(stream.getWritten()), 10);

    const mat = try Matrix(T).init(rows, cols, allocator);

    for (0..rows * cols) |i| {

        const bytes = try reader.readBytesNoEof(20); try reader.skipBytes(1, .{});

        var j: u32 = 0; while (bytes[j] == ' ') : (j += 1) {}

        mat.data[i] = try std.fmt.parseFloat(T, bytes[j..bytes.len]);
    }

    return mat;
}

/// Subtract two matrices element-wise. The output matrix is stored in the matrix C.
pub fn sub(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (@typeInfo(T) != .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] - B.data[i];
    };

    if (@typeInfo(T) == .Struct) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].sub(B.data[i]);
    };
}

/// Transpose the matrix A. The output matrix is stored in the matrix C.
pub fn transpose(comptime T: type, C: *Matrix(T), A: Matrix(T)) void {
    for (0..A.rows) |i| for (0..A.cols) |j| {
        C.ptr(j, i).* = A.at(i, j);
    };
}
