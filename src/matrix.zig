//! This module provides a matrix class and functions to manipulate matrices.

const std = @import("std"); const cwp = @import("cwrapper.zig");

const hlp = @import("helper.zig");
const inp = @import("input.zig" );
const mth = @import("math.zig"  );
const vec = @import("vector.zig");

const StridedArray = @import("strided_array.zig").StridedArray;
const Tensor       = @import("tensor.zig"       ).Tensor      ;
const Vector       = @import("vector.zig"       ).Vector      ;

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
        pub fn complex(self: Matrix(T)) !Matrix(std.math.Complex(T)) {
            var other = try Matrix(std.math.Complex(T)).init(self.rows, self.cols, self.allocator);

            for (0..self.rows) |i| for (0..self.cols) |j| {
                other.ptr(i, j).* = std.math.Complex(T).init(self.at(i, j), 0);
            };

            return other;
        }

        /// Calculate the product of the diagonal elements of the matrix.
        pub fn diagProd(self: Matrix(T)) T {
            var result: T = 1;

            for (0..@min(self.rows, self.cols)) |i| {
                result *= self.at(i, i);
            }

            return result;
        }

        /// Check if two matrices are equal. The function returns true if the matrices are equal, false otherwise.
        pub fn eq(self: Matrix(T), other: Matrix(T)) bool {
            if (self.rows != other.rows or self.cols != other.cols) return false;

            for (0..self.data.len) |i| {
                if (self.data[i] != other.data[i]) return false;
            }

            return true;
        }

        /// Fill the matrix with a given value.
        pub fn fill(self: Matrix(T), value: T) void {
            @memset(self.data, value);
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
                self.data[i] = start + hlp.asfloat(T, i) * (end - start) / hlp.asfloat(T, self.rows * self.cols - 1);
            }
        }

        /// Copy the data to another matrix.
        pub fn memcpy(self: Matrix(T), dest: Matrix(T)) void {
            @memcpy(dest.data, self.data);
        }

        /// Multiply the matrix with a scalar.
        pub fn muls(self: Matrix(T), s: T) void {
            if (comptime !hlp.istruct(T)) for (0..self.data.len) |i| {
                self.data[i] *= s;
            };

            if (comptime hlp.istruct(T)) for (0..self.data.len) |i| {
                self.data[i] = self.data[i].mul(s);
            };
        }

        /// Calculates the 1-norm of the matrix.
        pub fn onorm(self: Matrix(T)) T {
            var result: T = 0;

            for (0..self.cols) |j| {

                const colsum = self.column(j).absSum();

                if (colsum > result) {
                    result = colsum;
                }
            }

            return result;
        }

        /// Print the matrix to the given device.
        pub fn print(self: Matrix(T), device: anytype) !void {
            var buffer: [32768]u8 = undefined; var writer = device.writer(&buffer); var writer_interface = &writer.interface;

            try writer_interface.print("{d} {d}\n", .{self.rows, self.cols});

            for (self.data, 1..) |e, i| {

                if (comptime hlp.istruct(T)) {
                    try writer_interface.print("({d:20.14}, {d:20.14}){s}", .{e.re, e.im, if(i % self.cols == 0) "\n" else " "});
                }

                else try writer_interface.print("{d:20.14}{s}", .{e, if(i % self.cols == 0) "\n" else " "});
            }

            try writer_interface.flush();
        }

        /// Access the element at the given row and column as a mutable reference.
        pub fn ptr(self: Matrix(T), i: usize, j: usize) *T {
            return &self.data[i * self.cols + j];
        }

        /// Fill the matrix with normally distributed random numbers.
        pub fn randn(self: Matrix(T), mean: T, stdev: T, seed: usize) void {
            var rng = std.Random.DefaultPrng.init(seed); const random = rng.random();

            for (0..self.data.len) |i| {
                self.data[i] = mean + stdev * random.floatNorm(T);
            }
        }

        /// Fill the matrix with uniformly distributed random numbers.
        pub fn randu(self: Matrix(T), low: T, high: T, seed: usize) void {
            var rng = std.Random.DefaultPrng.init(seed); const random = rng.random();

            for (0..self.data.len) |i| {
                self.data[i] = low + (high - low) * random.float(T);
            }
        }

        /// Reshape the matrix to a new number of rows and columns. The function returns an error if the new shape is not compatible with the current data.
        pub fn reshape(self: *Matrix(T), nm: usize, nn: usize) !void {
            if (nm * nn != self.rows * self.cols) return error.IncompatibleShape;

            self.rows = nm; self.cols = nn;
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

        /// Symmetrize the matrix as A = (A + A^T) / 2.
        pub fn symmetrize(self: Matrix(T)) !void {
            if (self.rows != self.cols) return error.NotSquareMatrixCannotSymmetrize;

            for (0..self.rows) |i| {
                for (i + 1..self.cols) |j| {

                    const val = (self.at(i, j) + self.at(j, i)) / 2;

                    self.ptr(i, j).* = val;
                    self.ptr(j, i).* = val;
                }
            }
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

            try self.print(file);
        }
    };
}

/// The main function for matrix manipulation.
pub fn run(comptime T: type, opt: inp.MatrixOptions(T), print: bool, allocator: std.mem.Allocator) !void {
    var A: ?Matrix(T) = null; var B: ?Matrix(T) = null; var C: ?Matrix(T) = null; var outputs = std.ArrayList(Matrix(T)){}; defer outputs.deinit(allocator);

    var timer = try std.time.Timer.start();

    if (opt.random != null) {
        for (0..opt.random.?.count) |i| {

            A = try Matrix(T).init(opt.random.?.dims[0], opt.random.?.dims[1], allocator); defer A.?.deinit();

            if (std.mem.eql(u8, opt.random.?.distribution, "normal")) {
                A.?.randn(opt.random.?.parameters[0], opt.random.?.parameters[1], opt.random.?.seed + i);
            }

            else return error.UnknownDistribution;

            if (opt.outputs != null and opt.outputs.?.len > i) {
                try outputs.append(allocator, try A.?.clone());
            }

            if (print and opt.print) {
                try hlp.print(std.fs.File.stdout(), "\nGENERATED MATRIX #{d:2}\n", .{i}); try A.?.print(std.fs.File.stdout());
            }
        }

        if (print) {
            try hlp.print(std.fs.File.stdout(), "\nGENERATING RANDOM MATRICES: {D}\n", .{timer.read()});
        }
    }

    if (opt.multiply != null and opt.inputs != null and opt.inputs.?.len > 0) {

        A = try read(T, opt.inputs.?[opt.inputs.?.len - 1], allocator); defer A.?.deinit();

        for (0..opt.inputs.?.len - 1) |i| {

            B = try read(T, opt.inputs.?[opt.inputs.?.len - i - 2], allocator); defer B.?.deinit();

            C = try Matrix(T).init(B.?.rows, A.?.cols, allocator); defer C.?.deinit();

            try cwp.Blas(T).dgemm(&C.?, B.?, false, A.?, false); C.?.memcpy(A.?);
        }

        if (print and opt.print) {
            try hlp.print(std.fs.File.stdout(), "\nRESULT\n", .{}); try A.?.print(std.fs.File.stdout());
        }

        if (opt.outputs != null) {
            try outputs.append(allocator, try A.?.clone());
        }

        if (print) {
            try hlp.print(std.fs.File.stdout(), "\nCALCULATING MATRIX PRODUCT: {D}\n", .{timer.read()});
        }
    }

    timer = try std.time.Timer.start();

    for (outputs.items, 0..) |e, i| {
        try e.write(opt.outputs.?[i]);
    }

    if (print and opt.outputs != null) {
        try hlp.print(std.fs.File.stdout(), "SAVING THE OUTPUT MATRICES: {D}\n", .{timer.read()});
    }

    if (opt.outputs != null) for (outputs.items) |*e| {
        e.*.deinit();
    };
}

/// Add two matrices element-wise. The output matrix is stored in the matrix C.
pub fn add(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (comptime !hlp.istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] + B.data[i];
    };

    if (comptime hlp.istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].add(B.data[i]);
    };
}

/// Add a scalar to a matrix element-wise. The output matrix is stored in the matrix C.
pub fn adds(comptime T: type, C: *Matrix(T), A: Matrix(T), s: T) void {
    if (comptime !hlp.istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] + s;
    };

    if (comptime hlp.istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].add(s);
    };
}

/// Divide a matrix by a scalar. The output matrix is stored in the matrix C.
pub fn divs(comptime T: type, C: *Matrix(T), A: Matrix(T), s: T) void {
    if (comptime !hlp.istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] / s;
    };

    if (comptime hlp.istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].div(s);
    };
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

/// Multiply two matrices element-wise. The output matrix is stored in the matrix C.
pub fn mul(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (comptime !hlp.istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] * B.data[i];
    };

    if (comptime hlp.istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].mul(B.data[i]);
    };
}

/// Multiply a matrix by a scalar. The output matrix is stored in the matrix C.
pub fn muls(comptime T: type, C: *Matrix(T), A: Matrix(T), s: T) void {
    if (comptime !hlp.istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] * s;
    };

    if (comptime hlp.istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].mul(s);
    };
}

/// Read a matrix from a file. The function returns an error if the file cannot be opened or if the matrix cannot be read.
pub fn read(comptime T: type, path: []const u8, allocator: std.mem.Allocator) !Matrix(T) {
    const file = try std.fs.cwd().openFile(path, .{}); defer file.close();

    var buffer: [32768]u8 = undefined; var reader = file.reader(&buffer); var reader_interface = &reader.interface;

    const rowstr = try reader_interface.takeDelimiterExclusive(' ' ); const rows = try std.fmt.parseInt(usize, rowstr, 10);
    const colstr = try reader_interface.takeDelimiterExclusive('\n'); const cols = try std.fmt.parseInt(usize, colstr, 10);

    const mat = try Matrix(T).init(rows, cols, allocator); var i: usize = 0;

    while (true) {

        const line = reader_interface.takeDelimiterExclusive('\n') catch {break;};

        var it = std.mem.tokenizeAny(u8, line, " "); 

        while (it.next()) |element| {
            mat.data[i] = try std.fmt.parseFloat(T, element); i += 1;
        }
    }

    return mat;
}

/// Subtract two matrices element-wise. The output matrix is stored in the matrix C.
pub fn sub(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    if (comptime !hlp.istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i] - B.data[i];
    };

    if (comptime hlp.istruct(T)) for (0..C.data.len) |i| {
        C.data[i] = A.data[i].sub(B.data[i]);
    };
}
