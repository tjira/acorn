//! System module for the quantum chemistry program.

const std = @import("std");

const A2AU  = @import("constant.zig").A2AU ;
const SM2AN = @import("constant.zig").SM2AN;

const inp = @import("input.zig");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;
const uncr    = @import("helper.zig").uncr;

/// System type.
pub fn System(comptime T: type) type {
    return struct {
        coords: Matrix(T), atoms: Vector(T), charge: i32,

        /// Deinitialize the system.
        pub fn deinit(self: System(T)) void {
            self.coords.deinit(); self.atoms.deinit();
        }

        /// Clone the system. The function returns an error if the allocation of the new system fails.
        pub fn clone(self: System(T)) !System(T) {
            return System(T){.coords = try self.coords.clone(), .atoms = try self.atoms.clone(), .charge = self.charge};
        }

        /// Get the atom to basis function map.
        pub fn getBasisMap(self: System(T), basis: std.ArrayList(T), allocator: std.mem.Allocator) !std.ArrayList(usize) {
            var map = std.ArrayList(usize).init(allocator); var i: usize = 0;

            while (i < basis.items.len) {
                
                const size: usize = @intFromFloat(basis.items[i    ]);
                const   am: usize = @intFromFloat(basis.items[i + 1]);

                for (0..self.coords.rows) |j| if (self.coords.at(j, 0) == basis.items[i + 2] and self.coords.at(j, 1) == basis.items[i + 3] and self.coords.at(j, 2) == basis.items[i + 4]) {
                    for (0..2 * am + 1) |_| try map.append(j);
                };

                i += 2 * size + 5;
            }

            return map;
        }

        /// Get the coordinates of the atom at the specified index.
        pub fn getCoords(self: System(T), i: usize) [3]T {
            return .{self.coords.at(i, 0), self.coords.at(i, 1), self.coords.at(i, 2)};
        }

        /// Get the number of electrons in the system.
        pub fn getElectrons(self: System(T)) usize {
            return @as(usize, @intFromFloat(self.atoms.sum() - asfloat(T, self.charge)));
        }

        /// Calculate the nuclear repulsion energy.
        pub fn nuclearRepulsion(self: System(T)) T {
            var VNN: T = 0;

            for (0..self.atoms.rows) |i| for (0..i) |j| {

                const dx = self.coords.at(i, 0) - self.coords.at(j, 0);
                const dy = self.coords.at(i, 1) - self.coords.at(j, 1);
                const dz = self.coords.at(i, 2) - self.coords.at(j, 2);

                const r = std.math.sqrt(dx * dx + dy * dy + dz * dz);

                VNN += self.atoms.at(i) * self.atoms.at(j) / r;
            };

            return VNN;
        }

        /// Print the matrix of the system.
        pub fn print(self: System(T), device: anytype) !void {
            var buffered = std.io.bufferedWriter(device); var writer = buffered.writer();

            for (self.atoms.data, 0..) |atom, i| {

                try writer.print("{d}", .{@as(u32, @intFromFloat(atom))});

                for (0..3) |j| {
                    try writer.print(" {d:20.14}", .{self.coords.at(i, j)});
                }

                try writer.print("\n", .{});
            }

            try buffered.flush();
        }
    };
}

/// Returns the system file given the option struct and the (possibly null) system file.
pub fn load(comptime T: type, opt: inp.HartreeFockOptions(T).System, file: ?[]const u8, allocator: std.mem.Allocator) !System(T) {
    if (file == null and (opt.coords == null or opt.atoms == null)) return error.SystemNotFullySpecified;

    var system: System(T) = undefined;

    if (file == null and opt.coords != null) {

        system = System(T){
            .coords = try Matrix(T).init(opt.coords.?.len, 3, allocator),
            .atoms  = try Vector(T).init(opt.atoms .?.len,    allocator),
            .charge = opt.charge,
        };

        for (0..opt.atoms.?.len) |i| {
            system.atoms.ptr(i).* = @as(T, @floatFromInt(opt.atoms.?[i]));
        }

        for (0..opt.coords.?.len) |i| for (0..3) |j| {
            system.coords.ptr(i, j).* = opt.coords.?[i][j] * A2AU;
        };
    }

    if (file != null) system = try read(T, file.?, opt.charge, allocator);

    return system;
}

/// Read the system from the specified file.
pub fn read(comptime T: type, path: []const u8, charge: i32, allocator: std.mem.Allocator) !System(T) {
    const file = try std.fs.cwd().openFile(path, .{}); defer file.close(); var buffer: [64]u8 = undefined;

    var buffered = std.io.bufferedReader(file.reader()); var reader = buffered.reader();
    var stream   = std.io.fixedBufferStream(&buffer);  const writer =   stream.writer();

    stream.reset(); try reader.streamUntilDelimiter(writer, '\n', 64);

    const natom = try std.fmt.parseInt(u32, uncr(stream.getWritten()), 10);

    stream.reset(); try reader.streamUntilDelimiter(writer, '\n', 64);

    var coords = try Matrix(T).init(natom, 3, allocator); coords.fill(0);
    var atoms  = try Vector(T).init(natom,    allocator);  atoms.fill(0);

    for (0..natom) |i| {

        stream.reset(); try reader.streamUntilDelimiter(writer, '\n', 64);

        var it = std.mem.splitScalar(u8, uncr(stream.getWritten()), ' '); 

        while (it.next()) |token| if (token.len > 0 and atoms.at(i) == 0) {
            atoms.ptr(i).* = asfloat(T, SM2AN.get(token).?); break;
        };

        var j: i32 = 0;

        while (it.next()) |token| : (j += if (token.len > 0) 1 else 0) if (token.len > 0) {
            coords.ptr(i, @as(usize, @intCast(j))).* = try std.fmt.parseFloat(T, token) * A2AU;
        };
    }

    var nocc: u32 = 0; for (0..natom) |i| nocc += @intFromFloat(atoms.at(i));

    return System(T) {.coords = coords, .atoms  = atoms, .charge = charge};
}
