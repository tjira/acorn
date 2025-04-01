//! System module for the quantum chemistry program.

const std = @import("std");

const A2AU  = @import("constant.zig").A2AU ;
const SM2AN = @import("constant.zig").SM2AN;

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

const asfloat = @import("helper.zig").asfloat;
const uncr    = @import("helper.zig").uncr;

/// System type.
pub fn System(comptime T: type) type {
    return struct {
        coords: Matrix(T), atoms: Vector(T), nocc: u32,

        /// Deinitialize the system.
        pub fn deinit(self: System(T)) void {
            self.coords.deinit(); self.atoms.deinit();
        }

        /// Get the coordinates of the atom at the specified index.
        pub fn getCoords(self: System(T), i: usize) [3]T {
            return .{self.coords.at(i, 0), self.coords.at(i, 1), self.coords.at(i, 2)};
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

        /// Read the system from the specified file.
        pub fn read(path: []const u8, allocator: std.mem.Allocator) !System(T) {
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

            return System(T) {.coords = coords, .atoms  = atoms, .nocc = nocc / 2};
        }
    };
}
