const std = @import("std");

const ClassicalDynamics = @import("classicaldynamics.zig").ClassicalDynamics;
const Matrix            = @import("matrix.zig"           ).Matrix           ;
const Vector            = @import("vector.zig"           ).Vector           ;

pub fn SurfaceHopping(comptime T: type) type {
    return struct {
        pub fn landauZener(traj_data: std.ArrayList(ClassicalDynamics(T).IterationData), time_step: T, adiabatic: bool) !?std.AutoHashMap(u8, T) {
            if (traj_data.items.len < 3) return null;
            var hop_probs = std.AutoHashMap(u8, T).init(traj_data.getLast().r.data.allocator);
            for (0..traj_data.getLast().U.rows) |i| {
                for (0..i) |j| {
                    if (traj_data.getLast().s != i and traj_data.getLast().s != j) continue;
                    const jump_state = if (traj_data.getLast().s == i) j else i;
                    const Ediff = [3]T{
                        traj_data.items[traj_data.items.len - 1].U.at(i, i) - traj_data.items[traj_data.items.len - 1].U.at(j, j),
                        traj_data.items[traj_data.items.len - 2].U.at(i, i) - traj_data.items[traj_data.items.len - 2].U.at(j, j),
                        traj_data.items[traj_data.items.len - 3].U.at(i, i) - traj_data.items[traj_data.items.len - 3].U.at(j, j)
                    };
                    const dEdiff = [2]T{(Ediff[1] - Ediff[0]) / time_step, (Ediff[2] - Ediff[1]) / time_step};
                    const ddEdiff = [1]T{(Ediff[0] - 2 * Ediff[1] + Ediff[2]) / time_step / time_step};
                    if (adiabatic) {
                        try hop_probs.put(@intCast(jump_state), std.math.exp(-0.5 * std.math.pi * std.math.sqrt(std.math.pow(T, Ediff[0], 3) / ddEdiff[0])));
                    } else {
                        try hop_probs.put(@intCast(jump_state), 1 - std.math.exp(-2 * std.math.pi * std.math.pow(T, traj_data.getLast().U.at(i, j), 2) / (if (dEdiff[0] < 0) -dEdiff[0] else dEdiff[0])));
                    }
                }
            }
            return hop_probs;
        }
    };
}
