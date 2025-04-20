const std = @import("std");

const Vec3 = @import("vec.zig").Vec3;

pub const Random = std.Random;

// [min, max)
pub fn randomInt(comptime T: type, rand: Random, min: T, max: T) T {
    return rand.intRangeLessThan(T, min, max);
}

// [0, 1)
pub fn randomReal01(rand: Random) f64 {
    return rand.float(f64);
}

// [min, max)
pub fn randomReal(rand: Random, min: f64, max: f64) f64 {
    return min + randomReal01(rand) * (max - min);
}

pub fn randomPointInUnitSphere(rand: Random) Vec3 {
    while (true) {
        const p = Vec3.random(rand, -1.0, 1.0);
        if (p.norm() >= 1) continue;
        return p;
    }
}

pub fn randomPointInUnitDisk(rand: Random) Vec3 {
    while (true) {
        const p = Vec3{ .x = randomReal(rand, -1.0, 1.0), .y = randomReal(rand, -1.0, 1.0), .z = 0.0 };
        if (p.norm() >= 1) continue;
        return p;
    }
}

pub fn randomUnitVector(rand: Random) Vec3 {
    return randomPointInUnitSphere(rand).normalized();
}
