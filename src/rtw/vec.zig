const std = @import("std");

const rand = @import("rand.zig");
const Random = rand.Random;
const randomReal = rand.randomReal;
const randomReal01 = rand.randomReal01;

pub const Vec3 = struct {
    x: f64,
    y: f64,
    z: f64,

    pub fn norm(v: Vec3) f64 {
        return @sqrt(v.normSquared());
    }

    pub fn normSquared(v: Vec3) f64 {
        return v.x * v.x + v.y * v.y + v.z * v.z;
    }

    pub fn dot(u: Vec3, v: Vec3) f64 {
        return u.x * v.x + u.y * v.y + u.z * v.z;
    }

    pub fn cross(u: Vec3, v: Vec3) Vec3 {
        return Vec3{
            .x = u.y * v.z - u.z * v.y,
            .y = u.z * v.x - u.x * v.z,
            .z = u.x * v.y - u.y * v.x,
        };
    }

    pub fn normalized(v: Vec3) Vec3 {
        const n = v.norm();
        if (n == 0.0) {
            return v;
        } else {
            return v.div(n);
        }
    }

    pub fn add(u: Vec3, v: Vec3) Vec3 {
        return Vec3{
            .x = u.x + v.x,
            .y = u.y + v.y,
            .z = u.z + v.z,
        };
    }

    pub fn sub(u: Vec3, v: Vec3) Vec3 {
        return Vec3{
            .x = u.x - v.x,
            .y = u.y - v.y,
            .z = u.z - v.z,
        };
    }

    pub fn mul(v: Vec3, t: f64) Vec3 {
        return Vec3{
            .x = v.x * t,
            .y = v.y * t,
            .z = v.z * t,
        };
    }

    pub fn mulV(u: Vec3, v: Vec3) Vec3 {
        return Vec3{
            .x = u.x * v.x,
            .y = u.y * v.y,
            .z = u.z * v.z,
        };
    }

    pub fn div(v: Vec3, t: f64) Vec3 {
        return Vec3{
            .x = v.x / t,
            .y = v.y / t,
            .z = v.z / t,
        };
    }

    pub fn random01(rng: Random) Vec3 {
        return Vec3{
            .x = randomReal01(rng),
            .y = randomReal01(rng),
            .z = randomReal01(rng),
        };
    }

    pub fn random(rng: Random, min: f64, max: f64) Vec3 {
        return Vec3{
            .x = randomReal(rng, min, max),
            .y = randomReal(rng, min, max),
            .z = randomReal(rng, min, max),
        };
    }

    pub fn near_zero(v: Vec3) bool {
        const epsilon = 1e-8;
        return @fabs(v.x) < epsilon and @fabs(v.y) < epsilon and @fabs(v.z) < epsilon;
    }
};

pub const Point3 = Vec3;
pub const Color = Vec3;

pub fn rgb(r: f64, g: f64, b: f64) Color {
    return .{ .x = r, .y = g, .z = b };
}
