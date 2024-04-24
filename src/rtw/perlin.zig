const std = @import("std");
const ArrayList = std.ArrayList;

const rand = @import("rand.zig");
const Random = rand.Random;
const randomReal01 = rand.randomReal01;
const randomInt = rand.randomInt;
const Point3 = @import("vec.zig").Point3;
const Vec3 = @import("vec.zig").Vec3;

pub const Perlin = struct {
    const POINT_COUNT = 256;

    randomVec: ArrayList(Vec3),
    permX: ArrayList(usize),
    permY: ArrayList(usize),
    permZ: ArrayList(usize),

    pub fn init(allocator: std.mem.Allocator, rng: Random) !Perlin {
        var perlin = Perlin{
            .randomVec = try ArrayList(Vec3).initCapacity(allocator, POINT_COUNT),
            .permX = try ArrayList(usize).initCapacity(allocator, POINT_COUNT),
            .permY = try ArrayList(usize).initCapacity(allocator, POINT_COUNT),
            .permZ = try ArrayList(usize).initCapacity(allocator, POINT_COUNT),
        };

        var i: usize = 0;
        while (i < POINT_COUNT) : (i += 1) {
            try perlin.randomVec.append(Vec3.random(rng, -1, 1).normalized());
            try perlin.permX.append(i);
            try perlin.permY.append(i);
            try perlin.permZ.append(i);
        }
        permute(rng, perlin.permX.items, POINT_COUNT);
        permute(rng, perlin.permY.items, POINT_COUNT);
        permute(rng, perlin.permZ.items, POINT_COUNT);

        return perlin;
    }

    pub fn deinit(perlin: Perlin) void {
        perlin.permZ.deinit();
        perlin.permY.deinit();
        perlin.permX.deinit();
        perlin.randomVec.deinit();
    }

    pub fn noise(perlin: Perlin, p: Point3) f64 {
        const u = p.x - @floor(p.x);
        const v = p.y - @floor(p.y);
        const w = p.z - @floor(p.z);
        const u_ = u * u * (3 - 2 * u);
        const v_ = v * v * (3 - 2 * v);
        const w_ = w * w * (3 - 2 * w);

        const i = @as(i32, @intFromFloat(@floor(p.x)));
        const j = @as(i32, @intFromFloat(@floor(p.y)));
        const k = @as(i32, @intFromFloat(@floor(p.z)));

        var c: [2][2][2]Vec3 = undefined;

        var di: usize = 0;
        while (di < 2) : (di += 1) {
            var dj: usize = 0;
            while (dj < 2) : (dj += 1) {
                var dk: usize = 0;
                while (dk < 2) : (dk += 1) {
                    const ix = @as(usize, @intCast(i + @as(i32, @intCast(di)))) & 255;
                    const iy = @as(usize, @intCast(j + @as(i32, @intCast(dj)))) & 255;
                    const iz = @as(usize, @intCast(k + @as(i32, @intCast(dk)))) & 255;
                    c[di][dj][dk] = perlin.randomVec.items[
                        perlin.permX.items[ix] ^ perlin.permY.items[iy] ^ perlin.permZ.items[iz]
                    ];
                }
            }
        }

        return perlinInterp(c, u_, v_, w_);
    }

    pub fn turb(perlin: Perlin, p: Point3, depth: u8) f64 {
        var accum: f64 = 0.0;
        var p_ = p;
        var weight: f64 = 1.0;
        var i: u8 = 0;
        while (i < depth) : (i += 1) {
            accum += weight * perlin.noise(p_);
            weight *= 0.5;
            p_ = p_.mul(2.0);
        }
        return @abs(accum);
    }

    fn permute(rng: Random, p: []usize, n: usize) void {
        var i = n - 1;
        while (i > 0) : (i -= 1) {
            const target = randomInt(usize, rng, 0, i);
            const tmp = p[i];
            p[i] = p[target];
            p[target] = tmp;
        }
    }

    fn perlinInterp(c: [2][2][2]Vec3, u: f64, v: f64, w: f64) f64 {
        var accum: f64 = 0.0;
        var i: usize = 0;
        while (i < 2) : (i += 1) {
            var j: usize = 0;
            while (j < 2) : (j += 1) {
                var k: usize = 0;
                while (k < 2) : (k += 1) {
                    const ti = @as(f64, @floatFromInt(i));
                    const tj = @as(f64, @floatFromInt(j));
                    const tk = @as(f64, @floatFromInt(k));
                    const weight = Vec3{ .x = u - ti, .y = v - tj, .z = w - tk };
                    accum +=
                        (ti * u + (1.0 - ti) * (1.0 - u)) *
                        (tj * v + (1.0 - tj) * (1.0 - v)) *
                        (tk * w + (1.0 - tk) * (1.0 - w)) *
                        Vec3.dot(c[i][j][k], weight);
                }
            }
        }
        return accum;
    }
};
