const rand = @import("rand.zig");
const Random = rand.Random;
const randomReal01 = rand.randomReal01;
const randomInt = rand.randomInt;
const Point3 = @import("vec.zig").Point3;

pub const Perlin = struct {
    const POINT_COUNT = 256;

    randomNoise: [POINT_COUNT]f64,
    permX: [POINT_COUNT]usize,
    permY: [POINT_COUNT]usize,
    permZ: [POINT_COUNT]usize,

    pub fn init(rng: Random) Perlin {
        var perlin = Perlin{
            .randomNoise = undefined,
            .permX = undefined,
            .permY = undefined,
            .permZ = undefined,
        };

        var i: usize = 0;
        while (i < POINT_COUNT) : (i += 1) {
            perlin.randomNoise[i] = randomReal01(rng);
        }
        perlinGeneratePerm(rng, &perlin.permX);
        perlinGeneratePerm(rng, &perlin.permY);
        perlinGeneratePerm(rng, &perlin.permZ);

        return perlin;
    }

    pub fn noise(perlin: Perlin, p: Point3) f64 {
        const i = @intCast(usize, @floatToInt(i32, 4.0 * p.x) & 255);
        const j = @intCast(usize, @floatToInt(i32, 4.0 * p.y) & 255);
        const k = @intCast(usize, @floatToInt(i32, 4.0 * p.z) & 255);

        return perlin.randomNoise[perlin.permX[i] ^ perlin.permY[j] ^ perlin.permZ[k]];
    }

    fn perlinGeneratePerm(rng: Random, p: []usize) void {
        var i: usize = 0;
        while (i < POINT_COUNT) : (i += 1) {
            p[i] = i;
        }
        permute(rng, p, POINT_COUNT);
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
};
