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
        const u = p.x - @floor(p.x);
        const v = p.y - @floor(p.y);
        const w = p.z - @floor(p.z);

        const i = @floatToInt(i32, @floor(p.x));
        const j = @floatToInt(i32, @floor(p.y));
        const k = @floatToInt(i32, @floor(p.z));

        var c: [2][2][2]f64 = undefined;

        var di: usize = 0;
        while (di < 2) : (di += 1) {
            var dj: usize = 0;
            while (dj < 2) : (dj += 1) {
                var dk: usize = 0;
                while (dk < 2) : (dk += 1) {
                    const ix = @intCast(usize, i + @intCast(i32, di)) & 255;
                    const iy = @intCast(usize, j + @intCast(i32, dj)) & 255;
                    const iz = @intCast(usize, k + @intCast(i32, dk)) & 255;
                    c[di][dj][dk] = perlin.randomNoise[perlin.permX[ix] ^ perlin.permY[iy] ^ perlin.permZ[iz]];
                }
            }
        }

        return trilinearInterp(c, u, v, w);
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

    fn trilinearInterp(c: [2][2][2]f64, u: f64, v: f64, w: f64) f64 {
        var accum: f64 = 0.0;
        var i: usize = 0;
        while (i < 2) : (i += 1) {
            var j: usize = 0;
            while (j < 2) : (j += 1) {
                var k: usize = 0;
                while (k < 2) : (k += 1) {
                    const ti = @intToFloat(f64, i);
                    const tj = @intToFloat(f64, j);
                    const tk = @intToFloat(f64, k);
                    accum +=
                        (ti * u + (1.0 - ti) * (1.0 - u)) *
                        (tj * v + (1.0 - tj) * (1.0 - v)) *
                        (tk * w + (1.0 - tk) * (1.0 - w)) *
                        c[i][j][k];
                }
            }
        }
        return accum;
    }
};
