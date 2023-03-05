const std = @import("std");
const Image = @import("zigimg").Image;

const Color = @import("vec.zig").Color;
const Vec3 = @import("vec.zig").Vec3;
const rgb = @import("vec.zig").rgb;
const Random = @import("rand.zig").Random;
const Perlin = @import("perlin.zig").Perlin;

const TextureTag = enum {
    solid,
    checker,
    noise,
    image,
};

pub const Texture = union(TextureTag) {
    solid: SolidTexture,
    checker: CheckerTexture,
    noise: NoiseTexture,
    image: ImageTexture,

    pub fn makeSolid(color: Color) Texture {
        return .{ .solid = .{ .color = color } };
    }

    pub fn makeChecker(allocator: std.mem.Allocator, odd: Color, even: Color) !Texture {
        return .{ .checker = try CheckerTexture.init(
            allocator,
            Texture.makeSolid(odd),
            Texture.makeSolid(even),
        ) };
    }

    pub fn makeNoise(allocator: std.mem.Allocator, scale: f64, rng: Random) !Texture {
        return .{ .noise = try NoiseTexture.init(allocator, rng, scale) };
    }

    pub fn makeImage(allocator: std.mem.Allocator, file_path: []const u8) !Texture {
        return .{ .image = try ImageTexture.init(allocator, file_path) };
    }

    pub fn value(tx: Texture, u: f64, v: f64, p: Vec3) Color {
        return switch (tx) {
            TextureTag.solid => |solidTx| solidTx.value(u, v, p),
            TextureTag.checker => |checkerTx| checkerTx.value(u, v, p),
            TextureTag.noise => |noiseTx| noiseTx.value(u, v, p),
            TextureTag.image => |imageTx| imageTx.value(u, v, p),
        };
    }
};

pub const SolidTexture = struct {
    color: Color,

    fn value(tx: SolidTexture, u: f64, v: f64, p: Vec3) Color {
        _ = u;
        _ = v;
        _ = p;
        return tx.color;
    }
};

pub const CheckerTexture = struct {
    allocator: std.mem.Allocator,
    odd: *const Texture,
    even: *const Texture,

    fn init(allocator: std.mem.Allocator, odd: Texture, even: Texture) !CheckerTexture {
        var odd_ = try allocator.create(Texture);
        var even_ = try allocator.create(Texture);
        odd_.* = odd;
        even_.* = even;
        return .{
            .allocator = allocator,
            .odd = odd_,
            .even = even_,
        };
    }

    fn deinit(tx: CheckerTexture) void {
        tx.allocator.destroy(tx.even);
        tx.allocator.destroy(tx.odd);
    }

    fn value(tx: CheckerTexture, u: f64, v: f64, p: Vec3) Color {
        const sines = @sin(10 * p.x) * @sin(10 * p.y) * @sin(10 * p.z);
        return if (sines < 0) tx.odd.value(u, v, p) else tx.even.value(u, v, p);
    }
};

pub const NoiseTexture = struct {
    perlin: Perlin,
    scale: f64,

    fn init(allocator: std.mem.Allocator, rng: Random, scale: f64) !NoiseTexture {
        return .{
            .perlin = try Perlin.init(allocator, rng),
            .scale = scale,
        };
    }

    fn deinit(tx: NoiseTexture) void {
        tx.perlin.deinit();
    }

    fn value(tx: NoiseTexture, u: f64, v: f64, p: Vec3) Color {
        _ = u;
        _ = v;
        return rgb(1, 1, 1).mul(0.5 * (1.0 + @sin(tx.scale * p.z + 10.0 * tx.perlin.turb(p, 7))));
    }
};

pub const ImageTexture = struct {
    image: Image,

    fn init(allocator: std.mem.Allocator, file_path: []const u8) !ImageTexture {
        const image = try Image.fromFilePath(allocator, file_path);
        return .{
            .image = image,
        };
    }

    fn deinit(tx: ImageTexture) void {
        tx.image.deinit();
    }

    fn value(tx: ImageTexture, u: f64, v: f64, p: Vec3) Color {
        _ = p;
        // Clamp input texture coordinates to [0, 1] x [1, 0]
        const u_ = std.math.clamp(u, 0.0, 1.0);
        const v_ = 1.0 - std.math.clamp(v, 0.0, 1.0); // Flip v to image coordinates
        const i = @floatToInt(usize, u_ * @intToFloat(f64, tx.image.width));
        const j = @floatToInt(usize, v_ * @intToFloat(f64, tx.image.height));
        // Clamp integer mapping, since actual coordinates should be less than 1.0
        const i_ = @min(i, tx.image.width - 1);
        const j_ = @min(j, tx.image.width - 1);
        const color_scale = 1.0 / 255.0;
        const pixels = tx.image.pixels.asBytes();
        const offset = j_ * tx.image.width * 4 + i_ * 4;
        const r = @intToFloat(f64, pixels[offset + 0]);
        const g = @intToFloat(f64, pixels[offset + 1]);
        const b = @intToFloat(f64, pixels[offset + 2]);
        const a = @intToFloat(f64, pixels[offset + 3]);
        if (a == 0) {
            // Ocean
            return rgb(0, 0, 1.0);
        } else {
            return rgb(color_scale * r, color_scale * g, color_scale * b);
        }
    }
};
