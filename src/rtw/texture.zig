const std = @import("std");

const Color = @import("vec.zig").Color;
const Vec3 = @import("vec.zig").Vec3;

const TextureTag = enum {
    solid,
    checker,
};

pub const Texture = union(TextureTag) {
    solid: SolidTexture,
    checker: CheckerTexture,

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

    pub fn value(tx: Texture, u: f64, v: f64, p: Vec3) Color {
        return switch (tx) {
            TextureTag.solid => |solidTx| solidTx.value(u, v, p),
            TextureTag.checker => |checkerTx| checkerTx.value(u, v, p),
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
