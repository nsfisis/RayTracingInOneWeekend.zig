const std = @import("std");
const debug = std.debug;
const math = std.math;

const Ray = @import("ray.zig").Ray;
const vec = @import("vec.zig");
const Vec3 = vec.Vec3;
const rgb = vec.rgb;
const Color = vec.Color;
const HitRecord = @import("hit_record.zig").HitRecord;
const Texture = @import("texture.zig").Texture;
const rand = @import("rand.zig");
const Random = rand.Random;
const randomPointInUnitSphere = rand.randomPointInUnitSphere;
const randomUnitVector = rand.randomUnitVector;
const randomReal01 = rand.randomReal01;

const MaterialTag = enum {
    diffuse,
    metal,
    dielectric,
    diffuse_light,
};

pub const Material = union(MaterialTag) {
    diffuse: DiffuseMaterial,
    metal: MetalMaterial,
    dielectric: DielectricMaterial,
    diffuse_light: DiffuseLightMaterial,

    pub fn scatter(mat: Material, r_in: Ray, record: HitRecord, attenuation: *Color, scattered: *Ray, rng: Random) bool {
        return switch (mat) {
            MaterialTag.diffuse => |diffuse_mat| diffuse_mat.scatter(r_in, record, attenuation, scattered, rng),
            MaterialTag.metal => |metal_mat| metal_mat.scatter(r_in, record, attenuation, scattered, rng),
            MaterialTag.dielectric => |dielectric_mat| dielectric_mat.scatter(r_in, record, attenuation, scattered, rng),
            MaterialTag.diffuse_light => |diffuse_light_mat| diffuse_light_mat.scatter(r_in, record, attenuation, scattered, rng),
        };
    }

    pub fn emitted(mat: Material, u: f64, v: f64, p: Vec3) Color {
        return switch (mat) {
            MaterialTag.diffuse => rgb(0, 0, 0),
            MaterialTag.metal => rgb(0, 0, 0),
            MaterialTag.dielectric => rgb(0, 0, 0),
            MaterialTag.diffuse_light => |diffuse_light_mat| diffuse_light_mat.emitted(u, v, p),
        };
    }
};

pub const DiffuseMaterial = struct {
    albedo: Texture,

    fn scatter(mat: DiffuseMaterial, r_in: Ray, record: HitRecord, attenuation: *Color, scattered: *Ray, rng: Random) bool {
        var scatter_direction = record.normal.add(randomUnitVector(rng));
        if (scatter_direction.near_zero()) {
            scatter_direction = record.normal;
        }
        scattered.* = .{ .origin = record.p, .dir = scatter_direction, .time = r_in.time };
        attenuation.* = mat.albedo.value(record.u, record.v, record.p);
        return true;
    }
};

pub const MetalMaterial = struct {
    albedo: Color,
    fuzz: f64,

    fn scatter(mat: MetalMaterial, r_in: Ray, record: HitRecord, attenuation: *Color, scattered: *Ray, rng: Random) bool {
        debug.assert(mat.fuzz <= 1.0);
        const reflected = reflect(r_in.dir.normalized(), record.normal);
        scattered.* = .{ .origin = record.p, .dir = reflected.add(randomPointInUnitSphere(rng).mul(mat.fuzz)), .time = r_in.time };
        attenuation.* = mat.albedo;
        return reflected.dot(record.normal) > 0.0;
    }
};

pub const DielectricMaterial = struct {
    // index of refraction.
    ir: f64,

    fn scatter(mat: DielectricMaterial, r_in: Ray, record: HitRecord, attenuation: *Color, scattered: *Ray, rng: Random) bool {
        const refraction_ratio = if (record.front_face) 1.0 / mat.ir else mat.ir;
        const unit_dir = r_in.dir.normalized();

        const cos_theta = @min(unit_dir.mul(-1.0).dot(record.normal), 1.0);
        const sin_theta = @sqrt(1.0 - cos_theta * cos_theta);
        // sin(theta') = refraction_ratio * sin(theta) <= 1
        const can_refract = refraction_ratio * sin_theta <= 1.0;

        const dir = if (can_refract and reflectance(cos_theta, refraction_ratio) < randomReal01(rng)) refract(unit_dir, record.normal, refraction_ratio) else reflect(unit_dir, record.normal);
        scattered.* = .{ .origin = record.p, .dir = dir, .time = r_in.time };
        attenuation.* = rgb(1.0, 1.0, 1.0);
        return true;
    }

    fn reflectance(cos: f64, refraction_idx: f64) f64 {
        const r0 = (1.0 - refraction_idx) / (1.0 + refraction_idx);
        const r1 = r0 * r0;
        return r1 + (1.0 - r1) * math.pow(f64, 1.0 - cos, 5.0);
    }
};

pub const DiffuseLightMaterial = struct {
    emit: Texture,

    fn scatter(mat: DiffuseLightMaterial, r_in: Ray, record: HitRecord, attenuation: *Color, scattered: *Ray, rng: Random) bool {
        _ = mat;
        _ = r_in;
        _ = record;
        _ = attenuation;
        _ = scattered;
        _ = rng;
        return false;
    }

    fn emitted(mat: DiffuseLightMaterial, u: f64, v: f64, p: Vec3) Color {
        return mat.emit.value(u, v, p);
    }
};

fn reflect(v: Vec3, n: Vec3) Vec3 {
    return v.sub(n.mul(2 * v.dot(n)));
}

fn refract(uv: Vec3, n: Vec3, etai_over_etat: f64) Vec3 {
    const cos_theta = @min(uv.mul(-1.0).dot(n), 1.0);
    const r_out_perpendicular = uv.add(n.mul(cos_theta)).mul(etai_over_etat);
    const r_out_parallel = n.mul(-@sqrt(@abs(1.0 - r_out_perpendicular.normSquared())));
    return r_out_perpendicular.add(r_out_parallel);
}
