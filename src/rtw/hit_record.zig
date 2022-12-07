const vec = @import("vec.zig");
const material = @import("material.zig");
const Vec3 = vec.Vec3;
const Point3 = vec.Point3;
const Material = material.Material;

pub const HitRecord = struct {
    // The point where the ray and the hittable hits.
    p: Point3,
    // The normal of the hittable at p.
    normal: Vec3,
    // The material at p.
    material: *const Material,
    // p = ray.at(t)
    t: f64,
    // The coordinate of the surface where the ray intersects.
    u: f64,
    v: f64,
    // True if the ray hits the hittable from the front face, i.e., outside of it.
    front_face: bool,
};
