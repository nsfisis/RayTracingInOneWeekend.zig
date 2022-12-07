const vec = @import("vec.zig");
const Vec3 = vec.Vec3;
const Point3 = vec.Point3;

pub const Ray = struct {
    origin: Vec3,
    dir: Vec3,
    time: f64,

    pub fn at(r: Ray, t: f64) Point3 {
        return r.origin.add(r.dir.mul(t));
    }
};
