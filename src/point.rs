use crate::Float;

#[derive(Debug, Clone, Copy)]
pub struct Point {
    pub x: Float,
    pub y: Float,
}

impl Point {
    pub const fn new(x: Float, y: Float) -> Self {
        Self { x, y }
    }
    pub fn dot(self, other: Self) -> Float {
        self.x * other.x + self.y * other.y
    }
    pub fn ortho(self) -> Self {
        Self { x: -self.y, y: self.x }
    }
}

impl core::ops::Add<Point> for Point {
    type Output = Point;
    fn add(self, other: Point) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl core::ops::Sub<Point> for Point {
    type Output = Point;
    fn sub(self, other: Point) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl core::ops::Neg for Point {
    type Output = Point;
    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}

impl core::ops::Mul<Point> for Float {
    type Output = Point;
    fn mul(self, point: Point) -> Point {
        Point {
            x: point.x * self,
            y: point.y * self,
        }
    }
}

impl core::ops::Mul<Float> for Point {
    type Output = Point;
    fn mul(self, factor: Float) -> Self {
        Self {
            x: self.x * factor,
            y: self.y * factor,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Point3 {
    pub x: Float,
    pub y: Float,
    pub z: Float,
}

impl Point3 {
    pub const fn new(x: Float, y: Float, z: Float) -> Self {
        Self { x, y, z }
    }

    pub fn dot(self, other: Self) -> Float {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    // [ x, y, z]
    // [ax,ay,az]
    // [bx,by,bz]
    pub fn cross(self, other: Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    pub fn lensq(self) -> Float {
        self.dot(self)
    }

    pub fn normalize(self) -> Self {
        self * (1.0 / self.lensq().sqrt())
    }

    pub fn d_dist_sq(self, other: Self) -> Self {
        use crate::dual::Dual;
        let dx = Dual::new(self.x - other.x, 1.0).square().dx;
        let dy = Dual::new(self.y - other.y, 1.0).square().dx;
        let dz = Dual::new(self.z - other.z, 1.0).square().dx;
        Self::new(dx, dy, dz) // gradient of distance squared
    }
}

impl core::ops::Add<Point3> for Point3 {
    type Output = Point3;
    fn add(self, other: Point3) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl core::ops::Sub<Point3> for Point3 {
    type Output = Point3;
    fn sub(self, other: Point3) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl core::ops::Neg for Point3 {
    type Output = Point3;
    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl core::ops::Mul<Point3> for Float {
    type Output = Point3;
    fn mul(self, point: Point3) -> Point3 {
        Point3 {
            x: point.x * self,
            y: point.y * self,
            z: point.z * self,
        }
    }
}

impl core::ops::Mul<Float> for Point3 {
    type Output = Point3;
    fn mul(self, factor: Float) -> Self {
        Self {
            x: self.x * factor,
            y: self.y * factor,
            z: self.z * factor,
        }
    }
}

impl core::iter::Sum<Point3> for Point3 {
    fn sum<I>(iter: I) -> Self
        where I: Iterator<Item = Point3>
    {
        let mut sum = Point3::new(0.0, 0.0, 0.0);
        for point in iter {
            sum = sum + point;
        }
        sum
    }
}