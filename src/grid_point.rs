use crate::point::Point;
use crate::Float;
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GridPoint {
    pub x: i64,
    pub y: i64,
}

impl GridPoint {
    pub const fn new(x: i64, y: i64) -> Self {
        Self { x, y }
    }

    pub fn dot(self, other: Self) -> i64 {
        self.x * other.x + self.y * other.y
    }

    pub const fn det(self, other: Self) -> i64 {
        self.x * other.y - self.y * other.x
    }

    pub fn to_euclid(self) -> Point {
        const M1: Float = 0.7886751345948128823; // (3+sqrt(3))/6
        const M2: Float = -0.2113248654051871177; // (3-sqrt(3))/6
        let x = self.x as Float;
        let y = self.y as Float;
        Point::new(M1 * x + M2 * y, M2 * x + M1 * y)
    }

    pub fn visit_neighbors<F>(self, mut f: F)
        where F: FnMut(Self)
    {
        let Self { x, y } = self;
        f(Self::new(x + 1, y));
        f(Self::new(x + 1, y + 1));
        f(Self::new(x, y + 1));
        f(Self::new(x - 1, y));
        f(Self::new(x - 1, y - 1));
        f(Self::new(x, y - 1));
    }

    pub fn rot60(self) -> Self {
        Self {
            x: self.x - self.y,
            y: self.x,
        }
    }

    pub fn rot120(self) -> Self {
        Self {
            x: -self.y,
            y: self.x - self.y,
        }
    }

    pub fn rotm120(self) -> Self {
        Self {
            x: -self.x + self.y,
            y: -self.x,
        }
    }
    /*
    pub fn rotm60(self) -> Self {
        Self {
            x: self.y,
            y: self.y - self.x,
        }
    }

    */
}

impl core::ops::Add<GridPoint> for GridPoint {
    type Output = GridPoint;
    fn add(self, other: GridPoint) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl core::ops::Sub<GridPoint> for GridPoint {
    type Output = GridPoint;
    fn sub(self, other: GridPoint) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl core::ops::Neg for GridPoint {
    type Output = GridPoint;
    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}

impl core::ops::Mul<GridPoint> for i64 {
    type Output = GridPoint;
    fn mul(self, point: GridPoint) -> GridPoint {
        GridPoint {
            x: point.x * self,
            y: point.y * self,
        }
    }
}

impl core::ops::Mul<i64> for GridPoint {
    type Output = GridPoint;
    fn mul(self, factor: i64) -> Self {
        Self {
            x: self.x * factor,
            y: self.y * factor,
        }
    }
}
