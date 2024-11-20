use crate::Float;

#[derive(Debug, Clone, Copy)]
pub struct Dual {
    pub x: Float,
    pub dx: Float,
}
impl Dual {
    pub fn new(x: Float, dx: Float) -> Self {
        Self { x, dx }
    }
    pub fn square(self) -> Self {
        self * self
    }
}

impl core::ops::Sub<Dual> for Dual {
    type Output = Dual;
    fn sub(self, other: Dual) -> Self {
        Self {
            x: self.x - other.x,
            dx: self.dx - other.dx,
        }
    }
}

impl core::ops::Neg for Dual {
    type Output = Dual;
    fn neg(self) -> Self {
        Self {
            x: -self.x,
            dx: -self.dx,
        }
    }
}

impl core::ops::Add<Dual> for Dual {
    type Output = Dual;
    fn add(self, other: Dual) -> Dual {
        Dual {
            x: self.x + other.x,
            dx: self.dx + other.dx,
        }
    }
}

impl core::ops::Mul<Dual> for Dual {
    type Output = Dual;
    fn mul(self, other: Dual) -> Dual {
        Dual {
            x: self.x * other.x,
            dx: self.x * other.dx + self.dx * other.x,
        }
    }
}
