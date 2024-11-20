pub mod grid_point;
pub mod point;
pub mod icosahedron;
pub mod dual;

use std::fmt::Write;
use grid_point::GridPoint;
use point::{Point, Point3};

type P = GridPoint;
type Float = f64;

#[derive(Debug)]
pub struct TriangleSegment {
    pub base_a: GridPoint,
    pub center: GridPoint,
    pub base_x: Point,
    pub base_y: Point,
    pub scale: f64,
    pub vertex_index: Vec<Vec<usize>>,
    pub vertices: Vec<GridPoint>,
    pub pentagon_vertex: GridPoint,
    pub in_hex: Vec<GridPoint>,
    pub edge_hex: Vec<GridPoint>,
}

impl TriangleSegment {
    pub fn to_normal(&self, point: GridPoint) -> Point {
        let point_f = point.to_euclid();
        let x = self.base_x.dot(point_f);
        let y = self.base_y.dot(point_f);
        Point::new(x, y) * self.scale
    }

    pub fn lookup_vertex(
        &self,
        point: GridPoint,
        neighbors: [u8; 3],
    ) -> (usize, usize, GridPoint) {
        // eprintln!("lookup: {point:?}");
        let max_y = self.vertex_index.len();
        let max_x = self.vertex_index[0].len();
        if point.x >= 0 && point.x < max_x as i64
            && point.y >= 0 && point.y < max_y as i64
        {
            let idx = self.vertex_index[point.y as usize][point.x as usize];
            if idx != !0 {
                return (3, idx, point);
            }
        }
        let centered = point - self.center;
        // eprintln!("centered: {centered:?}");
        let (neighbor, rotated) = if self.base_a.det(point) <= 0 {
            (0, centered)
        } else if point.det(self.base_a.rot60()) <= 0 {
            (2, centered.rot120())
        } else {
            (1, centered.rotm120())
        };
        // eprintln!("rotated: {rotated:?}, neighbor: {neighbor}");
        let inside = self.base_a - 2 * self.center - rotated;
        // eprintln!("inside: {inside:?}");
        let unrotated = match neighbors[neighbor as usize] {
            0 => inside,
            1 => inside.rot120(),
            _ => inside.rotm120(),
        };
        let result = unrotated + self.center;
        // eprintln!("result: {result:?}");
        let idx = self.vertex_index[result.y as usize][result.x as usize];
        (neighbor, idx, result)
    }

    pub fn ico_vertices_2(&self) -> Vec<Point3> {
        use icosahedron::{ICO_FACES, ICO_VERTICES};
        let mut tri_vert = Vec::new();
        // let (mut mx, mut nx): (Float, Float) = (-1000.0, 1000.0);
        // let (mut my, mut ny): (Float, Float) = (-1000.0, 1000.0);
        for &vertex in &self.vertices {
            let vertex = self.to_normal(vertex);
            // mx = mx.max(vertex.x); nx = nx.min(vertex.x);
            // my = my.max(vertex.y); ny = ny.min(vertex.y);
            tri_vert.push(vertex);
        }
        // println!("({nx} {mx}) ({ny} {my}) {tri_vert:?}");
        let mut vertices = Vec::new();
        for &idx in &ICO_FACES {
            let abc = idx.map(|ii| ICO_VERTICES[ii]);
            for &vertex in &tri_vert {
                let vert = map_to_sphere_2(abc, vertex);
                vertices.push(vert);
            }
        }
        vertices
    }

    pub fn ico_vertices(&self) -> Vec<Point3> {
        const OVER_SIN_60: Float = 1.1547005383792515290182975; // 2/sqrt(3)
        use icosahedron::{ICO_FACES, ICO_VERTICES};
        let mut tri_vert = Vec::new();
        for &vertex in &self.vertices {
            tri_vert.push(self.to_normal(vertex));
        }
        let mut vertices = Vec::new();
        for &idx in &ICO_FACES {
            let [a, b, c] = idx.map(|ii| ICO_VERTICES[ii]);
            let ax_x = b - a;
            let mid_ab = 0.5 * (b + a);
            let ax_y = OVER_SIN_60 * (c - mid_ab);
            for &vertex in &tri_vert {
                let vert = vertex.x * ax_x + vertex.y * ax_y + a;
                vertices.push(vert);
            }
        }
        vertices
    }

    pub fn hex_faces(&self) -> Vec<[usize; 6]> {
        use icosahedron::{
            ICO_NEIGHBORS, ICO_EDGE_INDEX,
            find_face_orientation,
        };
        let mut hex_faces = Vec::new();
        let vert_per_face = self.vertices.len();
        for (ii, &(dirs, neigs)) in ICO_NEIGHBORS.iter().enumerate() {
            for &center in &self.in_hex {
                let mut face = [0; 6];
                let mut jj = 0;
                center.visit_neighbors(|p| {
                    let (neigh, vert, _point) = self.lookup_vertex(p, dirs);
                    let f_idx = if neigh == 3 { ii } else { neigs[neigh] };
                    face[jj] = f_idx * vert_per_face + vert;
                    jj += 1;
                });
                hex_faces.push(face);
            }
            // break;
        }
        // hex_faces.clear();
        // eprintln!("faces: {}", hex_faces.len());
        for (ii, other) in ICO_EDGE_INDEX.iter().enumerate() {
            for &jj in other {
                if jj < ii {
                    continue;
                }
                let (ori, face_idx) = find_face_orientation(ii, jj);
                for &hex in &self.edge_hex {
                    let centered = hex - self.center;
                    let vertex = match ori {
                        0 => centered,
                        1 => centered.rot120(),
                        _ => centered.rotm120(),
                    } + self.center;
                    let mut face = [0; 6];
                    let mut kk = 0;
                    let (dirs, neigs) = ICO_NEIGHBORS[face_idx];
                    vertex.visit_neighbors(|p| {
                        let (neigh, vert, _point) = self.lookup_vertex(p, dirs);
                        let f_idx = if neigh == 3 { face_idx } else { neigs[neigh] };
                        face[kk] = f_idx * vert_per_face + vert;
                        kk += 1;
                    });
                    hex_faces.push(face);
                }
            }
        }
        // eprintln!("faces: {}", hex_faces.len());
        hex_faces
    }

    pub fn penta_faces(&self) -> Vec<[usize; 5]> {
        use icosahedron::{
            ICO_FACES, ICO_EDGE_INDEX,
            find_face_orientation,
        };
        let vert_per_face = self.vertices.len();
        let centered = self.pentagon_vertex - self.center;
        let mut penta_faces = Vec::with_capacity(12);
        for (ai, other) in ICO_EDGE_INDEX.iter().enumerate() {
            let mut bi = other[0];
            let mut pentagon = [0; 5];
            for jj in 0..5 {
                let (ori, face) = find_face_orientation(bi, ai);
                let vertex = match ori {
                    0 => centered,
                    1 => centered.rot120(),
                    _ => centered.rotm120(),
                } + self.center;
                let vi = self.vertex_index[vertex.y as usize][vertex.x as usize];
                let index = (ori as usize + 2) % 3;
                let ci = ICO_FACES[face][index];
                bi = ci;
                pentagon[jj] = vi + face * vert_per_face;
            }
            penta_faces.push(pentagon);
        }
        penta_faces
    }
}

fn map_to_sphere(abc: [Point3; 3], point: Point) -> Point3 {
    const OVER_SIN_60: Float = 1.1547005383792515290182975; // 2/sqrt(3)
    let Point { mut x, mut y } = point;
    x = x * 2.0 - 1.0;
    y *= OVER_SIN_60;
    let [mut a, mut b, mut c] = abc;
    for _ii in 0..21 {
        assert!(x >= -1.0 && x <= 1.0);
        assert!(x + 1.0 >= y);
        assert!(1.0 - x >= y);
        assert!(y >= 0.0 && y <= 1.0);
        let mab = (a + b).normalize();
        let mbc = (b + c).normalize();
        let mca = (c + a).normalize();
        /*
           c
          / \
         /   \
        a_____b
        */
        if y > 0.5 {
            y -= 0.5;
            a = mca;
            b = mbc;
        } else if x > y {
            x -= 0.5;
            a = mab;
            c = mbc;
        } else if -x > y {
            x += 0.5;
            b = mab;
            c = mca;
        } else {
            x = -x;
            y = 0.5 - y;
            a = mbc;
            b = mca;
            c = mab;
        }
        x = x * 2.0;
        y = y * 2.0;
    }
    a
}
fn midpoint_circle(a: Point3, b: Point3, norm: Point3) -> Point3 {
    let dist = norm.dot(a);
    let r_circ = (1.0 - dist * dist).sqrt();
    let center = dist * norm;
    let a = a - center;
    let b = b - center;
    let mid = (a + b).normalize() * r_circ;
    let res = mid + center;
    assert!((res.dot(res) - 1.0) < 0.01, "{res:?}");
    res
}

fn map_to_sphere_2(abc: [Point3; 3], point: Point) -> Point3 {
    const OVER_SIN_60: Float = 1.1547005383792515290182975; // 2/sqrt(3)
    let Point { mut x, mut y } = point;
    x = x * 2.0 - 1.0;
    y *= OVER_SIN_60;
    let [mut a, mut b, mut c] = abc;
    a = a.normalize();
    b = b.normalize();
    c = c.normalize();
    let mut nab = a.cross(b).normalize();
    let mut nbc = b.cross(c).normalize();
    let mut nca = c.cross(a).normalize();
    for _ii in 0..21 {
        assert!((a.dot(a) - 1.0) < 0.01, "{a:?}");
        assert!((b.dot(b) - 1.0) < 0.01, "{b:?}");
        assert!((c.dot(c) - 1.0) < 0.01, "{c:?}");
        assert!(x >= -1.0 && x <= 1.0);
        assert!(x + 1.0 >= y);
        assert!(1.0 - x >= y);
        assert!(y >= 0.0 && y <= 1.0);
        let mab = midpoint_circle(a, b, nab);
        let mbc = midpoint_circle(b, c, nbc);
        let mca = midpoint_circle(c, a, nca);
        // let up = (mab - mbc).cross(mab - mca).normalize();
        let up = (a - b).cross(a - c).normalize();
        let nnab = (mbc - mca).cross(up).normalize();
        let nnbc = (mca - mab).cross(up).normalize();
        let nnca = (mab - mbc).cross(up).normalize();
        /*
           c
          / \
         /   \
        a_____b
        */
        if y > 0.5 {
            y -= 0.5;
            a = mca;
            b = mbc;
            nab = nnab;
        } else if x > y {
            x -= 0.5;
            a = mab;
            c = mbc;
            nca = nnca;
        } else if -x > y {
            x += 0.5;
            b = mab;
            c = mca;
            nbc = nnbc;
        } else {
            x = -x;
            y = 0.5 - y;
            a = mbc;
            b = mca;
            c = mab;
            nab = nnab;
            nbc = nnbc;
            nca = nnca;
        }
        x = x * 2.0;
        y = y * 2.0;
    }
    a
}
/*
q2 = np.sqrt(1/2)
q3 = np.sqrt(1/3)
s = np.matrix([[q3, 0], [0, 1]])
r = np.matrix([[q2, -q2], [q2, q2]])
ir = np.linalg.inv(r)
# transform normal grid into skewed square grid
k = np.matmul(r, np.matmul(s, ir))
# and back
ik = np.linalg.inv(k)
# M1 = (3+sqrt(3))/6
# M2 = (3-sqrt(3))/6
# np.matrix([[M1, -M2], [-M2, M1]])
# np.matrix([[ 0.78867513, -0.21132487], [-0.21132487,  0.78867513]])
c60 = 0.5
s60 = np.sin(np.pi / 3)
# rotate by 60 degrees
r60 = np.matrix([[c60, -s60], [s60, c60]])
# rotate by 60 degrees IN the square grid
'''
linear_space_1 (normal hexagonal grid)
linear_space_2 (skewed square grid with integer points)
linear_space_1 -> linear_space_2
x' = M_{1->2} * x
x = M_{2->1} * x' = inverse(M_{1->2}) * x
M_{rotate} (rotate vector by some amount) = [
    [cos(angle), -sin(angle)],
    [sin(angle), cos(angle)]
]

x_rotated = M_{rotate} * x
x'_rotated = M_{1->2} * M_{rotate} * inverse(M_{1->2}) * x'
M_{rotate_in_2} = M_{1->2} * M_{rotate} * inverse(M_{1->2})
'''
r60g = np.matmul(ik, np.matmul(r60, k))
r120 = np.matmul(r60, r60)
r120g = np.matmul(ik, np.matmul(r120, k))

r120i = np.matrix([[0, -1], [1, -1]])
r60i = np.matrix([[1, -1], [1, 0]])
rm60i = np.matrix([[0,  1], [-1, 1]])
rm120i = np.matrix([[-1,  1], [-1,  0]])
*/

// create vertices and faces on a single triangle
fn generate_tri_segment(n: i64, m: i64) -> Option<TriangleSegment> {
    // const DIR_O: P = P::new(-1, 1);
    let mut a = P::new(2 * n + m, n + 2 * m);
    const DIR_N: P = P::new(2, 1);
    const DIR_M: P = P::new(1, 2);
    assert!(a == DIR_N * n + DIR_M * m, "sanity check");
    let mut b = P::new(n - m, 2 * n + m);
    assert!(b == a.rot60(), "sanity check");
    if n < m {
        // rotate by -tau/6 in case the leftmost point has negative x
        assert!(a == (a - b).rot60(), "sanity check");
        (a, b) = (a - b, a);
    }
    let ab = b - a;
    let center = P::new((a.x + b.x) / 3, (a.y + b.y) / 3);
    let max_x = a.x;
    let max_y = b.y.max(a.y);
    // println!("a={a:?}, b={b:?}, max_x={max_x}, max_y={max_y}");
    let mut pentagon_vertex = None;
    // initialize index of vertices
    let mut vertex_index = vec![
        vec![!0; max_x as usize + 1];
        max_y as usize + 1
    ];
    let mut in_hex = Vec::new();
    let mut edge_hex = Vec::new();
    let mut vertices = Vec::new();
    let mut index = 0;

    // triangle rasterization algorithm, similar to Bresenham's line algorithm
    let mut left_x = 0;
    let mut ldet_b = 0;
    for y in 0..max_y+1 {
        // move x to the start of the edge
        while ldet_b < 0 {
            left_x += 1;
            ldet_b += b.y;
        }
        let mut x = left_x;
        let ap = P::new(x, y) - a;
        let mut det_ab = ab.x * ap.y - ab.y * ap.x;
        let mut det_a = a.x * y - a.y * x;
        let mut det_b = ldet_b;
        while det_a >= 0 && det_ab >= 0 {
            // println!("{x}, {y}, {det_a}, {det_b}, {det_ab}");
            let mod3 = (x + y) % 3;
            if mod3 == 0 {
                if det_b == 0 || det_ab == 0 {
                    //
                } else if det_a == 0 {
                    edge_hex.push(P::new(x, y));
                } else {
                    in_hex.push(P::new(x, y));
                }
            } else {
                if index == 0 {
                    pentagon_vertex = Some(P::new(x, y));
                }
                let p = P::new(x, y);

                let (line_start, line_end) = if det_a == 0 {
                    (P::new(0, 0), a)
                } else if det_b == 0 {
                    (b, P::new(0, 0))
                } else if det_ab == 0 {
                    (a, b)
                } else {
                    (p, P::new(0, 0))
                };
                let dstart = p - line_start;
                let dend = p - line_end;
                let do_store = dstart.dot(dstart) < dend.dot(dend);
                if do_store {
                    vertices.push(P::new(x, y));
                    vertex_index[y as usize][x as usize] = index;
                    index += 1;
                }
            }
            det_a -= a.y;
            det_b -= b.y;
            det_ab -= ab.y;
            x += 1;
        }
        ldet_b -= b.x;
    }

    let base_x = a.to_euclid();
    Some(TriangleSegment {
        base_a: a,
        center,
        base_x,
        base_y: base_x.ortho(),
        scale: 1.0 / base_x.dot(base_x),
        vertex_index,
        vertices,
        pentagon_vertex: pentagon_vertex?,
        in_hex,
        edge_hex,
    })
}

type LazyResult<T> = Result<T, Box<dyn std::error::Error>>;

fn print_obj<W: std::io::Write>(
    writer: &mut W, n: i64, m: i64
) -> LazyResult<()> {
    let seg = generate_tri_segment(n, m).unwrap();
    // eprintln!("{:?}", seg);
    writeln!(writer, "o Generalized Buckyball n={n} m={m}")?;
    let ico_vertices = seg.ico_vertices();
    // eprintln!("{:?}", ico_vertices);
    for &Point3 { x, y, z } in &ico_vertices {
        writeln!(writer, "v {x:.4} {y:.4} {z:.4}")?;
    }

    let hex_faces = seg.hex_faces();
    let penta_faces = seg.penta_faces();
    for face in hex_faces.iter().map(|x| &x[..])
        .chain(penta_faces.iter().map(|x| &x[..]))
    {
        let mut buf = "f".to_string();
        for idx in face.iter() {
            let _ = write!(&mut buf, " {}", idx + 1);
        }
        writeln!(writer, "{buf}")?;
    }

    Ok(())
}

fn print_tri_hierarchy<W: std::io::Write>(
    writer: &mut W, n: i64, m: i64
) -> LazyResult<()> {
    let seg = generate_tri_segment(n, m).unwrap();
    writeln!(writer, "o Generalized Buckyball n={n} m={m}")?;
    let ico_vertices = seg.ico_vertices_2();
    // eprintln!("{:?}", ico_vertices);
    for &Point3 { x, y, z } in &ico_vertices {
        writeln!(writer, "v {x:.4} {y:.4} {z:.4}")?;
    }

    let hex_faces = seg.hex_faces();
    let penta_faces = seg.penta_faces();
    for face in hex_faces.iter().map(|x| &x[..])
        .chain(penta_faces.iter().map(|x| &x[..]))
    {
        let mut buf = "f".to_string();
        for idx in face.iter() {
            let _ = write!(&mut buf, " {}", idx + 1);
        }
        writeln!(writer, "{buf}")?;
    }
    Ok(())
}

fn print_svg<W: std::io::Write>(
    writer: &mut W, n: i64, m: i64
) -> LazyResult<()> {
    let seg = generate_tri_segment(n, m).unwrap();
    // eprintln!("{:?}", seg);

    const SCALE: Float = 1000.0;
    const BASE: Point = Point::new(0.0, SCALE * 0.05);
    let a = BASE + SCALE * Point::new(1.0, 0.0);
    const SQRT_3_4: Float = 0.8660254037844386;
    let b = BASE + SCALE * Point::new(0.5, SQRT_3_4);
    let triangle = format!("M{:.2} {:.2}L{:.2} {:.2}L{:.2} {:.2}z\n",
        BASE.x, BASE.y, a.x, a.y, b.x, b.y);
    fn draw_hex(s: &mut String, seg: &TriangleSegment, hex: GridPoint) {
        let down = BASE + SCALE * seg.to_normal(hex + P::new(0, -1));
        let _ = write!(s, "M{:.2} {:.2}", down.x, down.y);
        hex.visit_neighbors(|p| {
            let n = BASE + SCALE * seg.to_normal(p);
            let _ = write!(s, "L{:.2} {:.2}", n.x, n.y).unwrap();
        });
        s.push_str("z\n");
    }
    let mut in_hex = String::new();
    for &hex in &seg.in_hex {
        draw_hex(&mut in_hex, &seg, hex);
    }
    let mut edge_hex = String::new();
    for &hex in &seg.edge_hex {
        draw_hex(&mut edge_hex, &seg, hex);
        let sa = (hex - seg.center).rot120() + seg.center;
        draw_hex(&mut edge_hex, &seg, sa);
        let sb = (sa - seg.center).rot120() + seg.center;
        draw_hex(&mut edge_hex, &seg, sb);
    }
    // let mut verts = String::new();
    let mut verts = String::new();
    for &hex in seg.in_hex.iter().chain(&seg.edge_hex) {
        let neighbors = [0, 1, 2]; // rotation of neighbors relative to us
        hex.visit_neighbors(|p| {
            let (nei, _idx, point) = seg.lookup_vertex(p, neighbors);
            // eprintln!("{p:?} {nei} {point:?}");
            let n = BASE + SCALE * seg.to_normal(point);
            let (color, r, opacity) = if nei == 3 {
                ("green", 4, 1.0)
            } else {
                ("blue", 8, 0.3)
            };
            let _ = writeln!(&mut verts,
                r#"<circle style="opacity: {}" r="{}" cx="{:.2}" cy="{:.2}" fill="{}" />"#,
                opacity, r, n.x, n.y, color
            );
        });
    }
    write!(writer, r#"<?xml version="1.0" encoding="UTF-8" standalone="no"?>
        <svg xmlns="http://www.w3.org/2000/svg"
            height="{SCALE}" width="{SCALE}">
            <path fill="none" stroke-width="2px" stroke="black" d="{triangle}" />
            <path fill="none" stroke-width="2px" stroke="red" d="{in_hex}" />
            <path fill="none" stroke-width="2px" stroke="blue" d="{edge_hex}" />
            {verts}
        </svg>
    "#)?;
    Ok(())
    // println!("{:?}", seg.to_normal(seg.base.rot60()));
}

fn optimize_vertices<W: std::io::Write>(
    writer: &mut W, n: i64, m: i64
) -> LazyResult<()> {
    use std::collections::BTreeMap;
    let seg = generate_tri_segment(n, m).unwrap();
    eprintln!("{:?}", seg);
    let mut ico_vertices = seg.ico_vertices();
    let mut edges_map = BTreeMap::<usize, Vec<usize>>::new();
    let hex_faces = seg.hex_faces();
    let penta_faces = seg.penta_faces();
    for face in hex_faces.iter().map(|x| &x[..])
        .chain(penta_faces.iter().map(|x| &x[..]))
    {
        let mut b = *face.last().unwrap();
        for &a in face {
            edges_map.entry(a).or_default().push(b);
            b = a;
        }
    }
    let mut edges = Vec::new();
    for ii in 0..ico_vertices.len() {
        let mut tmp = [0_usize; 3];
        tmp.copy_from_slice(&edges_map[&ii]);
        edges.push(tmp);
    }
    for vert in &mut ico_vertices[..] {
        *vert = vert.normalize();
    }
    let delta = 0.01;
    for _jj in 0..1000 {
        // let mut measure = 0.0;
        // for ii in 0..ico_vertices.len() {
        //     let vert = ico_vertices[ii];
        //     let diff = edges[ii].iter()
        //         .map(|&vv| (vert - ico_vertices[vv]).lensq())
        //         .sum::<Float>();
        //     measure += diff;
        // }
        // eprintln!("step={jj:4} measure={measure}");
        for ii in 0..ico_vertices.len() {
            let vert = ico_vertices[ii];
            let diff = edges[ii].iter()
                .map(|&vv| vert.d_dist_sq(ico_vertices[vv]))
                .sum::<Point3>();
            ico_vertices[ii] = (vert - diff * delta).normalize();
        }
    }

    writeln!(writer, "o Generalized Buckyball n={n} m={m}")?;
    for &Point3 { x, y, z } in &ico_vertices {
        writeln!(writer, "v {x:.4} {y:.4} {z:.4}")?;
    }
    for face in hex_faces.iter().map(|x| &x[..])
        .chain(penta_faces.iter().map(|x| &x[..]))
    {
        let mut buf = "f".to_string();
        for idx in face.iter() {
            let _ = write!(&mut buf, " {}", idx + 1);
        }
        writeln!(writer, "{buf}")?;
    }

    Ok(())
}

#[derive(PartialEq, Eq)]
enum Mode {
    TriHierarchy,
    WriteOptimized,
    WriteFlat,
    WriteTriangle
}

const MODE: Mode = Mode::TriHierarchy;
fn main() -> LazyResult<()> {
    if MODE == Mode::TriHierarchy {
        let (n, m) = (0, 20);
        let name = format!("hie_buckyball_{n:04}_{m:04}.obj");
        let mut file = std::fs::OpenOptions::new()
            .write(true).create(true).truncate(true)
            .open(&name)?;
        print_tri_hierarchy(&mut file, n, m)?;
    } else if MODE == Mode::WriteOptimized {
        let (n, m) = (0, 20);
        let name = format!("opt_general_buckyball_{n:04}_{m:04}.obj");
        let mut file = std::fs::OpenOptions::new()
            .write(true).create(true).truncate(true)
            .open(&name)?;
        optimize_vertices(&mut file, n, m)?;
    } else if MODE == Mode::WriteFlat {
        let max = 10;
        for n in 1..max {
            for m in 0..max - n {
                let name = format!("general_buckyball_{n:04}_{m:04}.obj");
                let mut file = std::fs::OpenOptions::new()
                    .write(true).create(true).truncate(true)
                    .open(&name)?;
                print_obj(&mut file, n, m)?;
            }
        }
    } else if MODE == Mode::WriteTriangle {
        for n in 1..4 {
            for m in 0..4 {
                let name = format!("general_buckyball_{n:04}_{m:04}.svg");
                let mut file = std::fs::OpenOptions::new()
                    .write(true).create(true).truncate(true)
                    .open(&name)?;
                print_svg(&mut file, n, m)?;
            }
        }
        // let (n, m) = (3, 0);
    }
    Ok(())
}
