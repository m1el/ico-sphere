use crate::Float;
use crate::point::Point3;
const PHI: Float = 1.6180339887498948482;
pub const ICO_VERTICES: [Point3; 12] = [
    Point3::new(0.0, 1.0, PHI),
    Point3::new(0.0, 1.0, -PHI),
    Point3::new(0.0, -1.0, PHI),
    Point3::new(0.0, -1.0, -PHI),

    Point3::new(PHI, 0.0, 1.0),
    Point3::new(-PHI, 0.0, 1.0),
    Point3::new(PHI, 0.0, -1.0),
    Point3::new(-PHI, 0.0, -1.0),

    Point3::new(1.0, PHI, 0.0),
    Point3::new(1.0, -PHI, 0.0),
    Point3::new(-1.0, PHI, 0.0),
    Point3::new(-1.0, -PHI, 0.0),
];

// vertex index -> array of vertex indices connected by an edge
pub const ICO_EDGE_INDEX: [[usize; 5]; 12] = [
    [2, 5, 10, 8, 4],
    [3, 6, 8, 10, 7],
    [0, 4, 9, 11, 5],
    [1, 7, 11, 9, 6],
    [0, 8, 6, 9, 2],
    [0, 2, 11, 7, 10],
    [1, 3, 9, 4, 8],
    [1, 10, 5, 11, 3],
    [0, 10, 1, 6, 4],
    [2, 4, 6, 3, 11],
    [0, 5, 7, 1, 8],
    [2, 9, 3, 7, 5],
];

// triplets of vertex indices which make a face
pub const ICO_FACES: [[usize; 3]; 20] = [
    [0, 4, 2],
    [0, 2, 5],
    [0, 8, 4],
    [0, 5, 10],
    [0, 10, 8],
    [1, 3, 6],
    [1, 7, 3],
    [1, 6, 8],
    [1, 10, 7],
    [1, 8, 10],
    [2, 4, 9],
    [2, 11, 5],
    [2, 9, 11],
    [3, 9, 6],
    [3, 7, 11],
    [3, 11, 9],
    [4, 8, 6],
    [4, 6, 9],
    [5, 7, 10],
    [5, 11, 7],
];

pub const ICO_NEIGHBORS: [([u8; 3], [usize; 3]); 20] = [
    ([2, 0, 0], [2, 10, 1]),
    ([2, 2, 0], [0, 11, 3]),
    ([2, 0, 0], [4, 16, 0]),
    ([2, 2, 0], [1, 18, 4]),
    ([2, 1, 0], [3, 9, 2]),
    ([2, 2, 0], [6, 13, 7]),
    ([2, 0, 0], [8, 14, 5]),
    ([2, 1, 0], [5, 16, 9]),
    ([2, 1, 0], [9, 18, 6]),
    ([2, 1, 0], [7, 4, 8]),
    ([1, 2, 0], [0, 17, 12]),
    ([2, 0, 1], [12, 19, 1]),
    ([2, 1, 0], [10, 15, 11]),
    ([2, 1, 1], [15, 17, 5]),
    ([1, 1, 0], [6, 19, 15]),
    ([2, 1, 0], [14, 12, 13]),
    ([1, 1, 0], [2, 7, 17]),
    ([2, 1, 1], [16, 13, 10]),
    ([2, 1, 1], [19, 8, 3]),
    ([1, 1, 0], [11, 14, 18]),
];

pub fn find_face_orientation(b: usize, a: usize) -> (u8, usize) {
    for (ii, &[oa, ob, oc]) in ICO_FACES.iter().enumerate() {
        if oa == a && ob == b {
            return (0, ii);
        }
        if ob == a && oc == b {
            return (1, ii);
        }
        if oc == a && oa == b {
            return (2, ii);
        }
    }
    unreachable!("neighbor must exist")
}

pub fn gen_face_neighbors() -> Vec<([u8; 3], [usize; 3])> {
    let mut neighbors = Vec::new();
    for &[a, b, c] in &ICO_FACES {
        let mut ori = [0; 3];
        let mut nei = [0; 3];
        (ori[0], nei[0]) = find_face_orientation(a, b);
        (ori[1], nei[1]) = find_face_orientation(b, c);
        (ori[2], nei[2]) = find_face_orientation(c, a);
        neighbors.push((ori, nei))
    }
    neighbors
}

fn smol(x: Float) -> bool {
    x.abs() < 1e-6
}

pub fn gen_ico_edge_index() -> Vec<Vec<usize>> {
    let count = ICO_VERTICES.len();
    let mut vertex_edges = vec![vec![]; count];
    for ii in 0..count - 1 {
        let a = ICO_VERTICES[ii];
        for jj in ii + 1..count {
            let b = ICO_VERTICES[jj];
            eprintln!("{} {} {}", ii, jj, (a-b).lensq());
            if smol((a-b).lensq() - 4.0) {
                vertex_edges[ii].push(jj);
                vertex_edges[jj].push(ii);
            }
        }
    }
    // eprintln!("{vertex_edges:?}");
    vertex_edges
}

pub fn gen_ico_faces() -> Vec<[usize; 3]> {
    let count = ICO_VERTICES.len();
    let mut faces = Vec::new();
    for ii in 0..count {
        let a = ICO_VERTICES[ii];
        for &jj in &ICO_EDGE_INDEX[ii] {
            if jj < ii {
                continue;
            }
            let b = ICO_VERTICES[jj];
            for &kk in &ICO_EDGE_INDEX[jj] {
                if kk < jj ||
                    !ICO_EDGE_INDEX[ii].contains(&kk)
                {
                    continue;
                }
                let c = ICO_VERTICES[kk];
                if (b - a).cross(c - a).dot(a) < 0.0 {
                    faces.push([ii, jj, kk]);
                } else {
                    faces.push([ii, kk, jj]);
                };
            }
        }
    }
    // eprintln!("{faces:?}");
    faces
}

pub fn gen_ico_sorted_edges() -> Vec<[usize; 5]> {
    let mut edges = Vec::new();
    for ii in 0..ICO_VERTICES.len() {
        let mut rv = [0; 5];
        let mut point: usize = ICO_EDGE_INDEX[ii][0];
        for jj in 0..5 {
            rv[jj] = point;
            let (ori, face) = find_face_orientation(point, ii);
            let index = (ori as usize + 2) % 3;
            point = ICO_FACES[face][index];
        }
        edges.push(rv);
    }
    edges
}
