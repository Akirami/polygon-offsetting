use std::collections::HashMap;

#[derive(Default, Debug, Clone)]
pub struct Offsets {
    pub contour: Vec<(f64, f64)>,
    pub offset_area: f64,
    pub perimeter: f64,
}

#[derive(Default, Debug, Copy, Clone)]
pub struct Vertex {
    pub x: f64,
    pub y: f64,
    pub is_intersect: bool
}

#[derive(Default, Debug, Copy, Clone)]
pub struct Edge {
    pub p1: usize,
    pub p2: usize,
    pub index: usize,
    pub outward_normal: Vertex,
}

#[derive(Default, Debug, Clone)]
pub struct Polygon {
    pub edges: Vec<Edge>,
    pub vertices: HashMap<usize, Vertex>,
    pub offset_margin: f64,
}

#[derive(Default, Clone, Debug)]
pub struct Segment {
    pub p1: (f64, f64),
    pub p2: (f64, f64),
}

#[inline]
pub fn compute_area(contours: &Vec<(f64, f64)>) -> f64 {
    let mut a = 0.0;
    if contours.len() == 0 { return 0.0 }
    for i in 0..contours.len() - 1 {
        a = a + (contours[i].0 * contours[i + 1].1) - (contours[i + 1].0 * contours[i].1);
    }
    (a * -0.5).abs()
}

pub fn compute_perimeter(contour2d: &Vec<(f64, f64)>) -> f64 {
    let mut perimeter2d = 0.;
    for i in 0..(contour2d.len() - 1) {
        let p1 = &contour2d[i];
        let p2 = &contour2d[i + 1];
        perimeter2d += ((p2.0 - p1.0).powi(2) + (p2.1 - p1.1).powi(2)).sqrt();
    }
    let p1 = &contour2d[contour2d.len() - 1];
    let p2 = &contour2d[0];
    perimeter2d += ((p2.0 - p1.0).powi(2) + (p2.1 - p1.1).powi(2)).sqrt();
    perimeter2d
}

#[inline]
pub fn get_dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
    ((p2.0 - p1.0).powi(2) + (p2.1 - p1.1).powi(2)).sqrt()
}

#[inline]
pub fn vector_sub(v1: (f64, f64), v2: (f64, f64)) -> (f64, f64) {
    (v1.0 - v2.0, v1.1 - v2.1)
}

#[inline]
pub fn vector_add(v1: (f64, f64), v2: (f64, f64)) -> (f64, f64) {
    (v1.0 + v2.0, v1.1 + v2.1)
}

#[inline]
pub fn vector_scale(v1: (f64, f64), len: f64) -> (f64, f64) {
    (v1.0 * len , v1.1 * len)
}

#[inline]
pub fn vector_dot(v1: (f64, f64), v2: (f64, f64)) -> f64 {
    (v1.0 * v2.0) + (v1.1 * v2.1)
}

#[inline]
pub fn normalize_vec(v1: (f64, f64)) -> (f64, f64) {
    let len: f64 = f64::sqrt(vector_dot(v1, v1));
    if len == 0.0 {
        return (v1.0, v1.1)
    };
    (v1.0 / len, v1.1 / len)
}

#[inline]
pub fn reverse_segments(sgmts: &Vec<Segment>) -> Vec<Segment> {
    let mut segments: Vec<Segment> = Vec::new();

    for s in sgmts.iter().rev() {
        segments.push( Segment {
            p1: s.p2,
            p2: s.p1,
        });
    }
    segments
}

// ================================================================================= 

impl Polygon {
    fn append_arc(&self,  center: &Vertex, radius: f64, start_vertex: &Vertex,  end_vertex: &Vertex, tolerance: f64) -> Vec<Vertex> {        
        // we compute start and end angles
        let mut start_angle = (start_vertex.y - center.y).atan2(start_vertex.x - center.x);
        let mut end_angle   = (end_vertex.y   - center.y).atan2(end_vertex.x   - center.x);
        if start_angle < 0.  { start_angle = 2. * std::f64::consts::PI + start_angle } 
        if end_angle   <= 0. { end_angle   = 2. * std::f64::consts::PI + end_angle } 

        // we compute oriented angle
        let mut angle : f64 = end_angle - start_angle;
        if angle.abs() > std::f64::consts::PI { 
            angle = -( angle / angle.abs()) *  (2. * std::f64::consts::PI - angle.abs());
        }

        // we compute number of segments to fit the tolerance
        let n_angle = (1. - (4. * radius * tolerance - 2. * tolerance.powi(2) ) / (radius.powi(2)) ).acos();
        if n_angle == 0.0 {
            let vect = vector_add((start_vertex.x, start_vertex.y), vector_sub((end_vertex.x, end_vertex.y), (start_vertex.x, start_vertex.y)));
            return vec![Vertex { x: vect.0, y: vect.1, is_intersect: false }]
        }
        let nseg: i64 = ((angle.abs() / n_angle.abs()).round() + 1.) as i64;

        // we define the angular step
        let angular_step = angle / nseg as f64;
        
        // we create the segments
        let mut dots = Vec::new();
        for i in 0..(nseg + 1) {
            let theta = start_angle + angular_step * i as f64;
            let x = radius * (theta).cos() + center.x;
            let y = radius * (theta).sin() + center.y;
            dots.push(Vertex {x: x, y: y, is_intersect: false});
        }
        dots
    }

    fn edges_intersection(&self, e1: &(Vertex, Vertex), e2: &(Vertex, Vertex), is_inters: bool) -> Option<Vertex> {
        let den = (e2.1.y - e2.0.y) * (e1.1.x - e1.0.x) - (e2.1.x - e2.0.x) * (e1.1.y - e1.0.y);
        if den > -0.0001 && den < 0.0001 {
            return None;  // lines are parallel or conincident
        }

        let ua = ((e2.1.x - e2.0.x) * (e1.0.y - e2.0.y) - (e2.1.y - e2.0.y) * (e1.0.x - e2.0.x)) / den;
        let ub = ((e1.1.x - e1.0.x) * (e1.0.y - e2.0.y) - (e1.1.y - e1.0.y) * (e1.0.x - e2.0.x)) / den;

        if ua < 0.0000001 || ub < 0.0000001 || ua > 0.9999999 || ub > 0.9999999 {
            return None;
        }

        let v_cross = Some(Vertex {
            x: e1.0.x + ua * (e1.1.x - e1.0.x),
            y: e1.0.y + ua * (e1.1.y - e1.0.y),
            is_intersect: is_inters
        });
        v_cross
    }

    fn create_offset_edge(&self, p1: &Vertex, p2: &Vertex, dx: f64, dy: f64) -> (Vertex, Vertex) {
        let mut v1: Vertex = Vertex::default();
        let mut v2: Vertex = Vertex::default();
        v1.x = p1.x + dx;
        v1.y = p1.y + dy;
        v2.x = p2.x + dx;
        v2.y = p2.y + dy;
        (v1, v2)
    }

    fn create_margin_polygon(&mut self, tolerance: f64) -> Polygon {
        let mut offset_edges: Vec<(Vertex, Vertex)> = Vec::new();
        let mut vertices: HashMap<usize, Vertex> = HashMap::new();
        let mut index: usize = 0;

        // compute and store offsets points
        self.edges.iter().for_each(|edge| {
            let p1 = self.vertices.get(&edge.p1).unwrap();
            let p2 = self.vertices.get(&edge.p2).unwrap();
            let dx = edge.outward_normal.x * self.offset_margin;
            let dy = edge.outward_normal.y * self.offset_margin;
            offset_edges.push(self.create_offset_edge(p1, p2, dx, dy));
        });

        
        for i in 0..offset_edges.len() {
            let this_edge = &offset_edges[i];
            let prev_edge = &offset_edges[(i + offset_edges.len() - 1) % offset_edges.len()];
            let vertex = self.edges_intersection(prev_edge, this_edge, false);

            // Check if offsets segments intersect or not
            if vertex.is_none() {
                let arc_center = &self.vertices.get(&i).unwrap();
                let arc_vertices = &self.append_arc(arc_center, self.offset_margin, &prev_edge.1, &this_edge.0, tolerance);

                // if not, adding arc segments, welds & boundaries matching our welded / bounded edges

                // adding index and vertex to our vertices container
                for av in arc_vertices.iter() {
                    vertices.insert(index, *av);
                    index += 1;
                }
            } else {
                vertices.insert(index, vertex.unwrap());
                index += 1;
            }
        }

        // Finally, recreate a polygon from vertices, welds ans boundaries
        let margin_polygon = self.create_polygon(vertices, self.offset_margin, false);
        margin_polygon
    }

    fn sort_by_squared_dist(&self, p1: &Vertex, cross: &Vec<usize>, poly: &Polygon) -> Vec<usize> {
        let mut cross_coord: Vec<(usize, Vertex)> = Vec::new();
        let mut vp: Vec<(f64, usize)> = Vec::new();
        let mut crossp = Vec::new();

        for c in cross {
            cross_coord.push((*c, *poly.vertices.get(&c).unwrap()));
        }

        for i in 0..cross_coord.len() {
            let dist = (p1.x - cross_coord[i].1.x).powi(2) + (p1.y - cross_coord[i].1.y).powi(2);
            vp.push(( dist, cross_coord[i].0));
        }

        vp.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        vp.iter().for_each(|v| {
            crossp.push(v.1);
        });
        crossp
    }

    fn detect_all_intersect(&mut self, margin_polygon: &mut Polygon) {
        let mut poly: Polygon = Polygon::default();
        let mut indices: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut vertices = margin_polygon.vertices.clone();
        
        let mut iteration = 0;
        margin_polygon.edges.iter_mut().for_each(|edge| {
            indices.insert(edge.index, Vec::new());
            edge.index = iteration;
            iteration += 1;
        });

        let mut max_vertices = 0;
        margin_polygon.vertices.iter().for_each(|(id, _)| {
            if id > &max_vertices { max_vertices = *id;}
        });

        // Store intersection points with index 
        margin_polygon.edges.iter().for_each(|edge1| {
            margin_polygon.edges.iter().for_each(|edge2| {
                if edge1.index > edge2.index {
                    let p1: Vertex = *margin_polygon.vertices.get(&edge1.p1).unwrap();
                    let p2: Vertex = *margin_polygon.vertices.get(&edge1.p2).unwrap();
                    let p3: Vertex = *margin_polygon.vertices.get(&edge2.p1).unwrap();
                    let p4: Vertex = *margin_polygon.vertices.get(&edge2.p2).unwrap();
                    let e1: (Vertex, Vertex) = (p1, p2);
                    let e2: (Vertex, Vertex) = (p3, p4);

                    let inters = self.edges_intersection(&e1, &e2, true);
                    if !inters.is_none() {
                        vertices.insert(max_vertices + 1, inters.unwrap());
                        max_vertices += 1;
                        
                        match indices.get_mut(&edge1.index) {
                            Some(v) => {
                                if !v.contains(&max_vertices) {
                                    v.push(max_vertices);
                                }
                            },
                            _ => {}
                        };
                        match indices.get_mut(&edge2.index) {
                            Some(v) => {
                                if !v.contains(&max_vertices) {
                                    v.push(max_vertices);
                                }
                            },
                            _ => {}
                        };
                    }
                }
            });
        });

        // Split segments into 2 new segments from intersection point
        margin_polygon.vertices = vertices.clone();
        let mut new_poly: Vec<Edge> = Vec::new();
        margin_polygon.edges.iter().for_each(|edge| {
            let idx = edge.index;
            let p1 = margin_polygon.vertices.get(&edge.p1).unwrap();

            match indices.get(&idx) {
                Some(cross) => {
                    if cross.len() > 0 {
                        let sorted_cross: Vec<usize> = self.sort_by_squared_dist(p1, &cross.to_vec(), &margin_polygon);
                        let new_edge: Edge = Edge {
                            p1: edge.p1,
                            p2: sorted_cross[0],
                            index: 0,
                            outward_normal: edge.outward_normal
                        };
                        new_poly.push(new_edge);

                        for i in 0..sorted_cross.len() - 1 {
                            let new_edge: Edge = Edge {
                                p1: sorted_cross[i],
                                p2: sorted_cross[i + 1],
                                index: 0,
                                outward_normal: edge.outward_normal
                            };
                            new_poly.push(new_edge);
                        }
                        let new_edge: Edge = Edge {
                            p1: sorted_cross[sorted_cross.len() - 1],
                            p2: edge.p2,
                            index: 0,
                            outward_normal: edge.outward_normal
                        };
                        new_poly.push(new_edge);
                    } else {
                        new_poly.push(*edge);
                    }
                },
                _ => {}
            }
        });
        let mut iteration = 0;
        new_poly.iter_mut().for_each(|edge| {
            edge.index = iteration;
            iteration += 1;
        });
        poly.edges = new_poly;
        poly.vertices = vertices;
        *margin_polygon = poly;
    }

    fn detect_regions(&self, polygon: &Polygon) -> Vec<Polygon> {
        // each region is described by a polygon
        let mut regions : Vec<Polygon> = Vec::new();

        // remaining is a vect of the indices of the vertices
        let mut remaining: Vec<usize> = Vec::new();

        polygon.vertices.iter().for_each(|(id, _)| remaining.push(*id));

        // we create a map: keys are the vertice index, values are the list of edges containing this vertice
        let mut map : HashMap<usize, Vec<usize>> = HashMap::new();
        polygon.vertices.iter().for_each(|(id, _)| { map.insert(*id, Vec::new()); });
        polygon.edges.iter().for_each(|edge| {
            match map.get_mut(&edge.p1) { Some(v) => {if !v.contains(&edge.index) {v.push(edge.index)}}, _ => {} };
        });

        // we follow the contour
        while remaining.len() > 0 {
            
            let mut current_region: Polygon = Polygon {
                edges: Vec::new(),
                vertices: HashMap::new(),
                offset_margin: 0.0,
            };
            
            let mut idx = remaining[0];
            let start_idx = idx;
            let mut prev_edge_indice: usize = 0;

            while idx != start_idx || current_region.vertices.len() == 0 {
                match polygon.vertices.get(&idx) {
                    Some(v) => {
                        // we add the current point to the vertices of the current region
                        current_region.vertices.insert(idx, *v);
                        remaining.retain(|r| {
                            *r != idx
                        });

                        match v.is_intersect {
                            true => {
                                let mut edges: Vec<usize> = map.get(&idx).unwrap().clone();
                                edges.retain(|e| {
                                    *e != prev_edge_indice + 1
                                });
                                let edge = edges[0];
                                let new_edge = Edge { 
                                    p1: idx, 
                                    p2: polygon.edges[edge].p2,
                                    outward_normal: polygon.edges[edge].outward_normal,
                                    index: 0
                                };
                                current_region.edges.push(new_edge);
                                idx = polygon.edges[edge].p2;
                                prev_edge_indice = edge;
                            },
                            false => {
                                let edge = map.get(&idx).unwrap()[0];
                                let new_edge = Edge { 
                                    p1: idx, 
                                    p2: polygon.edges[edge].p2,
                                    outward_normal: polygon.edges[edge].outward_normal,
                                    index: 0
                                };
                                current_region.edges.push(new_edge);
                                idx = polygon.edges[edge].p2;
                                prev_edge_indice = edge;
                            }
                        }
                    },
                    _ => panic!("Regions, vertice not found")
                }
            }
            let mut iteration = 0;
            current_region.edges.iter_mut().for_each(|e| {
                e.index = iteration;
                iteration += 1;
            });
            regions.push(current_region)
        }
        regions
    }

    fn outward_edge_normal(&self, v1: &Vertex, v2: &Vertex) -> Vertex {
        let dx = v2.x - v1.x;
        let dy = v2.y - v1.y;
        let edge_length = (dx * dx + dy * dy).sqrt();
        Vertex {
            is_intersect: false,
            x:  dy / edge_length, 
            y: -dx / edge_length
        }
    }

    fn contour_to_vertices(contour: Vec<(f64, f64)>) -> HashMap<usize, Vertex> {
         let mut vtxs: HashMap<usize, Vertex> = HashMap::new();

         for i in 0..contour.len() - 1 {
             let mut vertex: Vertex = Vertex::default();
             vertex.is_intersect = false;
             vertex.x = contour[i].0;
             vertex.y = contour[i].1;
             vtxs.insert(i, vertex);
         }
         vtxs
    }

    pub fn get_polygon_area<'a>(&self, poly: &'a mut Polygon) -> (f64, &'a Polygon) {
        let mut contours: Vec<(f64, f64)> = Vec::new();

        let edges = &mut poly.edges;
        edges.sort_by_key(|k| k.index);

        for i in 0..edges.len() {
            let v = poly.vertices.get(&edges[i].p1).unwrap();
            contours.push((v.x, v.y));
        }

        let mut a = 0.0;
        for i in 0..contours.len() - 1 {
            a = a + (contours[i].0 * contours[i + 1].1) - (contours[i + 1].0 * contours[i].1);
        }
        a = a + (contours[contours.len() - 1].0 * contours[0].1) - (contours[0].0 * contours[contours.len() - 1].1);
        ((a * -0.5).abs(), poly)
    }

    // convert tuples to a Vec of Segment
    pub fn tuples_to_segments(contour1: &Vec<(f64, f64)>) -> Vec<Segment> {
        let mut segments: Vec<Segment> = Vec::new();

        for i in 1..contour1.len() {
            let mut sgmt: Segment = Segment::default();
            sgmt.p1.0 = contour1[i - 1].0;
            sgmt.p1.1 = contour1[i - 1].1;
            if i == contour1.len() { 
                sgmt.p2.0 = contour1[0].0;
                sgmt.p2.1 = contour1[0].1;
                segments.push(sgmt);
                break ;
            }
            sgmt.p2.0 = contour1[i].0;
            sgmt.p2.1 = contour1[i].1;
            segments.push(sgmt);
        }
        segments
    }

    fn create_polygon(
        &mut self,
        vertices: HashMap<usize, Vertex>,
        offset_size: f64,
        is_initial_polygon: bool
    ) -> Polygon {
        let mut polygon = Polygon::default();
        let mut edges: Vec<Edge> = Vec::new();
        polygon.vertices = vertices;
        polygon.offset_margin = offset_size;

        // Create edges segments with their index and normal
        for i in 0..polygon.vertices.len() {
            let mut edge = Edge {
                p1: i, 
                p2: (i + 1) % polygon.vertices.len(),
                index: i,
                outward_normal: Vertex::default(),
            };
            if edge.p1 == edge.p2 {
                continue ;
            }
            if is_initial_polygon {
                edge.outward_normal = self.outward_edge_normal(
                    polygon.vertices.get(&edge.p1).unwrap(),
                    polygon.vertices.get(&edge.p2).unwrap()
                );
            }
            edges.push(edge);
        }

        polygon.edges = edges;
        polygon
    }

    pub fn new(
        initial_contour: &Vec<(f64, f64)>,
        offset_size: f64,
    ) -> Result <Polygon, Box<dyn std::error::Error>> {

        let mut sub_blank_segments: Vec<Segment> = Polygon::tuples_to_segments(&initial_contour);
        sub_blank_segments.retain(|s| {
            get_dist(s.p1, s.p2) != 0.
        });

        let value = sub_blank_segments.iter().fold(
            0.,
            |acc, seg|
            acc + (seg.p2.0 - seg.p1.0) * (seg.p2.1 + seg.p1.1)
        );

        if value > 0. { sub_blank_segments = reverse_segments(&sub_blank_segments); }

        // compute extra margin for our welds
        let mut segments: Vec<Segment> = Vec::new();
        for sgmt in sub_blank_segments.iter_mut() {
            segments.push( Segment { p1: sgmt.p1, p2: sgmt.p2 });
        }

        // Offset 0.0 mean we don't compute offseting
        if offset_size == 0.0 {
            segments.clear();
            sub_blank_segments.iter().for_each(|sgmt| {
                segments.push( Segment { p1: sgmt.p1, p2: sgmt.p2 });
            });
        }

        let mut points: Vec<(f64, f64)> = Vec::new();
        for s in segments.iter() {
            points.push(s.p1);
        }
        points.push(points[0]);

        let mut initial_polygon: Polygon = Polygon::default();
        let vertices: HashMap<usize, Vertex> = Polygon::contour_to_vertices(points.to_vec());
        Ok(initial_polygon.create_polygon(vertices, offset_size, true))
    }

    pub fn offsetting(&mut self, tolerance: f64) -> Result <Offsets, Box<dyn std::error::Error>> {
        let mut offsets: Offsets = Offsets::default();
        let mut areas: Vec<(f64, &Polygon)> = Vec::new();

        if tolerance <= 0.0 {
            return Err("The tolerance can't be below or egal to 0".into())
        }

        // Offset 0.0, only return an Offsets struct without computing
        if self.offset_margin == 0.0 {
            let mut points: Vec<(f64, f64)> = Vec::new();
            self.edges.iter().for_each(|e| {
                let p1 = self.vertices.get(&e.p1).unwrap();
                points.push((p1.x, p1.y));
            });
            points.push(points[0]);

            let offset = Offsets {
                offset_area: compute_area(&points),
                perimeter: compute_perimeter(&points),
                contour: points,
            };
            return Ok(offset)
        }

        // Now we can compute our margin polygon
        let mut margin_polygon: Polygon = self.create_margin_polygon(tolerance);
        
        // Check for all intersections
        self.detect_all_intersect(&mut margin_polygon);

        // Intersections for new regions
        let mut regions = self.detect_regions(&margin_polygon);
        regions.iter_mut().for_each(|r| {
            areas.push(self.get_polygon_area(&mut *r));
        });
        areas.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        
        // And keep only the biggest
        let offset_polygon: &Polygon = areas[areas.len() - 1].1;
        offsets.offset_area = areas[areas.len() - 1].0;
              
        offset_polygon.edges.iter().for_each(|edge| {
            offsets.contour.push((
                offset_polygon.vertices.get(&edge.p1).unwrap().x,
                offset_polygon.vertices.get(&edge.p1).unwrap().y
            ));
        });
        offsets.contour.push(offsets.contour[0]); 
        offsets.perimeter = compute_perimeter(&offsets.contour);

        Ok(offsets)
    }
}






