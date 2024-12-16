use crate::{GreatCircle, GreatCircleArc, SphericalError, SphericalPoint, VEC_LEN_IS_ZERO};

/// Specifies the direction in which an edge is defined.
///
/// Flipping the direction will cause a polygon to be the complement of what you would expect to be your polygon (if you wanted to create a small triangle around the North Pole but flipped the edge orientation, you would define the polygon to be everywhere apart from the North Pole).
///
/// # Determining the direction
/// ## Algorithmic method
/// Imagine you are **inside** the sphere and are looking at the polygon. Now imagine a reference point inside the polygon and find the closest edge to it.
/// Now the direction is the same as the direction of the edge with reference to the point chosen.
///
/// ## More intuitive method
/// Imagine you are standing on the **inside** surface of the sphere, your head pointing in the direction of the centre of the sphere.
/// If you were to walk along the edge and the inside of the polygon was on your left choose `CounterClockwise`, else choose `Clockwise`.
#[derive(Clone, Copy)]
pub enum EdgeDirection {
    Clockwise,
    CounterClockwise,
}

/// A polygon on a unit sphere, given by its vertices and the edge direction
pub struct Polygon {
    vertices: Vec<SphericalPoint>,
    edges_direction: EdgeDirection,
}

impl Polygon {
    /// Creates a new polygon with the vertices and edges direction provided.
    ///
    /// # Important
    /// If the final vertex is not equal to the first one, it will be added automatically.
    ///
    /// Flipping the direction of the edges will cause a polygon to be the complement of what you would expect to be your polygon (if you wanted to create a small triangle around the North Pole but flipped the edge orientation, you would define the polygon to be everywhere apart from the North Pole).
    ///
    /// # Panics
    /// This function panics if no vertices were provided. That does not constitute a valid polygon and you should always provide at least two vertices.
    ///
    /// # Errors
    /// If any edge is defined by essentially equal or antipodal points, returns `SphericalError::AntipodalOrTooClosePoints` as in the case of identical or antipodal points the great circle (and therefore also the edge) is not uniquely defined.
    pub fn new(vertices_in: Vec<SphericalPoint>, edges_direction: EdgeDirection) -> Result<Self, SphericalError> {
        let mut vertices = vertices_in;
        if !vertices[0].approximately_equals(&vertices[vertices.len() - 1], VEC_LEN_IS_ZERO) {
            // The last vertex is not the same as the first one -> insert the first one to the back
            vertices.push(vertices[0]);
        }
        for i in 0..vertices.len() - 1 {
            if vertices[i].cartesian().cross(&vertices[i + 1].cartesian()).magnitude_squared() < VEC_LEN_IS_ZERO.powi(2) {
                return Err(SphericalError::AntipodalOrTooClosePoints);
            }
        }
        Ok(Self { vertices, edges_direction })
    }

    /// Returns a reference to the vertices list
    pub fn vertices(&self) -> &Vec<SphericalPoint> {
        &self.vertices
    }

    /// Returns the edges direction
    pub fn edges_direction(&self) -> EdgeDirection {
        self.edges_direction
    }

    /// Checks if the polygon contains the given point
    ///
    /// # Errors
    /// This function does not produce its own errors, but it will propagate inner errors out, see below. That should however never happen - if it does, it is a bug in an implementation in the library, so please report it should you encounter it.
    ///
    /// If any of the edges fails to be constructed as a `GreatCircleArc`, returns the corresponding error. This should however never happen, as that is checked when the polygon is constructed.
    ///
    /// Also, if any intersections fail the corresponding error will be returned. This should however also never happen.
    pub fn contains_point(&self, point: &SphericalPoint) -> Result<bool, SphericalError> {
        let tiebreaker_lim = 10e-5;
        // Algorithm description:
        // 1) Find the closest edge by finding an intersection with each of the edges with a great circle perpendicular to it. Use the clamped intersection, returning one of the endpoints in case of a miss.
        // 2) Determine if the closest edge is in the correct orientation

        // Step 1
        let mut closest_edge_i = 0;
        let mut closest_edge_dist_metric = f32::INFINITY;
        let mut tiebreaker = 0.0;
        for i in 0..self.vertices.len() - 1 {
            let edge = GreatCircleArc::new(self.vertices[i], self.vertices[i + 1])?;
            if edge.contains_point(point) {
                return Ok(true);
            }
            let (tiebreaker_dist, edge_distance_metric) = match edge.perpendicular_circle_through_point(point) {
                Ok(circle) => {
                    let closest_point = edge.closest_point_to_point_with_circle(&circle, point)?;
                    let unclamped_dist = GreatCircle::from_arc(&edge)
                        .intersect_great_circle(&circle)?
                        .iter()
                        .map(|p| p.minus_cotan_distance(point))
                        .min_by(|a, b| a.total_cmp(b))
                        .unwrap();
                    (unclamped_dist, closest_point.minus_cotan_distance(point))
                }
                Err(SphericalError::AntipodalOrTooClosePoints) => {
                    // The point is essentially the pole of the arc, so it is basically PI/2 radians away -> distance metric = -1/tan(PI/2) = 0
                    (f32::INFINITY, 0.0)
                }
                Err(err) => return Err(err),
            };
            if (edge_distance_metric - closest_edge_dist_metric).abs() < tiebreaker_lim {
                if tiebreaker < tiebreaker_dist {
                    closest_edge_i = i;
                    closest_edge_dist_metric = edge_distance_metric;
                    tiebreaker = tiebreaker_dist;
                }
            } else if edge_distance_metric < closest_edge_dist_metric {
                closest_edge_i = i;
                closest_edge_dist_metric = edge_distance_metric;
                tiebreaker = tiebreaker_dist;
            }
        }

        // Step 2
        let closest_edge_normal = self.vertices[closest_edge_i].cartesian().cross(&self.vertices[closest_edge_i + 1].cartesian());
        let cos_angle = closest_edge_normal.dot(&point.cartesian());
        let is_inside = match self.edges_direction {
            EdgeDirection::Clockwise => cos_angle >= 0.0,
            EdgeDirection::CounterClockwise => cos_angle <= 0.0,
        };

        Ok(is_inside)
    }

    /// Returns the intersections of the edges of the polygon with a given great circle arc
    ///
    /// # Errors
    /// This function does not generate its own errors, but may propagate the following:
    ///  - If any of the edges fails to be constructed as a [GreatCircleArc], returns the corresponding error (see [GreatCircleArc::new]). This should however never happen, as that is checked when the polygon is constructed.
    ///  - If an edge fails to be intersected with the arc, it returns the corresponding error (refer to [GreatCircleArc::intersect_great_circle_arc] for more details).
    pub fn great_circle_arc_intersections(&self, arc: &GreatCircleArc) -> Result<Vec<SphericalPoint>, SphericalError> {
        let mut intersections = Vec::new();

        for i in 0..self.vertices.len() - 1 {
            let edge = GreatCircleArc::new(self.vertices[i], self.vertices[i + 1])?;
            let ints = edge.intersect_great_circle_arc(arc)?;
            for int in ints {
                if intersections
                    .iter()
                    .any(|intersection: &SphericalPoint| intersection.approximately_equals(&int, crate::IDENTICAL_POINTS))
                {
                    continue;
                }
                intersections.push(int);
            }
        }

        Ok(intersections)
    }

    /// Checks if there exists an intersection of the edges of the polygon with the provided [GreatCircleArc]
    ///
    /// This function will in theory return errors less often than [Self::great_circle_arc_intersections] as it handles the cases when arcs are parallel
    ///
    /// # Errors
    /// This function does not generate its own errors, but may propagate the following:
    ///  - If any of the edges fails to be constructed as a [GreatCircleArc], returns the corresponding error (see [GreatCircleArc::new]). This should however never happen, as that is checked when the polygon is constructed.
    ///  - If an edge fails to be intersected with the arc, it returns the corresponding error (refer to [GreatCircleArc::intersect_great_circle_arc] for more details).
    pub fn intersects_great_circle_arc(&self, arc: &GreatCircleArc) -> Result<bool, SphericalError> {
        for i in 0..self.vertices.len() - 1 {
            let edge = GreatCircleArc::new(self.vertices[i], self.vertices[i + 1])?;
            if edge.intersects_great_circle_arc(arc)? {
                return Ok(true);
            }
        }
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    #[test]
    fn is_point_inside() {
        let polygon_1 = Polygon::new(
            vec![
                SphericalPoint::new(0.0, PI / 3.0),
                SphericalPoint::new(-2.0 * PI / 3.0, PI / 3.0),
                SphericalPoint::new(2.0 * PI / 3.0, PI / 3.0),
            ],
            EdgeDirection::CounterClockwise,
        )
        .expect("The polygon should be constructable");
        let north_pole = SphericalPoint::new(0.0, PI / 2.0);
        assert!(polygon_1.contains_point(&north_pole).expect("It should be possible to determine if the point is inside the polygon"));
        let south_pole = SphericalPoint::new(0.0, -PI / 2.0);
        assert!(!polygon_1.contains_point(&south_pole).expect("It should be possible to determine if the point is inside the polygon"));
    }

    #[test]
    fn intersects_arc() {
        let polygon_1 = Polygon::new(
            vec![
                SphericalPoint::new(0.0, PI / 3.0),
                SphericalPoint::new(-2.0 * PI / 3.0, PI / 3.0),
                SphericalPoint::new(2.0 * PI / 3.0, PI / 3.0),
            ],
            EdgeDirection::CounterClockwise,
        )
        .expect("The polygon should be constructable");

        let arc_1 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 3.0), SphericalPoint::new(-2.0 * PI / 3.0, PI / 3.0)).expect("The arc should be constructable");
        assert!(polygon_1
            .intersects_great_circle_arc(&arc_1)
            .expect("It should be possible to determine if the arc intersects the polygon"));

        let arc_2 = GreatCircleArc::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(0.0, PI / 3.0)).expect("The arc should be constructable");
        assert!(polygon_1
            .intersects_great_circle_arc(&arc_2)
            .expect("It should be possible to determine if the arc intersects the polygon"));

        let arc_3 = GreatCircleArc::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(0.0, PI / 2.0)).expect("The arc should be constructable");
        assert!(polygon_1
            .intersects_great_circle_arc(&arc_3)
            .expect("It should be possible to determine if the arc intersects the polygon"));

        let arc_4 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 4.0), SphericalPoint::new(PI / 6.0, PI / 5.0)).expect("The arc should be constructable");
        assert!(!polygon_1
            .intersects_great_circle_arc(&arc_4)
            .expect("It should be possible to determine if the arc intersects the polygon"));
    }

    #[test]
    fn arc_intersections() {
        let polygon_1 = Polygon::new(
            vec![
                SphericalPoint::new(0.0, PI / 3.0),
                SphericalPoint::new(-2.0 * PI / 3.0, PI / 3.0),
                SphericalPoint::new(2.0 * PI / 3.0, PI / 3.0),
            ],
            EdgeDirection::CounterClockwise,
        )
        .expect("The polygon should be constructable");

        let arc_1 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 3.0), SphericalPoint::new(-2.0 * PI / 3.0, PI / 3.0)).expect("The arc should be constructable");
        assert!(matches!(polygon_1.great_circle_arc_intersections(&arc_1), Err(SphericalError::IdenticalGreatCircles)));

        let arc_2 = GreatCircleArc::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(0.0, PI / 3.0)).expect("The arc should be constructable");
        assert_eq!(
            polygon_1
                .great_circle_arc_intersections(&arc_2)
                .expect("It should be possible to determine if the arc intersects the polygon")
                .len(),
            1
        );

        let arc_3 = GreatCircleArc::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(0.0, PI / 2.0)).expect("The arc should be constructable");
        assert_eq!(
            polygon_1
                .great_circle_arc_intersections(&arc_3)
                .expect("It should be possible to determine if the arc intersects the polygon")
                .len(),
            1
        );

        let arc_4 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 4.0), SphericalPoint::new(PI / 6.0, PI / 5.0)).expect("The arc should be constructable");
        assert!(polygon_1
            .great_circle_arc_intersections(&arc_4)
            .expect("It should be possible to determine if the arc intersects the polygon")
            .is_empty());

        let polygon_2 = Polygon::new(
            vec![
                SphericalPoint::new(0.0, 0.0),
                SphericalPoint::new(0.0, 0.5),
                SphericalPoint::new(1.2, 0.5),
                SphericalPoint::new(0.8, 0.25),
                SphericalPoint::new(1.2, 0.0),
            ],
            EdgeDirection::CounterClockwise,
        )
        .expect("The polygon should be constructable");

        let tolerance = 10e-4;
        let arc_5 = GreatCircleArc::new(SphericalPoint::new(1.0, 0.5), SphericalPoint::new(0.9, -0.15)).expect("The arc should be constructable");
        let intersections_2_5 = polygon_2
            .great_circle_arc_intersections(&arc_5)
            .expect("It should be possible to determine if the arc intersects the polygon");
        assert_eq!(intersections_2_5.len(), 3);
        dbg!(&intersections_2_5);
        let expected_i_1 = SphericalPoint::new(0.92165, 0.0);
        let expected_i_2 = SphericalPoint::new(0.94532, 0.16371);
        let expected_i_3 = SphericalPoint::new(0.97795, 0.37422);
        assert!(intersections_2_5.iter().any(|p| expected_i_1.approximately_equals(p, tolerance)));
        assert!(intersections_2_5.iter().any(|p| expected_i_2.approximately_equals(p, tolerance)));
        assert!(intersections_2_5.iter().any(|p| expected_i_3.approximately_equals(p, tolerance)));
    }
}
