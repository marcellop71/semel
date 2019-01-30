#ifndef __enums__
#define __enums__

enum semel_manifold {
  MANIFOLD_FLAT_UNBOUNDED = 0,
  MANIFOLD_SURFACE_CYLINDER,
  MANIFOLD_SURFACE_TORUS,
};

enum semel_distance {
  DISTANCE_L2 = 0,
  DISTANCE_COSINE,
  DISTANCE_SURFACE_CYLINDER,
  DISTANCE_SURFACE_TORUS,
};

enum semel_cohomology_algorithm {
  COHOMOLOGY_FIELD = 0,
  COHOMOLOGY_INTEGRAL = 1,
};

enum semel_filtration {
  FILTRATION_ALPHA = 0,
  FILTRATION_RIPS,
};

//enum semel_complex_method {
//  COMPLEX_VIETORIS_RIPS,
//  COMPLEX_CECH,
//  COMPLEX_CECH_JUNG_APPROX,
//  COMPLEX_DELAUNAY_CECH,
//  COMPLEX_DELAUNAY_CECH_JUNG_APPROX,
//  COMPLEX_DELAUNAY,
//};

#endif
