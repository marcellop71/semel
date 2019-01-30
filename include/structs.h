#ifndef __structs__
#define __structs__

#include "semel.h"

struct point_cloud {
  uint32_t d; // dimension of embedding space
  uint32_t n; // number of points
  double * p; // point coordinates
  void * pn; // point names

  uint32_t n_;

  uint32_t manifold; // manifold
  uint32_t distance; // distance
  uint32_t nb;
  double * b;

  uint64_t nA; // number of entry in the distance matrix
  double * A; // distance matrix (strictly upper triangular)
  double Amax; // max distance
};

// simplex is JSL, value is const uint64_t = 1
// simplices are indexed by a simplexwise construction of the complex
// id is simplexwise order
struct complex_simplicial {
  void * S; // custom simplices: JL: index -> simplex
  uint8_t typeS;
  uint8_t filtration_type;

  void * G; // custom graph

  void * filtration; // JL: radius -> { JSL: serialization -> simplex }

  void * simplices; // JL: id -> simplex
  void * sx2id; // JSL: serialization -> id: value
  void * faces; // JL: id -> { JL: id -> face index: value }
  void * cofaces; // JL: id -> { JL: id -> coface index: value }

  void * dimension; // JL: id -> dimension: value
  void * age; // JL: id -> age: value
  void * radius; // JL: id -> radius: value

  //uint32_t nr;
  //uint64_t * r; // radius is an index (spatial and temporal), serialization is such that ordering makes the complex constructible

  uint32_t dim_sk; // dimensions
  void * scd [MAX_DIM]; // by dim: counter -> id
  void * sci [MAX_DIM]; // by dim: id -> counter
  uint64_t scc [MAX_DIM]; // by dim: simplices counter
  uint64_t scc_; // simplices counter
};

// chain/cochain is JL, value is fq_t
struct complex_cochain {
  complex_simplicial_t C;

  fq_ctx_t fqctx;

  void * D; // dim -> {}

  void * elementary_boundaries; // JL: id -> chain
  void * elementary_coboundaries; // JL: id -> cochain

  void * cocycles; // { cochain }
  void * intervals; // { a -> b }

  //void * cocycles_vanished; // { cochain }
  void * cocycles_harmonic;

  uint32_t algorithm;                  // COHOMOLOGY_FIELD or COHOMOLOGY_INTEGRAL
  uint32_t l_prime;                    // l-adic prime (0 = report all torsion)

  void * elementary_coboundaries_z;    // JL: id -> JL(id -> fmpz_t*)
  void * cocycles_z;                   // JL: id -> JL(id -> fmpz_t*)
  void * intervals_torsion;            // JL: interval_id -> uint64_t (torsion order)

  void * local_system;                 // JL: edge_key -> int64_t weight
                                       // edge_key = (min_vert << 32) | max_vert

  uint32_t sheaf_rank;                 // r (0 = no sheaf)
  void * sheaf_transport;              // JL: edge_key -> int64_t* (r*r row-major)
                                       // edge_key = (min_vert << 32) | max_vert
};

struct time_series {
  uint32_t d; // dimension of x = (x_0,...,x_{dx-1})
  uint32_t n; // length of x
  double * x; // data
};

int32_t point_cloud_init(
    point_cloud_t * P,
    uint32_t d, uint32_t n,
    uint32_t manifold, uint32_t distance,
    uint32_t nb, double * b);

int32_t point_cloud_free(
    point_cloud_t * P);

int32_t complex_simplicial_init(
    complex_simplicial_t * C);

int32_t complex_simplicial_free(
    complex_simplicial_t * C);

int32_t complex_cochain_init(
    complex_cochain_t * Cc, uint32_t fp, uint32_t fd);

int32_t complex_cochain_init_z(
    complex_cochain_t * Cc, uint32_t l);

int32_t complex_cochain_free(
    complex_cochain_t * Cc);

int32_t time_series_init(
    time_series_t * S,
    uint32_t d, uint32_t n);

int32_t time_series_free(
    time_series_t * S);

#endif
