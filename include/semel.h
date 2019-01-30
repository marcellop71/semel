#ifndef __semel__
#define __semel__

uint64_t _get_point_cloud_index(
    );

int32_t _add_point_cloud_raw(
    uint64_t k,
    uint32_t manifold, uint32_t distance,
    uint32_t n, uint32_t d,
    double * points,
    uint32_t nb, double * border);

int32_t _unload_point_cloud(
    uint64_t k);

int32_t _unload_point_cloud_all(
    );

int32_t _add_complex(
    uint64_t k, uint32_t fp, uint32_t fd);

int32_t _add_complex_ladic(
    uint64_t k, uint32_t l);

int32_t _unload_complex(
    uint64_t k);

int32_t _unload_complex_all(
    );

int32_t _add_simplices_raw(
    uint64_t k, uint8_t typeS,
    uint32_t n_simplices, uint32_t verts_per_simplex,
    uint64_t * indices);

int32_t _generate_delaunay(
    uint64_t k);

int32_t _build_complex(
    uint64_t k, uint32_t dim_sk, double rmax, uint8_t filtration_type);

char * _export_cocycles(
    uint64_t k);

char * _export_complex(
    uint64_t k);

char *  _export_coboundary(
    uint64_t k, uint32_t d);

int32_t _Laplace_eigen(
    uint64_t k, uint32_t n, uint32_t nev,
    uint32_t ** na, double *** a, double *** v);

char * _cohomology(
    uint64_t k, uint8_t measures_only);

int32_t _get_scales(
    uint64_t k,
    uint64_t * nr, uint64_t ** r, double * precision);

int32_t _add_graph_raw(
    uint64_t k,
    uint32_t n_edges, uint64_t * endpoints, double * weights);

int32_t _set_local_system(
    uint64_t k,
    uint32_t n_edges, uint64_t * endpoints, int64_t * weights);

int32_t _set_sheaf(
    uint64_t k, uint32_t rank,
    uint32_t n_edges, uint64_t * endpoints, int64_t * matrices);

char * _sheaf_cohomology(
    uint64_t k, uint8_t measures_only);

int32_t _sheaf_Laplace_eigen(
    uint64_t k, uint32_t n, uint32_t nev,
    uint32_t ** na, double *** a, double *** v);

int32_t sheaf_unroll(
    complex_cochain_t * Cc_orig, complex_cochain_t * Cc_u);

int32_t sheaf_unroll_free(
    complex_cochain_t * Cc_u);

int32_t _add_time_series_raw(
    uint64_t k, uint32_t d, uint32_t n, double * data);

int32_t _unload_time_series(
    uint64_t k);

int32_t _unload_time_series_all(
    );

int32_t _time_delay_embeddings(
    uint64_t k, uint32_t tau, uint32_t d,
    uint64_t * kP);

int32_t _nsrps(
    uint64_t ns, uint32_t * s);

int32_t persistent_cohomology(
    complex_cochain_t * Cc);

int32_t persistent_cohomology_z(
    complex_cochain_t * Cc);

uint32_t elementary_dimension(
    void * alpha, complex_cochain_t * C);

int32_t elementary_simplices_boundary(
    complex_cochain_t * Cc);

int32_t elementary_simplices_coboundary(
    complex_cochain_t * Cc);

int32_t elementary_simplices_coboundary_z(
    complex_cochain_t * Cc);

void * elementary_coboundary_of(
    uint64_t sid, complex_cochain_t * Cc);

void * elementary_coboundary_of_z(
    uint64_t sid, complex_cochain_t * Cc);

uint32_t elementary_delta(
    void * alpha, uint32_t method, complex_cochain_t * Cc,
    void ** d_alpha);

uint32_t elementary_delta_z(
    void * alpha, complex_cochain_t * Cc,
    void ** d_alpha);

int32_t coboundary_operator(
    complex_cochain_t * Cc);

int32_t coboundary_operator_as_gsl_matrix(
    uint32_t k, complex_cochain_t * Cc,
    uint32_t * mD, uint32_t * nD, gsl_matrix ** D);

int32_t build_cochain_complex(
    complex_cochain_t * Cc,
    semel_ctx_log_t * log);

int32_t add_faces(
    void * s, complex_simplicial_t * C, point_cloud_t * P, uint64_t rmax,
    void ** filtration);

int32_t complex_simplicial_topdown(
    complex_simplicial_t * C, point_cloud_t * P, uint64_t rmax);

int32_t add_cofaces(
    void * s, void * N, uint64_t rmax, void ** L, point_cloud_t * P,
    complex_simplicial_t * C);

int32_t expand_neighborhood_graph(
    uint32_t np, uint64_t rmax, void ** L, point_cloud_t * P,
    complex_simplicial_t * C);

int32_t annotated_lower_neighborhoods(
    uint64_t rmax, uint32_t nA, double * A,
    void *** L);

int32_t get_array_of_indices(
    void * a,
    uint32_t * nidx, uint64_t ** idx);

int32_t scales_from_distance_matrix(
    uint32_t nA, double * A,
    uint32_t * nr, uint64_t ** r);

int32_t complex_simplicial_bottomup(
    complex_simplicial_t * C, point_cloud_t * P, uint64_t rmax);

int32_t complex_simplicial_build(
    complex_simplicial_t * C);

int32_t complex_simplicial_build_auto(
    complex_simplicial_t * C,
    point_cloud_t * P,
    uint32_t dim_sk, double rmax,
    uint8_t filtration_type,
    semel_ctx_log_t * log);

char * export_cocycles(
    complex_cochain_t * Cc);

char * export_complex_simplicial(
    complex_simplicial_t * C);

char * export_cohomology(
    complex_cochain_t * Cc, uint8_t measures_only);

char * export_cohomology_z(
    complex_cochain_t * Cc, uint8_t measures_only);

char * export_cocycles_z(
    complex_cochain_t * Cc);

char * export_coboundary_operator(
    uint32_t k, complex_cochain_t * Cc);

double vector_dist(
    uint32_t d, double * p, double * q,
    uint32_t manifold, uint32_t distance, uint32_t nb, double * b);

double vector_norm_squared(
    uint32_t d, double * p);

double vector_dot(
    uint32_t d, double * p, double * q);

int32_t vector_add(
    uint32_t d, double * p, double * q, double * x);

int32_t vector_sub(
    uint32_t d, double * p, double * q, double * x);

int32_t vector_cross(
    uint32_t d, double * p, double * q, double * cp);

int32_t distance_matrix(
    uint32_t d, uint32_t n, double * p,
    uint32_t manifold, uint32_t distance, uint32_t nb, double * b,
    double ** A, double * Amax);

uint32_t num_points_from_card_of_distance_matrix(
    uint32_t a);

double get_dist(
    uint32_t i, uint32_t j, point_cloud_t * P);

double max_pairwise_dist_G(
    uint32_t d, uint32_t * pi, complex_simplicial_t * C);

double max_pairwise_dist(
    uint32_t d, uint32_t * pi, point_cloud_t * P);

double circumradius_jung(
    uint32_t d, uint32_t * pi, point_cloud_t * P);

double circumradius(
    uint32_t d, uint32_t * pi, point_cloud_t * P);

int32_t circumsphere(
    uint32_t dim, uint32_t np, double * p,
    double * cc, double * cr);

int32_t miniball(
    uint32_t dim,
    uint32_t np, double * p,
    void ** tau, void ** v,
    double * cc, double * cr);

int32_t smallest_ball(
    uint32_t dim, uint32_t np, double * p,
    double * cc, double * cr);

int32_t Laplace_eigen(
    complex_cochain_t * Cc, uint32_t n,
    uint32_t ** na, double *** a, double *** v);

int32_t Laplace_eigen_sparse(
    complex_cochain_t * Cc, uint32_t n, uint32_t nev,
    uint32_t ** na, double *** a, double *** v);

int32_t point_cloud_init_distance_matrix(
    point_cloud_t * P);

int32_t simplex_serialize(
    void * s, uint8_t * s_);

uint32_t simplex_dimension(
    void * s);

int32_t simplex_to_array(
    void * s, uint32_t * ns, point_cloud_t * P,
    uint32_t ** sa);

int32_t array_to_simplex(
    uint32_t ns, uint32_t * sa, void ** s);

double simplex_max_pairwise_dist_G(
    void * s, complex_simplicial_t * C);

double simplex_max_pairwise_dist(
    void * s, point_cloud_t * P);

double simplex_circumradius_jung(
    void * s, point_cloud_t * P);

double simplex_circumradius(
    void * s, point_cloud_t * P);

int32_t simplex_smallestball(
    void * s,
    uint32_t dim, uint32_t np, double * p,
    double * cc, double * cr);

uint64_t simplex_index_in_coface(
    void * s, void * s_coface);

int64_t simplex_sign_in_coface(
    void * s, void * s_coface);

int32_t simplex_faces(
    void * s, void ** d_s);

int32_t simplices_faces(
    complex_simplicial_t * C);

int32_t simplices_cofaces(
    complex_simplicial_t * C);

uint64_t semel_lap(
		void);

uint64_t bitreverse(
    uint64_t n);

uint32_t low(
    uint32_t nx, uint64_t * x);

int32_t identity(
    uint32_t n, uint64_t * V);

int32_t dense_compress(
    uint32_t n, uint64_t * A, uint64_t * B);

int32_t dense_decompress(
    uint32_t n, uint64_t * A, uint64_t * B);

int32_t dense_anti_transpose(
    uint32_t n, uint64_t * A, uint64_t * B);

int32_t dense_compressed_to_upper_triangle(
    uint32_t n, uint64_t * A, uint64_t * B);

int32_t nsrps(
    uint64_t nx, uint32_t * x,
    uint64_t * c, semel_ctx_log_t * log);

int32_t entropy_dist(
    uint64_t np, double * p, double alpha, double * e);

int32_t dist_words_seq_bin(
    uint64_t nx, uint8_t * x, uint64_t w,
    uint64_t * np, double ** p);

int32_t nrps_seq_bin(
    uint64_t nx, uint8_t * x);

int32_t _entropy_dist(
    uint64_t np, double * p, double alpha, double * e);

int32_t _dist_words_seq_bin(
    uint64_t nx, uint8_t * x, uint64_t w,
    uint64_t * np, double ** p);

int32_t _nrps_seq_bin(
    uint64_t nx, uint8_t * x);

int32_t ds_logistic(
    double lambda_,
    uint64_t np, double ** p);

int32_t Judy_free_1_a(
    void * ja);

int32_t Judy_free_2_a(
    void * ja);

int32_t Judy_free_1_b(
    void * ja, fq_ctx_t * fqctx);

int32_t Judy_free_2_b(
    void * ja, fq_ctx_t * fqctx);

int32_t Judy_free_1_z(
    void * ja);

int32_t Judy_free_2_z(
    void * ja);

fq_t * fq_one_(
    fq_ctx_t * fqctx);

fmpz_t * fmpz_one_(void);

int32_t uint64_to_256(
    uint64_t x, uint8_t * sx);

uint64_t uint256_to_64(
    uint8_t * sx);

#endif
