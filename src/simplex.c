#include "common.h"
#include "semel.h"

int32_t simplex_to_array(
    void * s, uint32_t * ns, point_cloud_t * P,
    uint32_t ** sa)
{
  uint64_t * p, * q;
  uint64_t i;
  uint8_t v [C_1_KB];

  uint32_t d = simplex_dimension(s);
  *ns = d + 1;
  *sa = (uint32_t *) malloc((*ns) * sizeof(uint32_t));

  v[0] = 0; JSLF(p, s, v); i = 0;
  while (p != NULL)
  {
    if (P != NULL)
    {
      JSLG(q, P->pn, v);
      if (q != NULL)
        { (*sa)[i] = (uint32_t) *q; }
      else
        { (*sa)[i] = 0; }
    }
    else
      { (*sa)[i] = (uint32_t) uint256_to_64(&v[0]); }

    JSLN(p, s, v); i++;
  }

  return 0;
}

uint32_t simplex_dimension(
    void * s)
{
  uint64_t c;
  uint64_t * p;
  uint8_t v [C_1_KB];
  v[0] = 0; JSLF(p, s, v); c = 0;
  while (p != NULL)
    { JSLN(p, s, v); c++; }
  return ((uint32_t) (c - 1));
}

uint64_t simplex_id(
    void * s, complex_simplicial_t * C)
{
  uint64_t * p;
  uint8_t sx [MAX_LEN_SIMPLEX_SERIALIZATION];
  simplex_serialize(s, &sx[0]);
  JSLG(p, C->sx2id, sx);
  uint64_t id = (uint64_t) -1;
  if (p != NULL)
    { id = *p; }
  return id;
}

double simplex_max_pairwise_dist_G(
    void * s, complex_simplicial_t * C)
{
  uint32_t ns = 0; uint32_t * sa = NULL;
  simplex_to_array(s, &ns, NULL, &sa);
  double r = max_pairwise_dist_G(ns - 1, sa, C);
  free(sa);
  return r;
}

double simplex_max_pairwise_dist(
    void * s, point_cloud_t * P)
{
  uint32_t ns = 0; uint32_t * sa = NULL;
  simplex_to_array(s, &ns, P, &sa);
  double r = max_pairwise_dist(ns - 1, sa, P);
  free(sa);
  return r;
}

double simplex_circumradius(
    void * s, point_cloud_t * P)
{
  uint32_t ns = 0; uint32_t * sa = NULL;
  simplex_to_array(s, &ns, P, &sa);
  double r = circumradius(ns - 1, sa, P);
  free(sa);
  return r;
}

uint64_t simplex_index_in_coface(
    void * s, void * s_coface)
{
  uint64_t * pt, * ps;
  uint8_t v [C_1_KB];
  uint64_t i;

  v[0] = 0; JSLF(pt, s_coface, v); i = 0;
  while (pt != NULL)
  {
    JSLG(ps, s, v);
    if (ps == NULL)
      { break; }
    JSLN(pt, s_coface, v); i++;
  }

  return i;
}

int64_t simplex_sign_in_coface(
    void * s, void * s_coface)
{
  uint64_t i = simplex_index_in_coface(s, s_coface);
  int64_t sign = ((i % 2) == 0) ? 1 : -1;
  return sign;
}

// vertices in the simplex are (totally) ordered
// to each face (a simplex) is attached the actual id of the vertex excluded
// (its ordinal, i.e. the ordinal of the vertex which is excluded,
// is already tracked in the index)

int32_t simplex_faces(
    void * s, void ** s_faces)
{
  uint64_t * p, * q;
  uint64_t i, v_k;
  uint8_t v [C_1_KB];

  uint64_t s_dim = simplex_dimension(s);
  for(uint64_t k = 0; k <= s_dim; k++)
  {
    v_k = 0;
    //JLBC(p, s, k+1, v);
    //if (p != NULL)
    //  { v_k = *p; }

    void * face = NULL;

    v[0] = 0; JSLF(p, s, v); i = 0;
    while (p != NULL)
    {
      if (i != k)
        { JSLI(q, face, v); *q = v_k; }
      JSLN(p, s, v); i++;
    }

    if (face != NULL)
      { JLI(p, *s_faces, k); *p = (uint64_t) face; }
  }

  return 0;
}

int32_t simplices_faces(
    complex_simplicial_t * C)
{
  uint64_t v, k, id_face;
  uint64_t * v_face, * p, * q;
  int rc;

  v = 0; JLF(p, C->simplices, v);
  while (p != NULL)
  {
    void * v_faces = NULL;
    simplex_faces(*p, &v_faces);

    void * s = NULL;
    k = 0; JLF(v_face, v_faces, k);
    while (v_face != NULL)
    {
      id_face = simplex_id(*v_face, C);
      if (id_face != (uint64_t) -1)
        { JLI(q, s, id_face); *q = k; }

      JLN(v_face, v_faces, k);
    }

    k = 0; JLF(q, v_faces, k);
    while (q != NULL)
    {
      JSLFA(rc, *q);
      JLN(q, v_faces, k);
    }
    JLFA(rc, v_faces);

    JLI(q, C->faces, v); *q = (uint64_t) s;
    JLN(p, C->simplices, v);
  }

  return 0;
}

int32_t simplices_cofaces(
    complex_simplicial_t * C)
{
  uint64_t v, v_face;
  uint64_t * k, * p, * q, * q_;

  v = 0; JLF(p, C->faces, v);
  while (p != NULL)
  {
    v_face = 0; JLF(k, *p, v_face);
    while (k != NULL)
    {
      JLG(q, C->cofaces, v_face);
      if (q == NULL)
      {
        void * s = NULL;
        JLI(q_, s, v); *q_ = *k;
        JLI(q_, C->cofaces, v_face); *q_ = (uint64_t) s;
      }
      else
        { JLI(q_, *q, v); *q_ = *k; }

      JLN(k, *p, v_face);
    }

    JLN(p, C->faces, v);
  }

  return 0;
}
