#ifndef __ctx__
#define __ctx__

#include "semel.h"

struct semel_ctx_log {
  zlog_category_t * flow;
  zlog_category_t * ts;
  zlog_category_t * metrics;
  zlog_category_t * data;
};

struct semel_ctx {
  semel_ctx_log_t log;
  void * P; // pointclouds
  void * Cc; // complexes
  void * S; // timeseries
};

int32_t ctx_log_init(
    semel_ctx_t * c);

int32_t ctx_log_free(
    semel_ctx_t * c);

int32_t ctx_init(
    semel_ctx_t * c);

int32_t ctx_free(
    semel_ctx_t * c);

#endif
