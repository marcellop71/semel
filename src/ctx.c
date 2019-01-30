#include "common.h"

int32_t ctx_log_init(
    semel_ctx_t * ctx)
{
  const char * zlog_conf = getenv("SEMEL_ZLOG_CONF");
  if (zlog_conf == NULL) zlog_conf = ZLOG_CONF_PATH;
  int rc = zlog_init(zlog_conf);
  if (rc) { fprintf(stderr, "zlog init failed\n"); return -1; }
  ctx->log.flow = zlog_get_category("flow");
  ctx->log.ts = zlog_get_category("ts");
  ctx->log.metrics = zlog_get_category("metrics");
  ctx->log.data = zlog_get_category("data");

  //zlog_info(ctx->log.flow, "log started");

  return 0;
}

int32_t ctx_log_free(
    semel_ctx_t * ctx)
{
  //zlog_info(ctx->log.flow, "log end");

  zlog_fini();

  return 0;
}

int32_t ctx_init(
    semel_ctx_t * ctx)
{
  //zlog_info(ctx->log.flow, "init context");

  ctx->P = NULL;
  ctx->Cc = NULL;
  ctx->S = NULL;

  return 0;
}

int32_t ctx_free(
    semel_ctx_t * ctx)
{
  int rc;
  JLFA(rc, ctx->P);
  JLFA(rc, ctx->Cc);
  JLFA(rc, ctx->S);

  return 0;
}
