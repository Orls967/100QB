#ifndef MD_IO_H
#define MD_IO_H

#include "md_types.h"
#include <stdio.h>

/* Context helpers */
void md_ctx_init(md_run_ctx *ctx);
void md_ctx_record_stage(md_run_ctx *ctx, const char *id, const char *name,
                         double elapsed_sec, int files_written);
void md_ctx_record_file(md_run_ctx *ctx, const char *rel_path, const char *kind,
                        uint64_t hash, size_t bytes, int rows, int cols);
void md_ctx_add_criterion(md_run_ctx *ctx, const char *id, const char *status,
                          const char *measured, const char *threshold);
void md_ctx_add_flag(md_run_ctx *ctx, const char *stage, const char *file, const char *reason);

/* Compute FNV-1a hash and byte size of a file on disk */
bool md_hash_file(const char *path, uint64_t *hash_out, size_t *bytes_out);

/* Register a written CSV file: computes hash, counts rows/cols, runs automated vacuous audit */
void md_register_csv(md_run_ctx *ctx, const char *stage_id, const char *path);

/* Register a written SVG file: computes hash and validates basic XML well-formedness */
bool md_register_svg(md_run_ctx *ctx, const char *stage_id, const char *path);

/* Audit a single CSV file for vacuous patterns:
 * - non-parameter numeric columns that are constant across >= 5 rows
 * - NaN / Inf values
 * - unphysical negative values on strictly positive physical columns (temp, gr, msd, etc.)
 * - overly round suspicious constants across non-parameter outputs */
void md_audit_csv_vacuous(md_run_ctx *ctx, const char *stage_id, const char *path,
                          int *rows_out, int *cols_out);

/* Write manual JSON results file to ../../../outputs/day009a/results.json */
bool md_write_results_json(const md_run_ctx *ctx, const char *path);

#endif /* MD_IO_H */
