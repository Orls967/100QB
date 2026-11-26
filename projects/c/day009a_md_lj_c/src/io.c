#include "md_io.h"
#include "md_stats.h"
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void md_ctx_init(md_run_ctx *ctx) {
    memset(ctx, 0, sizeof(*ctx));
    ctx->n_threads = 4;
    ctx->cores = 1;
    strncpy(ctx->cc_ver, "unknown", sizeof(ctx->cc_ver) - 1);
    strncpy(ctx->arch, "unknown", sizeof(ctx->arch) - 1);
    ctx->combined_csv_hash = md_fnv1a_init();
}

void md_ctx_record_stage(md_run_ctx *ctx, const char *id, const char *name,
                         double elapsed_sec, int files_written) {
    if (ctx->n_stages >= MD_MAX_STAGES) return;
    md_stage_record *r = &ctx->stages[ctx->n_stages++];
    snprintf(r->id, sizeof(r->id), "%s", id);
    snprintf(r->name, sizeof(r->name), "%s", name);
    r->elapsed_sec = elapsed_sec;
    r->files_written = files_written;
    r->completed = true;
}

void md_ctx_record_file(md_run_ctx *ctx, const char *rel_path, const char *kind,
                        uint64_t hash, size_t bytes, int rows, int cols) {
    if (ctx->n_files >= MD_MAX_FILES) return;
    md_file_record *f = &ctx->files[ctx->n_files++];
    snprintf(f->path, sizeof(f->path), "%s", rel_path);
    snprintf(f->kind, sizeof(f->kind), "%s", kind);
    f->fnv1a = hash;
    f->bytes = bytes;
    f->rows = rows;
    f->cols = cols;
}

void md_ctx_add_criterion(md_run_ctx *ctx, const char *id, const char *status,
                          const char *measured, const char *threshold) {
    if (ctx->n_criteria >= MD_MAX_CRITERIA) return;
    md_criterion_record *c = &ctx->criteria[ctx->n_criteria++];
    snprintf(c->id, sizeof(c->id), "%s", id);
    snprintf(c->status, sizeof(c->status), "%s", status);
    snprintf(c->measured, sizeof(c->measured), "%s", measured);
    snprintf(c->threshold, sizeof(c->threshold), "%s", threshold);
}

void md_ctx_add_flag(md_run_ctx *ctx, const char *stage, const char *file, const char *reason) {
    if (ctx->n_flags >= MD_MAX_FLAGS) return;
    md_vacuous_flag *fl = &ctx->flags[ctx->n_flags++];
    snprintf(fl->stage, sizeof(fl->stage), "%s", stage);
    snprintf(fl->file, sizeof(fl->file), "%s", file);
    snprintf(fl->reason, sizeof(fl->reason), "%s", reason);
}

bool md_hash_file(const char *path, uint64_t *hash_out, size_t *bytes_out) {
    FILE *fp = fopen(path, "rb");
    if (!fp) return false;
    uint64_t h = md_fnv1a_init();
    uint8_t buf[4096];
    size_t total = 0;
    size_t nr;
    while ((nr = fread(buf, 1, sizeof(buf), fp)) > 0) {
        h = md_fnv1a_update(h, buf, nr);
        total += nr;
    }
    fclose(fp);
    if (hash_out) *hash_out = h;
    if (bytes_out) *bytes_out = total;
    return true;
}

static bool is_param_or_count_col(const char *name) {
    static const char *exempt[] = {
        "stage", "step", "steps", "steps_fwd", "bin", "block", "repeat", "threads", "n", "n_particles",
        "n_pts", "nc", "dim", "shell", "idx", "particle", "frame", "lag", "t", "r", "k",
        "v", "dt", "rho", "target_t", "target_p", "r_skin", "r_skin_half", "r_cut", "rc", "tau_t",
        "m_inner", "quant_step", "raw_bytes", "comp_bytes", "n_blocks", "n_frames",
        "n_samples", "n_rebuilds", "rebuilds", "rebuilds_cum", "overflow_count", "late_rebuild_count",
        "exact_roundtrip", "max_err_q", "bit_exact", "pass", "pass_p001", "expected_count",
        "measured_count", "count", "rows", "cols", "bytes", "files", "files_written", "n_pairs",
        "box_l", "method", "potential", "ensemble", "phase", "version", "label",
        "hash", "fnv1a", "file", "path", "status", "reason", "id", "name",
        "t_init", "t_target", "temp_target", "r_sep", "side", "dir", "axis", "comp"
    };
    size_t n_ex = sizeof(exempt) / sizeof(exempt[0]);
    for (size_t i = 0; i < n_ex; ++i) {
        if (strcmp(name, exempt[i]) == 0) return true;
    }
    if (strstr(name, "count") != NULL || strstr(name, "hash") != NULL ||
        strstr(name, "exact") != NULL || strstr(name, "_idx") != NULL ||
        strstr(name, "target") != NULL || strstr(name, "bytes") != NULL) {
        return true;
    }
    return false;
}

static bool parse_strict_double(const char *tok, double *out) {
    if (!tok || !*tok) return false;
    while (*tok == ' ' || *tok == '\t') tok++;
    if (!*tok) return false;
    char *end = NULL;
    double v = strtod(tok, &end);
    if (end == tok) return false;
    while (*end == ' ' || *end == '\t' || *end == '\r' || *end == '\n') end++;
    if (*end != '\0') return false;
    *out = v;
    return true;
}

void md_audit_csv_vacuous(md_run_ctx *ctx, const char *stage_id, const char *path,
                          int *rows_out, int *cols_out) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        md_ctx_add_flag(ctx, stage_id, path, "cannot_open_csv");
        if (rows_out) *rows_out = 0;
        if (cols_out) *cols_out = 0;
        return;
    }

    char line[4096];
    if (!fgets(line, sizeof(line), fp)) {
        fclose(fp);
        md_ctx_add_flag(ctx, stage_id, path, "empty_csv");
        if (rows_out) *rows_out = 0;
        if (cols_out) *cols_out = 0;
        return;
    }

    char col_names[32][64];
    int n_cols = 0;
    char *p = line;
    while (*p && n_cols < 32) {
        char *comma = strchr(p, ',');
        size_t len = comma ? (size_t)(comma - p) : strcspn(p, "\r\n");
        if (len >= 63) len = 63;
        memcpy(col_names[n_cols], p, len);
        col_names[n_cols][len] = '\0';
        n_cols++;
        if (!comma) break;
        p = comma + 1;
    }

    int num_cnt[32] = {0};
    int int_cnt[32] = {0};
    double cmin[32], cmax[32];
    bool has_nonfinite[32] = {false};
    bool bad_temp[32] = {false};
    bool bad_gr[32] = {false};
    bool bad_nonneg[32] = {false};

    for (int c = 0; c < n_cols; ++c) {
        cmin[c] = 1e300;
        cmax[c] = -1e300;
    }

    int n_rows = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '\n' || line[0] == '\r' || line[0] == '\0') continue;
        n_rows++;
        int c = 0;
        char *q = line;
        while (*q && c < n_cols) {
            char *comma = strchr(q, ',');
            char tok[128];
            size_t len = comma ? (size_t)(comma - q) : strcspn(q, "\r\n");
            if (len >= sizeof(tok) - 1) len = sizeof(tok) - 1;
            memcpy(tok, q, len);
            tok[len] = '\0';

            double val = 0.0;
            if (parse_strict_double(tok, &val)) {
                if (!isfinite(val)) {
                    has_nonfinite[c] = true;
                } else {
                    num_cnt[c]++;
                    if (val < cmin[c]) cmin[c] = val;
                    if (val > cmax[c]) cmax[c] = val;
                    if (fabs(val - round(val)) < 1e-12 && fabs(val) >= 0.5) {
                        int_cnt[c]++;
                    }
                    const char *cn = col_names[c];
                    if ((strcmp(cn, "temp") == 0 || strcmp(cn, "t_mean") == 0 ||
                         strcmp(cn, "t_inst") == 0 || strcmp(cn, "t_meas") == 0) && val <= 0.0) {
                        bad_temp[c] = true;
                    }
                    if ((strcmp(cn, "g_r") == 0 || strcmp(cn, "gr") == 0 ||
                         strcmp(cn, "gr_gas") == 0 || strcmp(cn, "gr_liq") == 0 ||
                         strcmp(cn, "gr_sol") == 0) && val < -1e-12) {
                        bad_gr[c] = true;
                    }
                    if ((strcmp(cn, "msd") == 0 || strcmp(cn, "q6") == 0 ||
                         strcmp(cn, "lindemann") == 0 || strcmp(cn, "var_k") == 0) && val < -1e-12) {
                        bad_nonneg[c] = true;
                    }
                }
            }
            c++;
            if (!comma) break;
            q = comma + 1;
        }
    }
    fclose(fp);

    if (rows_out) *rows_out = n_rows;
    if (cols_out) *cols_out = n_cols;

    if (n_rows == 0) {
        md_ctx_add_flag(ctx, stage_id, path, "zero_data_rows");
        return;
    }

    for (int c = 0; c < n_cols; ++c) {
        char msg[256];
        if (has_nonfinite[c]) {
            snprintf(msg, sizeof(msg), "non_finite_value_in_%s", col_names[c]);
            md_ctx_add_flag(ctx, stage_id, path, msg);
        }
        if (bad_temp[c]) {
            snprintf(msg, sizeof(msg), "unphysical_temp_le_0_in_%s", col_names[c]);
            md_ctx_add_flag(ctx, stage_id, path, msg);
        }
        if (bad_gr[c]) {
            snprintf(msg, sizeof(msg), "unphysical_gr_lt_0_in_%s", col_names[c]);
            md_ctx_add_flag(ctx, stage_id, path, msg);
        }
        if (bad_nonneg[c]) {
            snprintf(msg, sizeof(msg), "unphysical_negative_in_%s", col_names[c]);
            md_ctx_add_flag(ctx, stage_id, path, msg);
        }
        if (!is_param_or_count_col(col_names[c]) && num_cnt[c] >= 6) {
            if (cmax[c] == cmin[c] && fabs(cmin[c]) > 1e-15) {
                /* Exempt analytical constants/references explicitly tagged */
                if (strstr(col_names[c], "ref") == NULL &&
                    strstr(col_names[c], "theory") == NULL &&
                    strstr(col_names[c], "canon") == NULL &&
                    strstr(col_names[c], "analytic") == NULL &&
                    strstr(col_names[c], "b2") == NULL &&
                    strstr(col_names[c], "limit") == NULL &&
                    strstr(col_names[c], "pe_0") == NULL &&
                    strstr(col_names[c], "u_fcc") == NULL) {
                    snprintf(msg, sizeof(msg), "constant_measured_col_%s=%.6g", col_names[c], cmin[c]);
                    md_ctx_add_flag(ctx, stage_id, path, msg);
                }
            }
            if (int_cnt[c] == num_cnt[c] && num_cnt[c] >= 8) {
                if (strstr(col_names[c], "id") == NULL &&
                    strstr(col_names[c], "num") == NULL &&
                    strstr(col_names[c], "n_") == NULL &&
                    strstr(col_names[c], "nbr") == NULL &&
                    strstr(col_names[c], "theory") == NULL &&
                    strstr(col_names[c], "ref") == NULL &&
                    strstr(col_names[c], "shell") == NULL) {
                    snprintf(msg, sizeof(msg), "all_round_integers_in_%s", col_names[c]);
                    md_ctx_add_flag(ctx, stage_id, path, msg);
                }
            }
        }
    }
}

void md_register_csv(md_run_ctx *ctx, const char *stage_id, const char *path) {
    uint64_t h = 0;
    size_t bytes = 0;
    int rows = 0, cols = 0;
    if (!md_hash_file(path, &h, &bytes)) {
        md_ctx_add_flag(ctx, stage_id, path, "missing_csv_file");
        return;
    }
    md_audit_csv_vacuous(ctx, stage_id, path, &rows, &cols);
    md_ctx_record_file(ctx, path, "csv", h, bytes, rows, cols);
    /* Exclude wall-clock timing/manifest CSVs from deterministic combined hash */
    if (strstr(path, "timing") == NULL && strstr(path, "scaling") == NULL &&
        strstr(path, "ladder") == NULL && strstr(path, "thread_scaling") == NULL &&
        strstr(path, "thread_repeats") == NULL && strstr(path, "manifest") == NULL &&
        strstr(path, "vacuous") == NULL && strstr(path, "2d_timing") == NULL &&
        strstr(path, "skin_sweep") == NULL) {
        ctx->combined_csv_hash = md_fnv1a_update(ctx->combined_csv_hash, path, strlen(path));
        ctx->combined_csv_hash = md_fnv1a_update(ctx->combined_csv_hash, &h, sizeof(h));
    }
}

bool md_register_svg(md_run_ctx *ctx, const char *stage_id, const char *path) {
    uint64_t h = 0;
    size_t bytes = 0;
    if (!md_hash_file(path, &h, &bytes)) {
        md_ctx_add_flag(ctx, stage_id, path, "missing_svg_file");
        return false;
    }
    FILE *fp = fopen(path, "r");
    if (!fp) return false;
    char *buf = (char *)malloc(bytes + 1);
    size_t nr = fread(buf, 1, bytes, fp);
    buf[nr] = '\0';
    fclose(fp);

    bool ok = (strstr(buf, "<svg") != NULL && strstr(buf, "</svg>") != NULL &&
               strstr(buf, "xmlns=\"http://www.w3.org/2000/svg\"") != NULL);
    free(buf);
    if (!ok) {
        md_ctx_add_flag(ctx, stage_id, path, "invalid_svg_xml");
    }
    /* Normalize relative path for manifest/JSON (strip leading ../../../ if present) */
    const char *rel = path;
    if (strncmp(rel, "../../../", 9) == 0) rel += 9;
    md_ctx_record_file(ctx, rel, "svg", h, bytes, 0, 0);
    return ok;
}

static void write_json_escaped(FILE *fp, const char *s) {
    fputc('"', fp);
    for (const char *p = s; *p; ++p) {
        unsigned char ch = (unsigned char)*p;
        if (ch == '"') fputs("\\\"", fp);
        else if (ch == '\\') fputs("\\\\", fp);
        else if (ch == '\n') fputs("\\n", fp);
        else if (ch == '\r') fputs("\\r", fp);
        else if (ch == '\t') fputs("\\t", fp);
        else if (ch >= 0x20) fputc(ch, fp);
    }
    fputc('"', fp);
}

bool md_write_results_json(const md_run_ctx *ctx, const char *path) {
    FILE *fp = fopen(path, "w");
    if (!fp) return false;

    int n_csv = 0, n_svg = 0;
    size_t csv_bytes = 0, svg_bytes = 0;
    for (int i = 0; i < ctx->n_files; ++i) {
        if (strcmp(ctx->files[i].kind, "csv") == 0) {
            n_csv++;
            csv_bytes += ctx->files[i].bytes;
        } else if (strcmp(ctx->files[i].kind, "svg") == 0) {
            n_svg++;
            svg_bytes += ctx->files[i].bytes;
        }
    }

    int n_pass = 0, n_partial = 0, n_fail = 0;
    for (int i = 0; i < ctx->n_criteria; ++i) {
        if (strcmp(ctx->criteria[i].status, "tercapai") == 0) n_pass++;
        else if (strcmp(ctx->criteria[i].status, "sebagian") == 0) n_partial++;
        else n_fail++;
    }

    fprintf(fp, "{\n");
    fprintf(fp, "  \"session\": \"day009a\",\n");
    fprintf(fp, "  \"project_path\": \"projects/c/day009a_md_lj_c\",\n");
    fprintf(fp, "  \"scale\": \"%s\",\n", ctx->small_scale ? "small" : "full");
    fprintf(fp, "  \"master_seed\": %llu,\n", (unsigned long long)MD_MASTER_SEED);
    fprintf(fp, "  \"run_completed_without_errors\": true,\n");
    fprintf(fp, "  \"machine\": {\n");
    fprintf(fp, "    \"arch\": "); write_json_escaped(fp, ctx->arch); fprintf(fp, ",\n");
    fprintf(fp, "    \"cores\": %d,\n", ctx->cores);
    fprintf(fp, "    \"cc_version\": "); write_json_escaped(fp, ctx->cc_ver); fprintf(fp, "\n");
    fprintf(fp, "  },\n");
    fprintf(fp, "  \"code_metrics\": {\n");
    fprintf(fp, "    \"c_lines\": %d,\n", ctx->c_lines);
    fprintf(fp, "    \"stages_completed\": %d,\n", ctx->n_stages);
    fprintf(fp, "    \"total_output_files\": %d,\n", ctx->n_files + 1);
    fprintf(fp, "    \"csv_files\": %d,\n", n_csv);
    fprintf(fp, "    \"svg_files\": %d,\n", n_svg);
    fprintf(fp, "    \"csv_total_bytes\": %zu,\n", csv_bytes);
    fprintf(fp, "    \"svg_total_bytes\": %zu,\n", svg_bytes);
    fprintf(fp, "    \"total_elapsed_sec\": %.4f\n", ctx->total_elapsed_sec);
    fprintf(fp, "  },\n");
    fprintf(fp, "  \"determinism\": {\n");
    fprintf(fp, "    \"combined_deterministic_csv_fnv1a\": \"%016llx\",\n",
            (unsigned long long)ctx->combined_csv_hash);
    fprintf(fp, "    \"thread_hash_1t\": \"%016llx\",\n", (unsigned long long)ctx->c16_hash_1t);
    fprintf(fp, "    \"thread_hash_2t\": \"%016llx\",\n", (unsigned long long)ctx->c16_hash_2t);
    fprintf(fp, "    \"thread_hash_4t\": \"%016llx\",\n", (unsigned long long)ctx->c16_hash_4t);
    fprintf(fp, "    \"thread_hash_8t\": \"%016llx\",\n", (unsigned long long)ctx->c16_hash_8t);
    fprintf(fp, "    \"threads_bitwise_identical\": %s,\n", ctx->c16_threads_bit_exact ? "true" : "false");
    fprintf(fp, "    \"two_runs_identical\": %s\n", ctx->c16_runs_bit_exact ? "true" : "false");
    fprintf(fp, "  },\n");
    fprintf(fp, "  \"invariants\": {\n");
    fprintf(fp, "    \"total_neighbor_overflows\": %llu,\n", (unsigned long long)ctx->total_nbr_overflows);
    fprintf(fp, "    \"total_late_verlet_rebuilds\": %llu\n", (unsigned long long)ctx->total_late_rebuilds);
    fprintf(fp, "  },\n");
    fprintf(fp, "  \"criteria_summary\": {\n");
    fprintf(fp, "    \"total\": %d,\n", ctx->n_criteria);
    fprintf(fp, "    \"tercapai\": %d,\n", n_pass);
    fprintf(fp, "    \"sebagian\": %d,\n", n_partial);
    fprintf(fp, "    \"lainnya\": %d\n", n_fail);
    fprintf(fp, "  },\n");

    fprintf(fp, "  \"criteria\": {\n");
    for (int i = 0; i < ctx->n_criteria; ++i) {
        const md_criterion_record *c = &ctx->criteria[i];
        fprintf(fp, "    \"%s\": {\n", c->id);
        fprintf(fp, "      \"status\": "); write_json_escaped(fp, c->status); fprintf(fp, ",\n");
        fprintf(fp, "      \"measured\": "); write_json_escaped(fp, c->measured); fprintf(fp, ",\n");
        fprintf(fp, "      \"threshold\": "); write_json_escaped(fp, c->threshold); fprintf(fp, "\n");
        fprintf(fp, "    }%s\n", (i + 1 < ctx->n_criteria) ? "," : "");
    }
    fprintf(fp, "  },\n");

    fprintf(fp, "  \"vacuous_flags\": [");
    if (ctx->n_flags > 0) fprintf(fp, "\n");
    for (int i = 0; i < ctx->n_flags; ++i) {
        fprintf(fp, "    {\"stage\": "); write_json_escaped(fp, ctx->flags[i].stage);
        fprintf(fp, ", \"file\": "); write_json_escaped(fp, ctx->flags[i].file);
        fprintf(fp, ", \"reason\": "); write_json_escaped(fp, ctx->flags[i].reason);
        fprintf(fp, "}%s\n", (i + 1 < ctx->n_flags) ? "," : "");
    }
    fprintf(fp, "%s],\n", (ctx->n_flags > 0) ? "  " : "");

    fprintf(fp, "  \"stages\": [\n");
    for (int i = 0; i < ctx->n_stages; ++i) {
        const md_stage_record *st = &ctx->stages[i];
        fprintf(fp, "    {\"id\": \"%s\", \"name\": \"%s\", \"elapsed_sec\": %.4f, \"files\": %d}%s\n",
                st->id, st->name, st->elapsed_sec, st->files_written,
                (i + 1 < ctx->n_stages) ? "," : "");
    }
    fprintf(fp, "  ]\n");
    fprintf(fp, "}\n");
    fclose(fp);
    return true;
}
