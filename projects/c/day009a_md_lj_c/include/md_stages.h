#ifndef MD_STAGES_H
#define MD_STAGES_H

#include "md_types.h"

/* Single built-in self-test suite (make test) */
int md_run_selftest(void);

/* Stage groups (60 stages S01..S60) */
void md_run_stages_a(md_run_ctx *ctx); /* S01..S10 */
void md_run_stages_b(md_run_ctx *ctx); /* S11..S18 */
void md_run_stages_c(md_run_ctx *ctx); /* S19..S26 */
void md_run_stages_d(md_run_ctx *ctx); /* S27..S34 */
void md_run_stages_e(md_run_ctx *ctx); /* S35..S44 */
void md_run_stages_f(md_run_ctx *ctx); /* S45..S54 */
void md_run_stages_g(md_run_ctx *ctx); /* S55..S58 */
void md_run_stages_h(md_run_ctx *ctx); /* S59..S60 + criteria + results.json */

/* Run full pipeline or single stage */
int md_run_all(md_run_ctx *ctx);
int md_run_thread_check(int n_threads);

#endif /* MD_STAGES_H */
