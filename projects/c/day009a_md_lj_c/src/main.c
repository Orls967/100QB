#include "md_io.h"
#include "md_stages.h"
#include "md_stats.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_usage(const char *prog) {
    printf("Usage:\n");
    printf("  %s selftest\n", prog);
    printf("  %s all [--scale full|small] [--lines N] [--cc STR] [--arch STR] [--cores N]\n", prog);
    printf("  %s threads <1|2|4|8>\n", prog);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "selftest") == 0) {
        return md_run_selftest();
    }

    if (strcmp(argv[1], "threads") == 0) {
        int nt = (argc >= 3) ? atoi(argv[2]) : 4;
        return md_run_thread_check(nt);
    }

    if (strcmp(argv[1], "all") == 0) {
        md_run_ctx ctx;
        md_ctx_init(&ctx);
        for (int i = 2; i < argc; ++i) {
            if (strcmp(argv[i], "--scale") == 0 && i + 1 < argc) {
                ctx.small_scale = (strcmp(argv[++i], "small") == 0);
            } else if (strcmp(argv[i], "--lines") == 0 && i + 1 < argc) {
                ctx.c_lines = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--cc") == 0 && i + 1 < argc) {
                snprintf(ctx.cc_ver, sizeof(ctx.cc_ver), "%s", argv[++i]);
            } else if (strcmp(argv[i], "--arch") == 0 && i + 1 < argc) {
                snprintf(ctx.arch, sizeof(ctx.arch), "%s", argv[++i]);
            } else if (strcmp(argv[i], "--cores") == 0 && i + 1 < argc) {
                ctx.cores = atoi(argv[++i]);
            }
        }
        return md_run_all(&ctx);
    }

    print_usage(argv[0]);
    return 1;
}
