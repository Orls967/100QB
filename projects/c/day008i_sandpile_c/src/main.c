#include "sandpile.h"
#include <string.h>

int main(int argc, char **argv) {
    const char *cmd = (argc > 1) ? argv[1] : "selftest";
    Metrics m;
    memset(&m, 0, sizeof(m));
    if (strcmp(cmd, "selftest") == 0) return cmd_selftest(&m);
    if (strcmp(cmd, "abelian") == 0) return cmd_abelian(&m, "out");
    if (strcmp(cmd, "group") == 0) return cmd_group(&m, "out");
    fprintf(stderr, "unknown or pending subcommand: %s\n", cmd);
    return 1;
}
