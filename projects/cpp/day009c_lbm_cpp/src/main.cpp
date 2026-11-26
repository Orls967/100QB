#include "lbm/io.hpp"
#include "lbm/stages.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::printf("Usage: lbm <selftest|all|group|threads> [options]\n");
        return 1;
    }

    const std::string cmd = argv[1];
    if (cmd == "selftest") {
        return lbm::run_selftest() ? 0 : 1;
    }

    if (cmd == "threads") {
        const int nt = (argc >= 3) ? std::atoi(argv[2]) : 4;
        const std::uint64_t h = lbm::benchmark_thread_parity_hash(nt);
        std::printf("threads=%d fnv1a=%s\n", nt, lbm::hex64(h).c_str());
        return 0;
    }

    lbm::RunContext ctx;
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--scale" && i + 1 < argc) ctx.scale = argv[++i];
        else if (a == "--threads" && i + 1 < argc) ctx.threads = std::atoi(argv[++i]);
        else if (a == "--lines" && i + 1 < argc) ctx.cpp_lines = static_cast<std::size_t>(std::strtoull(argv[++i], nullptr, 10));
        else if (a == "--compiler" && i + 1 < argc) ctx.compiler_version = argv[++i];
        else if (a == "--arch" && i + 1 < argc) ctx.arch = argv[++i];
        else if (a == "--cores" && i + 1 < argc) ctx.cores = std::atoi(argv[++i]);
    }

    if (cmd == "group") {
        const std::string grp = (argc >= 3) ? argv[2] : "A";
        if (grp == "A") lbm::run_group_a(ctx);
        else if (grp == "B") lbm::run_group_b(ctx);
        else if (grp == "C") lbm::run_group_c(ctx);
        else if (grp == "D") lbm::run_group_d(ctx);
        else if (grp == "E") lbm::run_group_e(ctx);
        else if (grp == "F") lbm::run_group_f(ctx);
        else if (grp == "G") lbm::run_group_g(ctx);
        else if (grp == "H") lbm::run_group_h(ctx);
        return 0;
    }

    if (cmd == "all") {
        std::printf("[lbm] Running selftest first...\n");
        ctx.selftest_ok = lbm::run_selftest();
        if (!ctx.selftest_ok) return 1;

        std::printf("[lbm] Running Group A...\n");
        lbm::run_group_a(ctx);
        std::printf("[lbm] Running Group B...\n");
        lbm::run_group_b(ctx);
        std::printf("[lbm] Running Group C...\n");
        lbm::run_group_c(ctx);
        std::printf("[lbm] Running Group D...\n");
        lbm::run_group_d(ctx);
        std::printf("[lbm] Running Group E...\n");
        lbm::run_group_e(ctx);
        std::printf("[lbm] Running Group F...\n");
        lbm::run_group_f(ctx);
        std::printf("[lbm] Running Group G...\n");
        lbm::run_group_g(ctx);
        std::printf("[lbm] Running Group H...\n");
        lbm::run_group_h(ctx);

        std::printf("[lbm] Completed %zu stages in %.2fs | Files=%zu (+1 JSON) | Det-CSV-FNV1a=%s | Flags=%zu\n",
                    ctx.stages.size(), ctx.total_elapsed_sec, ctx.files.size(),
                    lbm::hex64(ctx.det_csv_hash).c_str(), ctx.vacuous_flags.size());
        return 0;
    }

    std::printf("Unknown command: %s\n", cmd.c_str());
    return 1;
}
