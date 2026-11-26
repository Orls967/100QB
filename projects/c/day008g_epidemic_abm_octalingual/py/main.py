import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from analysis import write_s01_to_s09, write_s10_to_s37
from engine import sim_run
from population import RunParams, Sim, load_config
from rng import fnv1a_bytes
from svg import write_s38_svgs


def manifest_scan(dir_path: str, prefix: str, out):
    for nm in sorted(os.listdir(dir_path)):
        if nm.startswith("."):
            continue
        p = os.path.join(dir_path, nm)
        if os.path.isdir(p):
            manifest_scan(p, f"{prefix}{nm}/", out)
        elif nm != "files_manifest.csv":
            with open(p, "rb") as f:
                b = f.read()
            lc = b.count(b"\n")
            out.write(f"{prefix}{nm},{lc},{len(b)},{fnv1a_bytes(b):016x}\n")


def main():
    out_dir = "out/py"
    os.makedirs(f"{out_dir}/figures", exist_ok=True)
    cfg = load_config("config.kv")
    if len(sys.argv) > 1 and sys.argv[1] == "s01":
        write_s01_to_s09(cfg, out_dir)
        return

    tot = 0.0
    for st in range(1, 10):
        print(f"[py] tahap {st}/39 (s{st:02d})")
    t0 = time.perf_counter()
    write_s01_to_s09(cfg, out_dir)
    tot += time.perf_counter() - t0

    n_main, beta_cal = 3000, 400
    print("[py] tahap 10/39 (s10)")
    t0 = time.perf_counter()
    s_no = Sim(cfg, RunParams(n_main, cfg["main_days"], beta_cal, interv_on=0), True)
    sim_run(s_no)
    tot += time.perf_counter() - t0

    print("[py] tahap 11/39 (s11)")
    t0 = time.perf_counter()
    s_main = Sim(cfg, RunParams(n_main, cfg["main_days"], beta_cal, interv_on=1), True)
    sim_run(s_main)
    tot += time.perf_counter() - t0

    for st in range(12, 35):
        print(f"[py] tahap {st}/39 (s{st:02d})")
    t0 = time.perf_counter()
    write_s10_to_s37(s_no, s_main, cfg, out_dir, beta_cal)
    tot += time.perf_counter() - t0

    print("[py] tahap 35/39 (s35)")
    vals = []
    for _ in range(5):
        sb = Sim(cfg, RunParams(2000, 60, beta_cal, interv_on=1), False)
        tb = time.perf_counter()
        sim_run(sb)
        dt = max(time.perf_counter() - tb, 1e-9)
        vals.append(120000.0 / dt)
    vals.sort()
    med = vals[2]
    dev = sorted(abs(v - med) for v in vals)
    with open(f"{out_dir}/s35_bench.csv", "w") as f:
        f.write(f"person_days,pd_per_sec_med,pd_per_sec_mad\n120000,{med:.2f},{dev[2]:.2f}\n")

    print("[py] tahap 36/39 (s36)\n[py] tahap 37/39 (s37)\n[py] tahap 38/39 (s38)")
    t0 = time.perf_counter()
    write_s38_svgs(s_no, s_main, out_dir)
    tot += time.perf_counter() - t0

    print("[py] tahap 39/39 (s39)")
    with open(f"{out_dir}/s39_summary.json", "w") as f:
        f.write(f'{{"lang":"py","N":{n_main},"days":{cfg["main_days"]},"N_sw":500,"toolchain":"Python 3.11.18","total_sec":{tot:.3f}}}\n')
    with open(f"{out_dir}/files_manifest.csv", "w") as mf:
        mf.write("file,lines,bytes,fnv1a\n")
        manifest_scan(out_dir, "", mf)


if __name__ == "__main__":
    main()
