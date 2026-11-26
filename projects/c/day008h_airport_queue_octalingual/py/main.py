import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from analysis import write_s14_to_s45
from config import default_rp, load_config
from engine import sim_run
from network import Sim
from rng import fnv1a_bytes
from stats import write_s01_to_s13
from svg import write_s46_svgs


def manifest_scan(directory: str, prefix: str, out_f) -> None:
    for nm in sorted(os.listdir(directory)):
        if nm.startswith("."):
            continue
        p = os.path.join(directory, nm)
        if os.path.isdir(p):
            manifest_scan(p, f"{prefix}{nm}/", out_f)
        elif nm != "files_manifest.csv":
            with open(p, "rb") as bf:
                b = bf.read()
            nl = b.count(b"\n")
            out_f.write(f"{prefix}{nm},{nl},{len(b)},{fnv1a_bytes(b):016x}\n")


def main() -> None:
    out_dir = "out/py"
    os.makedirs(f"{out_dir}/figures", exist_ok=True)
    cfg = load_config("config.kv")
    if len(sys.argv) > 1 and sys.argv[1] == "s01":
        write_s01_to_s13(cfg, out_dir)
        return

    tot = 0.0
    for st in range(1, 14):
        print(f"[py] tahap {st}/48 (s{st:02d})")
    t0 = time.monotonic()
    write_s01_to_s13(cfg, out_dir)
    tot += time.monotonic() - t0

    main_days = 7
    print("[py] tahap 14/48 (s14)")
    t0 = time.monotonic()
    rpm = default_rp(cfg, main_days)
    rpm.audit_on = 1
    s_main = Sim(cfg, rpm)
    sim_run(s_main)
    tot += time.monotonic() - t0

    for st in range(15, 42):
        print(f"[py] tahap {st}/48 (s{st:02d})")
    t0 = time.monotonic()
    write_s14_to_s45(cfg, out_dir, s_main)
    tot += time.monotonic() - t0

    print("[py] tahap 42/48 (s42)")
    vals = []
    for _ in range(5):
        sb = Sim(cfg, default_rp(cfg, 3))
        tb = time.monotonic()
        sim_run(sb)
        dt = max(time.monotonic() - tb, 1e-9)
        vals.append(sb.sys_sum_l / dt)
    vals.sort()
    med = vals[2]
    dev = sorted(abs(v - med) for v in vals)
    with open(f"{out_dir}/s42_bench.csv", "w", encoding="utf-8") as f:
        f.write(f"slot_pax_med,slot_pax_per_sec_med,slot_pax_per_sec_mad\n1,{med:.2f},{dev[2]:.2f}\n")

    print("[py] tahap 43/48 (s43)\n[py] tahap 44/48 (s44)\n[py] tahap 45/48 (s45)\n[py] tahap 46/48 (s46)")
    t0 = time.monotonic()
    write_s46_svgs(s_main, out_dir)
    tot += time.monotonic() - t0

    print("[py] tahap 47/48 (s47)\n[py] tahap 48/48 (s48)")
    py_ver = f"python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    with open(f"{out_dir}/s48_summary.json", "w", encoding="utf-8") as f:
        f.write(f'{{"lang":"py","days":{main_days},"sw_days":1,"toolchain":"{py_ver}","total_sec":{tot:.3f}}}\n')
    with open(f"{out_dir}/files_manifest.csv", "w", encoding="utf-8") as mf:
        mf.write("file,lines,bytes,fnv1a\n")
        manifest_scan(out_dir, "", mf)


if __name__ == "__main__":
    main()
