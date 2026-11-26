import pathlib, sys, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import analysis, engine, entities, rng, svg

out_dir = "out/py"
pathlib.Path(f"{out_dir}/figures").mkdir(parents=True, exist_ok=True)
cfg = entities.load_config("config.kv")
cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
if cmd == "s01":
    analysis.write_s01_to_s09(cfg, out_dir)
    sys.exit(0)

tot = 0.0
for st in range(1, 10): print(f"[py] tahap {st}/40 (s{st:02d})")
t0 = time.perf_counter(); analysis.write_s01_to_s09(cfg, out_dir); tot += time.perf_counter() - t0
print("[py] tahap 10/40 (s10)")
t0 = time.perf_counter(); sim = entities.Sim(cfg, True, out_dir)
engine.sim_run_days(sim, cfg["num_days"], 8000, 0, -1, 1); sim.close_logs(); tot += time.perf_counter() - t0
for st in range(11, 39): print(f"[py] tahap {st}/40 (s{st:02d})")
t0 = time.perf_counter(); analysis.write_s11_to_s38(sim, cfg, out_dir); tot += time.perf_counter() - t0
print("[py] tahap 39/40 (s39)")
t0 = time.perf_counter(); svg.write_s39_svgs(sim, out_dir); tot += time.perf_counter() - t0
print("[py] tahap 40/40 (s40)")
(pathlib.Path(out_dir) / "s40_summary.json").write_text(
    f'{{"lang":"py","days":5,"ev_per_day":8000,"ev_total":{cfg["num_days"]*8000},"toolchain":"python 3.11.18","total_sec":{tot:.3f}}}\n'
)
mf = ["file,lines,bytes,fnv1a"]
for p in sorted(pathlib.Path(out_dir).rglob("*")):
    if p.is_dir() or p.name == "files_manifest.csv": continue
    rel = p.relative_to(out_dir).as_posix(); b = p.read_bytes()
    nl = b.count(b"\n")
    mf.append(f"{rel},{nl},{len(b)},{rng.fnv1a_bytes(b):016x}")
(pathlib.Path(out_dir) / "files_manifest.csv").write_text("\n".join(mf) + "\n")
