#!/usr/bin/env python3
# NOTE: 本脚本源自项目内部流水线, 含两台计算机的绝对路径。
# 复现时请全局替换:
#   /data/wcf/MD_trend_analysis  -> 本地仓库根 (analysis 数据目录的父目录)
#   /data/md_extend              -> 远程机工作根
#   /home/rcdb/anaconda3/envs/AmberTools23 -> 远程 AMBERHOME (AmberTools23)
#   /data/Tools/Amber22          -> 本地 AMBERHOME (Amber22)
#   rcdb@100.76.253.60 / rcdb@100.97.18.74:10021 -> 计算机地址(需自备 ssh key)
"""统一参考点重算 RMSD 并拼接 原始100ns + 延伸段 (v2)。

参考点 = 各体系原始生产轨迹第 1 帧 (AF3 起始姿态, dry)。
原始轨迹树 (由 tasks_manifest prmtop 路径溯源):
  SH3: /data/wcf/AF3_prediction/IgGM_2d4d2_sh3_op_260126_part3_100ns_amber_try5/
  HCG: /data/wcf/protein_filter_lib/examples/hcg_r2/
本地延伸: /data/wcf/MD_trend_analysis/md_extend/runs_nvt/{SYSTEM}/{name}
远程延伸: rcdb@100.76.253.60:10021 /data/md_extend/runs/... (cpptraj@AmberTools23)

输出: extend_analysis/{SYSTEM}/{name}/unified_rmsd.dat
  用法: 16_unified_rmsd.py [--local] [--remote] [--plot]
"""
import json, subprocess, sys, tempfile
from pathlib import Path
import numpy as np

ROOT = Path('/data/wcf/MD_trend_analysis')
EA = ROOT / 'extend_analysis'
CPPTRAJ = '/data/Tools/Amber22/bin/cpptraj'
STRIP = ':WAT,K+,Cl-,Na+'
# 核心区掩码: 排除固有无序段与纯化标签后, 用于收敛判据 (见 handoff 十七节)
#   SH3: 去掉抗原脯氨酸富集 C 端 184-211 (RMSF 3.4-18.3 Å)
#   HCG: 去掉 VNAR 的 GGGS-10xHis-FLAG 标签 393-414 (RMSF 4-15 Å)
CORE_MASKS = {'SH3_VNAR': ':1-183', 'HCG_VNAR': ':1-392'}
SSH = ['ssh', '-p', '10021', '-i', str(ROOT/'md_extend'/'deploy_key'),
       '-o', 'StrictHostKeyChecking=no', 'rcdb@100.76.253.60']
SCPH = ['-P', '10021', '-i', str(ROOT/'md_extend'/'deploy_key'),
        '-o', 'StrictHostKeyChecking=no']
RCPPTRAJ = '/home/rcdb/anaconda3/envs/AmberTools23/bin/cpptraj'

def sh(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)

def run_cpptraj(script, cwd, exe=CPPTRAJ):
    with tempfile.NamedTemporaryFile('w', suffix='.in', delete=False) as f:
        f.write(script); fin = f.name
    r = sh([exe, '-i', fin], cwd=cwd)
    Path(fin).unlink()
    return r.returncode == 0, r.stdout + r.stderr

def load_systems():
    recs = []
    for j in ['trend_summary.json', 'trend_summary_remote.json']:
        p = EA / j
        if p.exists():
            for r in json.load(open(p)):
                r['host'] = 'remote' if 'remote' in j else 'local'
                recs.append(r)
    best = {}
    for r in recs:
        k = f"{r['system']}/{r['name']}"
        if k not in best or r['nseg'] > best[k]['nseg']: best[k] = r
    return best

def manifest_run_dirs():
    m = ROOT/'md_extend'/'tasks_manifest.json'
    out = {}
    if m.exists():
        d = json.load(open(m)); ts = d if isinstance(d, list) else d.get('tasks', [])
        for t in ts:
            p = t.get('prmtop')
            if p: out[t['name']] = str(Path(p).parent)
    return out

def find_orig_run(name, mf):
    if name in mf and Path(mf[name]).is_dir():
        return Path(mf[name])
    trees = [Path('/data/wcf/AF3_prediction/IgGM_2d4d2_sh3_op_260126_part3_100ns_amber_try5'),
             Path('/data/wcf/protein_filter_lib/examples/hcg_r2'),
             Path('/data/wcf/protein_filter_lib/examples/affinity_maturation_example/part3_amber_out')]
    for t in trees:
        hits = sorted([d for d in t.rglob(name) if d.is_dir()])
        for d in hits:
            if list(d.glob('md_total.nc')) or list(d.glob('md_1.nc')):
                return d
    return None

def pick_nc(d):
    cands = [p for p in sorted(d.glob('*.nc'), key=lambda p: p.stat().st_mtime, reverse=True)
             if not any(s in p.name.lower() for s in ('nowat', 'dry', 'strip', 'complex_only', 'dummy'))]
    # 优先 md_total.nc（完整生产轨迹）；原始 100ns 若曾中断续跑会拆成
    # md_1.nc(短段)+md_2.nc(续跑)，md_1.nc 只含前 ~1.6ns，切不可当作完整 100ns。
    cands.sort(key=lambda q: 0 if q.name == 'md_total.nc' else (1 if q.name == 'md_1.nc' else 2))
    return cands

def make_ref(orig, out_dir):
    """原始生产第1帧 -> dry pdb 公共参考 (自动跳过 dry 轨迹)。"""
    prmtop = orig/'system.prmtop'
    if not prmtop.exists():
        prmtop = next(iter(orig.glob('*.prmtop')))
    for nc in pick_nc(orig):
        ok, _ = run_cpptraj(f"""parm {prmtop}
trajin {nc} 1 1
strip {STRIP}
autoimage
trajout {out_dir}/init_dry.pdb pdb
""", out_dir)
        if ok and (out_dir/'init_dry.pdb').exists():
            return prmtop
        (out_dir/'init_dry.pdb').unlink(missing_ok=True)
    return None

def rms_local(parm, traj, out_dir, outfile, ref_dir=None, last=None, core_mask=None):
    """单次遍历同时算全复合物(outfile)与核心区(outfile 的 _core 版) RMSD。

    core_mask 用于排除固有无序段/纯化标签(见 CORE_MASKS), 供收敛判据使用。"""
    ref_dir = ref_dir or out_dir
    trajin = f"trajin {traj}" + (f" 1 {last}" if last else "")
    core_line = ""
    if core_mask:
        core_out = str(outfile).replace('.dat', '_core.dat')
        core_line = f"rms coreref ref initref {core_mask}@CA,C,N out {core_out}\n"
    ok, _ = run_cpptraj(f"""parm {parm}
parm {ref_dir}/init_dry.pdb [refp]
reference {ref_dir}/init_dry.pdb parm [refp] name initref
{trajin}
strip {STRIP}
autoimage
rms bb ref initref @CA,C,N out {outfile}
{core_line}run
quit
""", out_dir)
    return ok and outfile.exists()

def dat_xy(dat, dt=0.01, t0=0.0):
    a = np.loadtxt(dat, comments=('#', '@'))
    return a[:, 0]*dt + t0, a[:, 1]

def probe_orig_dt(orig: Path):
    """原始生产 ns/帧 = ntwx*dt(0.002ps)。从 md*.in 探测 ntwx, 失败回退 5000。"""
    import re
    for f in sorted(orig.glob('md*.in')):
        m = re.search(r'ntwx\s*=\s*(\d+)', f.read_text(errors='ignore'))
        if m:
            return int(m.group(1)) * 2e-6
    return 0.01  # 兜底 (5000*0.002)

NC_DUMP = '/data/Tools/Amber22/bin/ncdump'

# 跨服务器迁移体系: 本地先跑了 NPT 段, 之后在远程用 NVT 续跑。
# 其本地段是同一轨迹的**前缀**(非重复), 必须计入; 其余体系的本地目录若有残留
# 段多为失败尝试(如盒子爆炸), 故用显式白名单而非自动探测。
LOCAL_NPT_PREFIX = {'S85T_model', 'S86R_Y92D_model', 'S85I_G91D_H95Y_model'}

def nc_span_ns(nc: Path):
    """从 NetCDF time 变量读真实时长(ns) = t[-1] - t[0] (权威口径)。
    生产 TIME 时钟含前序平衡偏移时(如 .out 最大 140.8ns 而实际只有 100ns),
    用 .out 的 TIME 最大值会把 span 高估, 使 n_last 被按比例错误截断。"""
    try:
        from scipy.io import netcdf_file
        with netcdf_file(str(nc), 'r', mmap=True) as f:
            if 'time' not in f.variables:
                return None
            t = np.asarray(f.variables['time'][:], dtype=float)
    except Exception:
        return None
    if t.size < 2:
        return None
    return float(t[-1] - t[0]) / 1000.0

def nc_box_stable(nc: Path, tol: float = 1.5):
    """盒子体积稳定性检查(挡掉"盒子爆炸"的失败段)。None = 无法判断。"""
    try:
        from scipy.io import netcdf_file
        with netcdf_file(str(nc), 'r', mmap=True) as f:
            if 'cell_lengths' not in f.variables:
                return None
            L = np.asarray(f.variables['cell_lengths'][:], dtype=float)
    except Exception:
        return None
    if L.ndim != 2 or L.shape[0] < 2:
        return None
    vs = L.prod(axis=1); vs = vs[vs > 0]
    return bool(vs.max() / vs.min() <= tol) if vs.size > 1 else None

def seg_key(p: Path):
    """md_seg10.nc 需排在 md_seg9.nc 之后 (自然序)"""
    import re
    m = re.search(r'(\d+)', p.stem)
    return (int(m.group(1)) if m else 0)

def probe_span_frames(orig: Path, nc: Path):
    """返回 (原始段时长ns, 帧数)。时长首选 NetCDF time 变量(权威), 回退 .out TIME(PS)。"""
    import re, subprocess
    frames = 0
    if NC_DUMP and Path(NC_DUMP).exists():
        r = subprocess.run([NC_DUMP, '-h', str(nc)], capture_output=True, text=True)
        mfr = re.search(r'frame = UNLIMITED ;\s*// \((\d+) currently\)', r.stdout)
        if mfr: frames = int(mfr.group(1))
    span = nc_span_ns(nc)
    if span is None:                    # 回退: .out 最大值 (可能含平衡时钟偏移)
        tmax = 0.0
        for f in orig.glob('md*.out'):
            for m in re.finditer(r'TIME.PS. =\s*([0-9.]+)', f.read_text(errors='ignore')):
                tmax = max(tmax, float(m.group(1)))
        span = tmax / 1000.0
    if frames and span > 0.1:
        return span, frames
    return 100.0, frames

def write_unified(out_dir, ts, vs, notes, fname='unified_rmsd.dat', title='@CA,C,N'):
    with open(out_dir/fname, 'w') as f:
        f.write(f"# time_ns rmsd_A (公共参考=原始生产第1帧, {title})\n")
        for n in notes: f.write(n+'\n')
        for t, v in zip(ts, vs):
            for ti, vi in zip(t, v): f.write(f"{ti:.2f} {vi:.3f}\n")

# ---------- 本地 ----------
def process_local(name, rec, mf):
    if rec.get('host') == 'remote':
        return f"{name}: 远程体系, 交给--remote"
    sys_name = rec['system']
    cmask = CORE_MASKS.get(sys_name)
    out_dir = EA/sys_name/name; out_dir.mkdir(parents=True, exist_ok=True)
    run_dir = rec.get('run_dir')
    run_dir = Path(run_dir) if run_dir and Path(run_dir).is_dir() \
              else ROOT/'md_extend'/'runs_nvt'/sys_name/name
    if not list(run_dir.glob('md_seg*.nc')):
        return f"{name}: 本地无延伸轨迹, 跳过(远程处理)"
    orig = find_orig_run(name, mf)
    ts, vs, cs, cvs, notes = [], [], [], [], []
    t_old = 0.0
    if orig:
        prmtop = make_ref(orig, out_dir)
        if prmtop:
            for nc in pick_nc(orig):
                span, frames = probe_span_frames(orig, nc)
                if not frames:
                    (out_dir/'_old.dat').unlink(missing_ok=True); continue
                n_last = frames if span <= 110 else max(int(round(frames * 100.0 / span)), 1)
                if not rms_local(prmtop, nc, out_dir, out_dir/'_old.dat',
                                 last=n_last, core_mask=cmask):
                    (out_dir/'_old.dat').unlink(missing_ok=True); continue
                dt = span / max(frames, 1)      # 每帧 ns (与读取的 n_last 无关)
                t, v = dat_xy(out_dir/'_old.dat', dt=dt)
                ts.append(t); vs.append(v)
                if cmask and (out_dir/'_old_core.dat').exists():
                    tc, vc = dat_xy(out_dir/'_old_core.dat', dt=dt)
                    cs.append(tc); cvs.append(vc)
                t_old = t[-1]; break
        else:
            notes.append("# NO_ORIGINAL_REF")
    else:
        notes.append("# NO_ORIGINAL_RUN")
    for seg in sorted(run_dir.glob('md_seg*.nc'), key=seg_key):
        if rms_local(run_dir/'system.prmtop', seg, out_dir, out_dir/'_ext.dat',
                     core_mask=cmask):
            t, v = dat_xy(out_dir/'_ext.dat', dt=0.01, t0=t_old)  # ntwx=5000, dt=2fs
            ts.append(t); vs.append(v)
            if cmask and (out_dir/'_ext_core.dat').exists():
                tc, vc = dat_xy(out_dir/'_ext_core.dat', dt=0.01, t0=t_old)
                cs.append(tc); cvs.append(vc)
            t_old = t[-1]  # 段间时间轴累进, 避免多段重叠
        for f in ('_ext.dat', '_ext_core.dat'):
            (out_dir/f).unlink(missing_ok=True)
    if not ts: return f"{name}: FAIL"
    if t_old < 90: notes.append(f"# REBUILD_AT {t_old:.1f}ns (原始生产中断于 {t_old:.1f}ns)")
    elif abs(t_old - 100) > 2: notes.append(f"# ORIG_LEN {t_old:.1f}ns")
    write_unified(out_dir, ts, vs, notes)
    if cs and len(cs) == len(ts):
        write_unified(out_dir, cs, cvs, notes, 'unified_rmsd_core.dat',
                      f'core {cmask}@CA,C,N')
    return f"{name}: OK t_old={t_old:.1f}ns segs={len(ts)-bool(ts and t_old>0)}"

# ---------- 远程 ----------
REMOTE_DRIVER = r'''
import json, re, subprocess, sys
from pathlib import Path
CPP = "__CPP__"
STRIP = ":WAT,K+,Cl-,Na+"
COREMASK = "__COREMASK__"
def cp(inp, cwd):
    r = subprocess.run([CPP, "-i", inp], capture_output=True, text=True, cwd=cwd)
    return r.returncode == 0, r.stdout[-800:]
def sk(p):
    m = re.search(r"(\d+)", p.stem); return int(m.group(1)) if m else 0
plan = json.load(open("__PLAN__"))
for name, info in plan.items():
    run_dir = Path(info["run_dir"]); ref = Path(info["ref_dir"])
    if not ref.joinpath("init_dry.pdb").exists():
        print(name, "NO_REF"); continue
    outs = []
    for seg in sorted(run_dir.glob("md_seg*.nc"), key=sk):
        # 关键守卫: 只处理**已终结**的段。pmemd 运行时 .nc 即存在(且持续长大),
        # 若此时算, 会把半成品当完整段, 且下面的"存在即跳过"幂等会**永久复用**它。
        # 权威判据 = 对应 .out 含 STOP/Final (见 handoff 十八节教训);
        # .nc/.rst 的存在性或大小都不可用(.rst 因 ntwr=5000 每 10ps 就写)。
        segout = seg.with_suffix(".out")
        if not segout.exists() or not re.search(
                r"STOP|Final", segout.read_text(errors="ignore")[-20000:]):
            for stale in (ref/f"_ext_{seg.stem}.dat", ref/f"_ext_{seg.stem}_core.dat"):
                stale.unlink(missing_ok=True)      # 清掉可能残留的半成品缓存
            print(name, seg.name, "SKIP_UNFINISHED"); continue
        whole = ref/f"_ext_{seg.stem}.dat"
        core = ref/f"_ext_{seg.stem}_core.dat"
        # 幂等: 两个文件都在 **且比该段 .out 新** 才跳过。
        # 若缓存在段仍在运行时生成(旧行为的遗留), 它比 .out 旧 -> 强制重算。
        if (whole.exists() and (not COREMASK or core.exists())
                and whole.stat().st_mtime >= segout.stat().st_mtime):
            outs.append(seg.stem); continue
        coreline = (f"rms coreref ref initref {COREMASK}@CA,C,N out {core}\n"
                    if COREMASK else "")
        script = f"""parm {run_dir}/system.prmtop
parm {ref}/init_dry.pdb [refp]
reference {ref}/init_dry.pdb parm [refp] name initref
trajin {seg}
strip {STRIP}
autoimage
rms bb ref initref @CA,C,N out {whole}
{coreline}"""
        p = ref/f"_in_{seg.stem}.in"; p.write_text(script)
        ok, msg = cp(str(p), str(ref)); p.unlink(missing_ok=True)
        if not ok: print(name, seg.name, "FAIL", msg); continue
        outs.append(seg.stem)
    # 回传本次实际纳入的段清单: 本地据此清除残留的半成品 _ext 缓存
    # (scp 只覆盖不删除, 若不清则会继续合并上一次留存的半成品)
    ref.joinpath("_included.json").write_text(json.dumps(outs))
    print(name, "OK", outs)
'''

def process_remote(names, mf):
    best = load_systems()
    best = {k: r for k, r in best.items() if r.get('host') == 'remote'}
    if names:
        best = {k: r for k, r in best.items() if r['name'] in names}
    results = {}
    for k, rec in sorted(best.items()):
        name, sys_name = rec['name'], rec['system']
        out_dir = EA/sys_name/name; out_dir.mkdir(parents=True, exist_ok=True)
        if not (out_dir/'init_dry.pdb').exists():
            orig = find_orig_run(name, mf)
            if not orig or not make_ref(orig, out_dir):
                results[name] = 'NO_ORIG'; continue
        rdir = f"/data/md_extend/refs/{name}"
        subprocess.run(['ssh', *SSH[1:], 'mkdir', '-p', rdir], capture_output=True)
        subprocess.run(['scp', *SCPH, str(out_dir/'init_dry.pdb'), f'rcdb@100.76.253.60:{rdir}/'],
                       capture_output=True)
        # 远程 run 目录: runs 或 runs_nvt
        r = subprocess.run(['ssh', *SSH[1:],
                            f'for T in runs runs_nvt; do d=/data/md_extend/$T/{sys_name}/{name}; '
                            f'[ -d "$d" ] && echo $d; done'],
                           capture_output=True, text=True)
        rrdir = r.stdout.strip().splitlines()
        if not rrdir:
            results[name] = 'NO_REMOTE_RUN'; continue
        rrun = rrdir[-1]
        plan = {name: {'run_dir': rrun, 'ref_dir': rdir}}
        pf = out_dir/'_remote_plan.json'; pf.write_text(json.dumps(plan))
        subprocess.run(['scp', *SCPH, str(pf),
                        f'rcdb@100.76.253.60:{rdir}/plan.json'], capture_output=True)
        script = REMOTE_DRIVER.replace('__CPP__', RCPPTRAJ).replace(
            '__PLAN__', f'{rdir}/plan.json').replace(
            '__COREMASK__', (CORE_MASKS.get(sys_name) or ''))
        sf = out_dir/'_remote_driver.py'; sf.write_text(script)
        subprocess.run(['scp', *SCPH, str(sf), f'rcdb@100.76.253.60:{rdir}/driver.py'],
                       capture_output=True)
        r = subprocess.run(['ssh', *SSH[1:],
                            'python3', f'{rdir}/driver.py'], capture_output=True, text=True)
        print(r.stdout.strip()[-300:])
        # 拉回合并
        subprocess.run(['scp', *SCPH, f'rcdb@100.76.253.60:{rdir}/_ext_md_seg*.dat', str(out_dir) + '/'],
                       capture_output=True)
        subprocess.run(['scp', *SCPH, f'rcdb@100.76.253.60:{rdir}/_included.json', str(out_dir) + '/'],
                       capture_output=True)
        # 按远程回传的"实际纳入段"清理本地残留: scp 只覆盖不删除,
        # 上一次留存但本轮被跳过(段未终结)的 _ext 缓存必须删掉, 否则仍会被合并进来。
        # 注意命名: driver 回传 'md_seg5', 本地文件名为 '_ext_md_seg5(.dat|_core.dat)'。
        inc = out_dir/'_included.json'
        if inc.exists():
            keep = set(json.load(open(inc)))
            for f in list(out_dir.glob('_ext_md_seg*.dat')):
                stem = f.stem.replace('_core', '').replace('_ext_', '')
                if stem not in keep:
                    f.unlink(missing_ok=True)
                    print(f"  清除残留半成品缓存: {f.name}")
        cmask = CORE_MASKS.get(sys_name)
        ts, vs, cs, cvs, notes = [], [], [], [], []
        t_old = 0.0
        orig = find_orig_run(name, mf)
        if orig:
            if not (out_dir/'init_dry.pdb').exists():
                make_ref(orig, out_dir)
            if not (out_dir/'_old.dat').exists() and (out_dir/'init_dry.pdb').exists():
                prmtop = next(iter(orig.glob('system.prmtop'))) if (orig/'system.prmtop').exists() \
                         else next(iter(orig.glob('*.prmtop')))
                for nc in pick_nc(orig):
                    span, frames = probe_span_frames(orig, nc)
                    if not frames:
                        (out_dir/'_old.dat').unlink(missing_ok=True); continue
                    n_last = frames if span <= 110 else max(int(round(frames * 100.0 / span)), 1)
                    if rms_local(prmtop, nc, out_dir, out_dir/'_old.dat', last=n_last,
                                 core_mask=cmask):
                        break
                    (out_dir/'_old.dat').unlink(missing_ok=True)
            if (out_dir/'_old.dat').exists():
                dt = 100.0 / max(len(np.loadtxt(out_dir/"_old.dat", comments=("#", "@"))), 1)
                t, v = dat_xy(out_dir/'_old.dat', dt=dt)
                ts.append(t); vs.append(v); t_old = t[-1]
                if cmask and (out_dir/'_old_core.dat').exists():
                    tc, vc = dat_xy(out_dir/'_old_core.dat', dt=dt)
                    cs.append(tc); cvs.append(vc)
            else:
                notes.append("# ORIGINAL_SEG_FAILED")
        else:
            notes.append("# NO_ORIGINAL_RUN")
        # 本地 NPT 前缀段 (仅白名单体系): 同一轨迹的前缀, 接在原始段之后、远程 NVT 段之前
        if name in LOCAL_NPT_PREFIX:
            loc_dir = ROOT/'md_extend'/'runs_nvt'/sys_name/name
            n_loc = 0
            for seg in sorted(loc_dir.glob('md_seg*.nc'), key=seg_key):
                if nc_box_stable(seg) is False:            # 盒子爆炸的失败段 -> 跳过
                    notes.append(f"# SKIP_UNSTABLE_BOX {seg.stem}"); continue
                tmp = out_dir/f'_loc_{seg.stem}.dat'
                if rms_local(loc_dir/'system.prmtop', seg, out_dir, tmp, core_mask=cmask):
                    t, v = dat_xy(tmp, dt=0.01, t0=t_old)  # ntwx=5000 -> 0.01ns/帧
                    ts.append(t); vs.append(v)
                    tmpc = out_dir/f'_loc_{seg.stem}_core.dat'
                    if cmask and tmpc.exists():
                        tc, vc = dat_xy(tmpc, dt=0.01, t0=t_old)
                        cs.append(tc); cvs.append(vc)
                    t_old = t[-1]; n_loc += 1
                for f in (f'_loc_{seg.stem}.dat', f'_loc_{seg.stem}_core.dat'):
                    (out_dir/f).unlink(missing_ok=True)
            if n_loc:
                notes.append(f"# LOCAL_NPT_PREFIX {n_loc}segs")
        # 远程延伸段: 需排除 *_core.dat (核心区版本, 单独收集)
        exts = [e for e in sorted(out_dir.glob('_ext_md_seg*.dat'), key=seg_key)
                if not e.stem.endswith('_core')]
        for e in exts:
            t, v = dat_xy(e, t0=t_old); ts.append(t); vs.append(v)
            ec = e.with_name(e.stem + '_core.dat')
            if cmask and ec.exists():
                tc, vc = dat_xy(ec, t0=t_old); cs.append(tc); cvs.append(vc)
            t_old = t[-1]  # 段间时间轴累进, 避免多段重叠
        if not ts:
            results[name] = 'FAIL'; continue
        write_unified(out_dir, ts, vs, notes)
        if cs and len(cs) == len(ts):
            write_unified(out_dir, cs, cvs, notes, 'unified_rmsd_core.dat',
                          f'core {cmask}@CA,C,N')
        results[name] = f"OK {len(exts)}segs"
    print(json.dumps(results, ensure_ascii=False, indent=1))

if __name__ == '__main__':
    mf = manifest_run_dirs()
    if '--remote' in sys.argv:
        process_remote([a for a in sys.argv[2:] if not a.startswith('-')], mf)
    else:
        best = load_systems()
        import os
        only = os.environ.get('ONLY')
        if only:
            names = set(only.split(','))
            best = {k: r for k, r in best.items() if r['name'] in names}
        for k, r in sorted(best.items()):
            print(process_local(r['name'], r, mf))
