#!/usr/bin/env python3
"""
Step 6: MD extension runner — manages GPU scheduling and runs pmemd.cuda.

Design:
  - 3 GPUs (physical IDs 2, 4, 7) run tasks concurrently.
  - Each task = N x 100ns segments (4 for high-priority, 2 for medium).
  - Segments run sequentially within a GPU slot.
  - When a structure finishes all segments, the next pending task is dispatched
    to the freed GPU.
  - Each segment: pmemd.cuda -O -i md_extend_100ns.in -o md_segK.out -p system.prmtop
                   -c md_segK_prev.rst -r md_segK.rst -x md_segK.nc
  - Segment 1 restarts from the original md_1.rst (end of 100ns production).
  - After all segments done, trajectories are concatenated with cpptraj.

Usage:
  python3 06_run_md_extend.py                 # run all tasks
  python3 06_run_md_extend.py --status         # show current progress
  python3 06_run_md_extend.py --gpu 2 --task 0 # run specific task on specific GPU
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

MANIFEST = Path(os.environ.get(
    "MD_EXTEND_MANIFEST",
    "/data/wcf/MD_trend_analysis/md_extend/tasks_manifest.json"))
RUN_ROOT = Path(os.environ.get(
    "MD_EXTEND_RUN_ROOT", "/data/wcf/MD_trend_analysis/md_extend/runs"))
LOG_ROOT = Path(os.environ.get(
    "MD_EXTEND_LOG_ROOT", "/data/wcf/MD_trend_analysis/md_extend/logs"))

AMBERHOME = os.environ.get("MD_EXTEND_AMBERHOME", "/data/Tools/Amber22")
PMEMD = os.environ.get("MD_EXTEND_PMEMD", f"{AMBERHOME}/bin/pmemd.cuda")
MD_INPUT = Path(os.environ.get(
    "MD_EXTEND_MDIN", "/data/wcf/MD_trend_analysis/md_extend/md_extend_100ns.in"))
EQUIL_INPUT = Path(os.environ.get(
    "MD_EXTEND_EQUILIN", "/data/wcf/MD_trend_analysis/md_extend/md_equil_1ns.in"))
CPPTRAJ = f"{AMBERHOME}/bin/cpptraj"

GPUS = [2, 4, 7]

# Box-guard thresholds (Angstrom): swollen-box NPT instability guard
BOX_WARN_A = 150.0   # warn but continue
BOX_STOP_A = 250.0   # mark task unstable, stop running further segments

# Physical GPU index -> UUID (required because GPU 1/3 are dead,
# causing CUDA to re-number devices; UUIDs are unambiguous)
GPU_UUID = {
    2: "GPU-8277e634-68d1-03db-6e08-d95cf37c2e3f",
    4: "GPU-cef1d233-ecad-b58f-0e2f-099515f265cf",
    6: "GPU-8cac1199-badc-4af1-1454-ed227d4a52f7",
    7: "GPU-2c5fe720-975d-230c-a345-acd352747a1d",
}


def load_manifest():
    return json.loads(MANIFEST.read_text())


def task_dir(task):
    return RUN_ROOT / task["system"] / task["name"]


def check_box_size(rst_path):
    """Check if restart file has a valid box (< 300 Å) using parmed."""
    try:
        import parmed as pmd
        rst = pmd.load_file(str(rst_path))
        box = rst.box
        if box is None:
            return 0, True  # no box = assume OK
        box_a = float(box[0])
        return box_a, box_a < 300.0
    except Exception:
        return 0, True  # can't check = assume OK


def setup_task(task):
    """Create the run directory and prepare restart files.

    Strategy:
    1. Check if backup restart (original pmemd output with velocities) has valid box
    2. If valid: use it as md_seg0.rst with irest=1, ntx=5 (read velocities)
    3. If invalid (box > 300 Å): extract from trajectory with autoimage,
       run 1ns equilibration with irest=0, ntx=1, then use that for production
    """
    import shutil
    d = task_dir(task)
    d.mkdir(parents=True, exist_ok=True)

    # symlink prmtop
    prmtop_link = d / "system.prmtop"
    if prmtop_link.exists() or prmtop_link.is_symlink():
        prmtop_link.unlink()
    src_prmtop = Path(task["src_dir"]) / "system.prmtop"
    prmtop_link.symlink_to(src_prmtop)

    # Check backup restart (original pmemd output with velocities)
    backup = Path(task["src_dir"]) / "md_1.rst.overwritten_backup"
    rst_src = Path(task["restart"])  # current md_1.rst (cpptraj-regenerated)

    rst_dst = d / "md_seg0.rst"
    if rst_dst.exists() or rst_dst.is_symlink():
        rst_dst.unlink()

    needs_equil = False

    if backup.exists():
        box_a, box_ok = check_box_size(backup)
        if box_ok:
            # Use the backup (has velocities, correct box)
            shutil.copy2(backup, rst_dst)
            print(f"  [setup] {task['name']}: using backup restart (box={box_a:.1f} Å, irest=1)")
        else:
            # Backup box is too large (from previous buggy run)
            print(f"  [setup] {task['name']}: backup box={box_a:.1f} Å too large, need equilibration")
            needs_equil = True
    else:
        # No backup, use cpptraj-extracted restart
        needs_equil = True

    if needs_equil:
        # Extract last frame with autoimage to fix wrapped coordinates
        # Then run 1ns equilibration before production
        auto_rst = d / "md_seg0_autoimage.rst"
        if not auto_rst.exists():
            traj = Path(task["src_dir"]) / "md_total.nc"
            inp = f"""parm {src_prmtop}
trajin {traj} lastframe
autoimage
trajout {auto_rst} restart
run
quit
"""
            inp_file = d / "autoimage_extract.in"
            inp_file.write_text(inp)
            logf = LOG_ROOT / f"{task['system']}_{task['name']}_autoimage.log"
            logf.parent.mkdir(parents=True, exist_ok=True)
            with open(logf, "w") as lf:
                subprocess.run(
                    [CPPTRAJ, "-i", str(inp_file)],
                    stdout=lf, stderr=subprocess.STDOUT,
                    env={"AMBERHOME": AMBERHOME, **os.environ},
                )
            if auto_rst.exists() and auto_rst.stat().st_size > 1000:
                shutil.copy2(auto_rst, rst_dst)
                print(f"  [setup] {task['name']}: using autoimage restart (needs 1ns equilibration)")
            else:
                # Fallback: use cpptraj-extracted restart
                shutil.copy2(rst_src, rst_dst)
                print(f"  [setup] {task['name']}: fallback to cpptraj restart")

    # Write md input (irest=1, ntx=5 for production)
    inp_dst = d / "md_extend_100ns.in"
    inp_dst.write_text(MD_INPUT.read_text())

    # Write equilibration input if needed
    if needs_equil:
        equil_dst = d / "md_equil_1ns.in"
        equil_dst.write_text(EQUIL_INPUT.read_text())

    # write status
    status = d / "status.json"
    if not status.exists():
        status.write_text(json.dumps({
            "task": task["name"],
            "system": task["system"],
            "tier": task["tier"],
            "target_ns": task["target_ns"],
            "n_segments": task["n_segments"],
            "completed_segments": 0,
            "status": "pending",
            "gpu": None,
            "needs_equil": needs_equil,
        }, indent=2))

    return d


def run_segment(task, seg_num, gpu):
    """Run a single 100ns segment on the specified GPU.

    For seg1: if needs_equil, run 1ns equilibration first (irest=0,ntx=1),
    then use the equilibrated restart for production (irest=1,ntx=5).
    For seg2+: use pmemd's own restart output directly (no cpptraj extraction).
    """
    d = task_dir(task)
    inp = d / "md_extend_100ns.in"

    # input restart
    if seg_num == 1:
        crd = d / "md_seg0.rst"
    else:
        crd = d / f"md_seg{seg_num-1}.rst"

    out = d / f"md_seg{seg_num}.out"
    rst = d / f"md_seg{seg_num}.rst"
    nc = d / f"md_seg{seg_num}.nc"

    # skip if already completed
    if out.exists() and rst.exists() and not rst.is_symlink():
        with open(out, "r", errors="ignore") as f:
            tail = f.read()[-2000:]
        if "STOP" in tail or "Final" in tail:
            print(f"  [skip] segment {seg_num} already completed for {task['name']}")
            return 0

    env = os.environ.copy()
    env["AMBERHOME"] = AMBERHOME
    if os.environ.get("MD_EXTEND_NO_UUID"):
        env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    else:
        env["CUDA_VISIBLE_DEVICES"] = GPU_UUID.get(gpu, str(gpu))
    env["LD_LIBRARY_PATH"] = f"{AMBERHOME}/lib:{env.get('LD_LIBRARY_PATH','')}"

    # If seg1 and needs_equil, run 1ns equilibration first
    if seg_num == 1:
        status_path = d / "status.json"
        if status_path.exists():
            status = json.loads(status_path.read_text())
            if status.get("needs_equil", False):
                equil_inp = d / "md_equil_1ns.in"
                equil_out = d / "md_equil.out"
                equil_rst = d / "md_equil.rst"
                equil_nc = d / "md_equil.nc"
                if not (equil_out.exists() and equil_rst.exists()):
                    print(f"  [equil] {task['name']}: running 1ns equilibration (irest=0)")
                    equil_cmd = [
                        PMEMD, "-O",
                        "-i", str(equil_inp),
                        "-o", str(equil_out),
                        "-p", str(d / "system.prmtop"),
                        "-c", str(crd),
                        "-r", str(equil_rst),
                        "-x", str(equil_nc),
                    ]
                    equil_log = LOG_ROOT / f"{task['system']}_{task['name']}_equil_gpu{gpu}.log"
                    with open(equil_log, "w") as lf:
                        proc = subprocess.run(equil_cmd, env=env, stdout=lf, stderr=subprocess.STDOUT)
                    if proc.returncode != 0:
                        print(f"  [FAIL] equilibration failed rc={proc.returncode}")
                        return proc.returncode
                    # Use equilibrated restart for production
                    crd = equil_rst
                    print(f"  [equil] done, starting production from equilibrated state")
                else:
                    crd = equil_rst
                    print(f"  [equil] already done, continuing")

    cmd = [
        PMEMD, "-O",
        "-i", str(inp),
        "-o", str(out),
        "-p", str(d / "system.prmtop"),
        "-c", str(crd),
        "-r", str(rst),
        "-x", str(nc),
    ]

    logf = LOG_ROOT / f"{task['system']}_{task['name']}_seg{seg_num}_gpu{gpu}.log"
    logf.parent.mkdir(parents=True, exist_ok=True)

    print(f"  [run] {task['name']} seg{seg_num}/{task['n_segments']} on GPU{gpu}")

    with open(logf, "w") as lf:
        proc = subprocess.run(cmd, env=env, stdout=lf, stderr=subprocess.STDOUT)

    if proc.returncode != 0:
        print(f"  [FAIL] {task['name']} seg{seg_num} rc={proc.returncode}")
        return proc.returncode

    # update status (no cpptraj extraction — use pmemd's own restart)
    status_path = d / "status.json"
    status = json.loads(status_path.read_text())
    status["completed_segments"] = seg_num
    status["gpu"] = gpu

    # Box guard: detect NPT swelling before starting the next segment
    box_a, _ = check_box_size(rst)
    if box_a >= BOX_STOP_A:
        status["status"] = "unstable"
        status["box_a"] = box_a
        status_path.write_text(json.dumps(status, indent=2))
        print(f"  [GUARD] {task['name']} seg{seg_num} box={box_a:.0f} Å >= {BOX_STOP_A:.0f}, "
              f"marking UNSTABLE and stopping task")
        return 99
    if box_a >= BOX_WARN_A:
        status["box_warn"] = True
        print(f"  [WARN] {task['name']} seg{seg_num} box={box_a:.0f} Å >= {BOX_WARN_A:.0f}")

    status["status"] = "running" if seg_num < task["n_segments"] else "completed"
    status["box_a"] = box_a
    status_path.write_text(json.dumps(status, indent=2))

    print(f"  [ok]  {task['name']} seg{seg_num} done")
    return 0


def run_task(task, gpu):
    """Run all segments of a task sequentially on one GPU."""
    if task.get("remote"):
        print(f"\n[skip-remote] {task['system']}/{task['name']} "
              f"-- delegated to remote server")
        return 0
    d = setup_task(task)
    print(f"\n[start] {task['system']}/{task['name']} on GPU{gpu} "
          f"({task['n_segments']} segments, target={task['target_ns']}ns)")

    for seg in range(1, task["n_segments"] + 1):
        rc = run_segment(task, seg, gpu)
        if rc != 0:
            if rc == 99:
                # box guard already marked task unstable; keep that status
                print(f"  [ABORT] {task['name']} stopped by box guard at segment {seg}")
            else:
                print(f"  [ABORT] {task['name']} failed at segment {seg}")
                status_path = d / "status.json"
                status = json.loads(status_path.read_text())
                status["status"] = "failed"
                status_path.write_text(json.dumps(status, indent=2))
            return rc

    print(f"  [DONE] {task['name']} completed all {task['n_segments']} segments")
    return 0


def show_status():
    """Show current status of all tasks."""
    manifest = load_manifest()
    tasks = manifest["tasks"]

    print(f"\n{'='*90}")
    print(f"{'GPU':>4} | {'System':10s} | {'Name':35s} | {'Tier':7s} | {'Segs':5s} | {'Status':12s} | {'Drift':>8s}")
    print(f"{'-'*90}")

    for t in tasks:
        d = task_dir(t)
        status_path = d / "status.json"
        if status_path.exists():
            st = json.loads(status_path.read_text())
            completed = st.get("completed_segments", 0)
            total = t["n_segments"]
            status_str = f"{completed}/{total} done"
            if st.get("status") == "failed":
                status_str += " FAILED"
            elif st.get("status") == "completed":
                status_str = "COMPLETE"
        else:
            status_str = "not started"

        print(f"GPU{t['gpu']:>2} | {t['system']:10s} | {t['name']:35s} | {t['tier']:7s} | {t['n_segments']:>4} | {status_str:12s} | {t['drift']:+.3f}")

    print(f"{'='*90}")
    print(f"Total: {len(tasks)} tasks | GPUs: {manifest['gpus']}")


def reset_task(task):
    """Reset a failed/unstable task back to pending (keep completed segments)."""
    d = task_dir(task)
    status_path = d / "status.json"
    if not status_path.exists():
        print(f"  [reset] {task['name']}: not started, nothing to reset")
        return
    st = json.loads(status_path.read_text())
    if st.get("status") in ("failed", "unstable"):
        old = st["status"]
        st["status"] = "pending"
        status_path.write_text(json.dumps(st, indent=2))
        print(f"  [reset] {task['name']}: {old} -> pending "
              f"(completed_segments={st.get('completed_segments', 0)} kept)")
    else:
        print(f"  [reset] {task['name']}: status={st.get('status')}, untouched")


def main():
    parser = argparse.ArgumentParser(description="MD extension runner")
    parser.add_argument("--status", action="store_true", help="show task status and exit")
    parser.add_argument("--reset", action="store_true",
                        help="reset ALL failed/unstable tasks to pending and exit")
    parser.add_argument("--reset-task", type=int, metavar="IDX",
                        help="reset single task by index to pending and exit")
    parser.add_argument("--gpu", type=int, help="run on specific GPU (for single task)")
    parser.add_argument("--task", type=int, help="task index to run (0-based)")
    args = parser.parse_args()

    if args.status:
        show_status()
        return

    manifest = load_manifest()
    tasks = manifest["tasks"]

    if args.reset:
        for t in tasks:
            reset_task(t)
        return

    if args.reset_task is not None:
        reset_task(tasks[args.reset_task])
        return

    if args.gpu is not None and args.task is not None:
        t = tasks[args.task]
        run_task(t, args.gpu)
        return

    # Full dispatch: assign tasks to GPUs round-robin, run concurrently
    # Group tasks by GPU
    gpu_tasks: dict[int, list] = {g: [] for g in manifest.get("gpus", GPUS)}
    for t in tasks:
        gpu_tasks[t["gpu"]].append(t)

    # Launch each GPU's task queue as a background subprocess
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    LOG_ROOT.mkdir(parents=True, exist_ok=True)

    pids = {}
    for gpu, gpu_task_list in gpu_tasks.items():
        if not gpu_task_list:
            continue
        # Build a command that runs all tasks for this GPU sequentially.
        # Use ';' (not '&&') so one failed/unstable task does NOT block
        # the remaining queued tasks on this GPU.
        task_indices = [tasks.index(t) for t in gpu_task_list]
        script = f"cd {RUN_ROOT}"
        for idx in task_indices:
            script += (f" ; python3 {Path(__file__).resolve()} --reset-task {idx}"
                       f" && python3 {Path(__file__).resolve()} --gpu {gpu} --task {idx}")
        logf = LOG_ROOT / f"gpu{gpu}_dispatcher.log"

        with open(logf, "w") as lf:
            proc = subprocess.Popen(
                ["bash", "-c", script],
                stdout=lf,
                stderr=subprocess.STDOUT,
                env={**os.environ, "AMBERHOME": AMBERHOME},
            )
        pids[gpu] = proc.pid
        print(f"[dispatch] GPU{gpu} -> PID {proc.pid}, {len(gpu_task_list)} tasks: {[t['name'] for t in gpu_task_list]}")

    print(f"\n[launched] {len(pids)} GPU workers. Monitor with:")
    print(f"  python3 {Path(__file__).resolve()} --status")
    print(f"  tail -f {LOG_ROOT}/gpu*_dispatcher.log")
    print(f"  nvidia-smi -l 10")


if __name__ == "__main__":
    main()
