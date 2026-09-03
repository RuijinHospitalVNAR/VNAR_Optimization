#!/usr/bin/env bash
# AMBER 单结构 MD：tleap -> min1 -> min2 -> heat -> pressure -> equil -> md -> cpptraj -> MMPBSA
# 续跑时写入新段 md_2.*, md_3.* ... 不覆写前段，可后续用 cpptraj 拼接。
# 用法: run_single.sh --structure <pdb|cif> --output_dir <dir> [options]
#
# Options:
#   --gpu_id N            GPU 编号（仅用于日志，实际设备由 CUDA_VISIBLE_DEVICES 控制）
#   --resume              从已有阶段续跑
#   --production_ns N     Production 时长（ns），默认 100
#   --npt_ns N            NPT 预平衡总时长（ns），默认 0（使用最小默认值：heat 0.1ns + pressure 0.1ns + equil 0.1ns）
#   --forcefield NAME     力场名称，默认 ff19SB（支持 ff19SB, ff14SB, amber99sb）
#   --buffer ANGSTROM     Solvate buffer 距离，默认 8.0
#   --receptor_mask MASK  MMPBSA receptor mask（可选，如 :1-105）
#   --ligand_mask MASK    MMPBSA ligand mask（可选，如 :106-211）
#   --interval N          MMPBSA 取帧间隔，默认 5
#   --skip_mmpbsa         跳过 MMPBSA 后处理，立即返回（用于 GPU 复用场景）

set -e

AMBER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
structure=""
output_dir=""
gpu_id=0
resume=0
production_ns=100
npt_ns=0
forcefield="ff14SB"
buffer="8.0"
receptor_mask=""
ligand_mask=""
interval=5
skip_mmpbsa=0

show_usage() {
    echo "Usage: $0 --structure <pdb|cif> --output_dir <dir> [--gpu_id N] [--resume] [--production_ns N] [--npt_ns N] [--forcefield NAME] [--buffer ANGSTROM] [--receptor_mask MASK] [--ligand_mask MASK] [--interval N]"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --structure) structure="$2"; shift 2 ;;
        --output_dir) output_dir="$2"; shift 2 ;;
        --gpu_id) gpu_id="$2"; shift 2 ;;
        --resume) resume=1; shift 1 ;;
        --production_ns) production_ns="$2"; shift 2 ;;
        --npt_ns) npt_ns="$2"; shift 2 ;;
        --forcefield) forcefield="$2"; shift 2 ;;
        --buffer) buffer="$2"; shift 2 ;;
        --receptor_mask) receptor_mask="$2"; shift 2 ;;
        --ligand_mask) ligand_mask="$2"; shift 2 ;;
        --interval) interval="$2"; shift 2 ;;
        --skip_mmpbsa) skip_mmpbsa=1; shift 1 ;;
        -h|--help) show_usage; exit 0 ;;
        *) echo "Unknown option: $1"; show_usage; exit 1 ;;
    esac
done

if [[ -z "$structure" || -z "$output_dir" ]]; then
    echo "Error: --structure and --output_dir are required."
    show_usage
    exit 1
fi

output_dir="$(realpath "$output_dir")"
structure="$(realpath -s "$structure")"
mkdir -p "$output_dir"
cd "$output_dir"

echo "[GPU-MAP] run_single gpu_id_arg=${gpu_id} CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}"

# ---------- 力场/水模型映射 ----------
FORCEFIELD_LEAP="$forcefield"
WATERMODEL="tip3p"
WATERBOX="TIP3PBOX"
case "$forcefield" in
    ff19SB|ff19sb)
        FORCEFIELD_LEAP="ff19SB"
        WATERMODEL="opc"
        WATERBOX="OPCBOX"
        ;;
    ff14SB|ff14sb)
        FORCEFIELD_LEAP="ff14SB"
        WATERMODEL="tip3p"
        WATERBOX="TIP3PBOX"
        ;;
    amber99sb| Amber99SB|ff99SB|ff99sb)
        FORCEFIELD_LEAP="amber99sb"
        WATERMODEL="tip3p"
        WATERBOX="TIP3PBOX"
        ;;
esac
echo "[INFO] Forcefield: ${FORCEFIELD_LEAP}, Water: ${WATERMODEL}, Box: ${WATERBOX}, Buffer: ${buffer} Å"

# ---------- PDB 准备 ----------
PDB_NAME="system.pdb"
if [[ "$structure" == *.cif || "$structure" == *.mmcif ]]; then
    echo "[INFO] Converting CIF to PDB..."
    if head -n 2 "$structure" 2>/dev/null | grep -Eq '^(HEADER|ATOM|HETATM|REMARK|CRYST1)'; then
        cp -f "$structure" "$PDB_NAME" 2>/dev/null || ln -sf "$structure" "$PDB_NAME"
    else
        obabel -icif "$structure" -opdb -O "$PDB_NAME" 2>/dev/null || {
            python3 -c "
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO
p = MMCIFParser(QUIET=True)
s = p.get_structure('X', '$structure')
io = PDBIO()
io.set_structure(s)
io.save('$PDB_NAME')
"
        }
    fi
else
    cp -f "$structure" "$PDB_NAME" 2>/dev/null || ln -sf "$structure" "$PDB_NAME"
fi

if [[ ! -f "$PDB_NAME" ]]; then
    echo "[ERROR] Failed to prepare $PDB_NAME from $structure"
    exit 1
fi

echo "[INFO] pdb4amber (fix CYS -> CYX for disulfide bonds)..."
pdb4amber "$PDB_NAME" --add-missing-atoms -y -o "${PDB_NAME}.p4a" 2>/dev/null && mv -f "${PDB_NAME}.p4a" "$PDB_NAME" || true
awk 'BEGIN{OFMT="%.3f"}
$0~/^ATOM/||$0~/^HETATM/{
  el=substr($0,77,2); gsub(/ /,"",el);
  if(el=="H") next;
}
{print $0}' "$PDB_NAME" > "${PDB_NAME}.noH" && mv -f "${PDB_NAME}.noH" "$PDB_NAME"

if ! grep -qE '^(ATOM|HETATM)' "$PDB_NAME" 2>/dev/null; then
    echo "[ERROR] $PDB_NAME has no ATOM/HETATM"
    exit 2
fi

# ---------- tleap ----------
if [[ ! -f system.prmtop || ! -f system.inpcrd ]] || [[ $resume -eq 0 ]]; then
    # 二硫键条件: ff19SB需要加载额外参数
    DISULFIDE_PARAM=""
    if [[ "$forcefield" == "ff19SB" ]]; then
        DISULFIDE_PARAM="loadamberparams frcmod.ff19SB"
    fi
    
    sed -e "s/{{FORCEFIELD}}/${FORCEFIELD_LEAP}/g" \
        -e "s/{{WATERMODEL}}/${WATERMODEL}/g" \
        -e "s/{{WATERBOX}}/${WATERBOX}/g" \
        -e "s/{{BUFFER}}/${buffer}/g" \
        -e "s/{{DISULFIDE_PARAM}}/${DISULFIDE_PARAM}/g" \
        -e "s/STRUCTURE_PDB/${PDB_NAME}/g" \
        "${AMBER_DIR}/tleap.in.template" > tleap.in
    echo "[INFO] Running tleap..."
    tleap -f tleap.in
fi

if [[ ! -f system.prmtop || ! -f system.inpcrd ]]; then
    echo "[ERROR] tleap did not produce system.prmtop/system.inpcrd"
    exit 3
fi

PRMTOP=system.prmtop

# ---------- 输入文件准备 ----------
# min1/min2 不变
for f in min1.in min2.in; do
    [[ -f "${AMBER_DIR}/$f" && ! -f "$f" ]] && cp "${AMBER_DIR}/$f" .
done

# npt_ns 控制 heat/pressure/equil 步数
if [[ "$npt_ns" -gt 0 ]]; then
    heat_steps=$(( npt_ns * 100000 ))
    [[ "$heat_steps" -lt 50000 ]] && heat_steps=50000
    pressure_steps=$(( npt_ns * 200000 ))
    [[ "$pressure_steps" -lt 50000 ]] && pressure_steps=50000
    equil_steps=$(( npt_ns * 200000 ))
    [[ "$equil_steps" -lt 50000 ]] && equil_steps=50000

    echo "[INFO] NPT total ${npt_ns} ns -> heat ${heat_steps} steps, pressure ${pressure_steps} steps, equil ${equil_steps} steps"

    # 生成 heat.in（含动态升温程序）
    python3 - <<PY
nstlim = ${heat_steps}
p1 = int(nstlim * 0.2)
p2 = int(nstlim * 0.4)
p3 = int(nstlim * 0.6)
p4 = nstlim
with open("heat.in", "w") as f:
    f.write(f"""heat
 &cntrl
  imin=0,irest=0,ntx=1,
  nstlim={nstlim},dt=0.002,
  ntc=2,ntf=2,
  cut=8.0, ntb=1,
  ntpr=1000, ntwx=1000,
  ntwr=1000, ntwe=1000,
  ntt=3, gamma_ln=2.0,
  tempi=0.0, temp0=300.0, ig=-1,
  ntr=1, restraintmask=':1-99999 & !@H=',
  restraint_wt=2.0,
  nmropt=1,vlimit=20,
 &end
  &wt TYPE='TEMP0', istep1=0, istep2={p1},
  value1=0.0, value2=100.0,
  &end
  &wt TYPE='TEMP0', istep1={p1}, istep2={p2},
  value1=100.0, value2=200.0,
  &end
  &wt TYPE='TEMP0', istep1={p2}, istep2={p3},
  value1=200.0, value2=300.0,
  &end
  &wt TYPE='TEMP0', istep1={p3}, istep2={p4},
  value1=300.0, value2=300.0,
  &end
 &wt TYPE='END' ,
&end
""")
PY

    # pressure / equil 只改 nstlim
    sed "s/nstlim=50000/nstlim=${pressure_steps}/" "${AMBER_DIR}/pressure.in" > pressure.in
    sed "s/nstlim=50000/nstlim=${equil_steps}/" "${AMBER_DIR}/equil.in" > equil.in
else
    # 使用默认最短预平衡
    for f in heat.in pressure.in equil.in; do
        [[ -f "${AMBER_DIR}/$f" && ! -f "$f" ]] && cp "${AMBER_DIR}/$f" .
    done
fi

# md.in：根据 production_ns 动态调整 nstlim 和输出频率
DT=0.002
TARGET_PS=$(( production_ns * 1000 ))
nstlim_first=$(( TARGET_PS * 1000 / 2 ))  # ps / dt，dt=0.002 -> 1000ps / 0.002 = 500000

# 输出频率：ns 越长，频率越低，节省 I/O
if [[ "$production_ns" -le 1 ]]; then
    ntwx=1000; ntpr=1000
elif [[ "$production_ns" -le 10 ]]; then
    ntwx=2500; ntpr=2500
else
    ntwx=5000; ntpr=5000
fi
cat > md.in <<EOF
Production ${production_ns}ns
 &cntrl
  imin=0,irest=1,ntx=5,
  nstlim=${nstlim_first},dt=${DT},
  ntc=2,ntf=2,
  cut=8.0, ntb=2, ntp=1, taup=2.0,
  ntpr=${ntpr}, ntwx=${ntwx},
  ntwr=${ntwx}, ntwe=${ntwx},
  ntt=3, gamma_ln=2.0,
  tempi=300.0, temp0=300.0, ig=-1,
  iwrap=1
 &end
EOF

echo "[INFO] Production: ${production_ns} ns (${nstlim_first} steps), ntwx=${ntwx}"

# ---------- 运行函数 ----------
run_min1() {
    pmemd.cuda -O -i min1.in -o min1.out -p $PRMTOP -c system.inpcrd -r min1.rst -ref system.inpcrd
}
run_min2() {
    pmemd.cuda -O -i min2.in -o min2.out -p $PRMTOP -c min1.rst -r min2.rst
}
run_heat() {
    pmemd.cuda -O -i heat.in -o heat.out -p $PRMTOP -c min2.rst -r heat.rst -x heat.mdcrd -ref min2.rst -e heat.mden
}
run_pressure() {
    pmemd.cuda -O -i pressure.in -o pressure.out -p $PRMTOP -c heat.rst -r pres.rst -x pres.mdcrd -ref heat.rst -e pres.mden
}
run_equil() {
    pmemd.cuda -O -i equil.in -o equil.out -p $PRMTOP -c pres.rst -r equil.rst -x equil.mdcrd -ref pres.rst -e equil.mden
}

run_md() {
    local inpcrd="$1"
    local ref="$2"
    local prefix="$3"
    local nstlim_opt="$4"
    local mdin="md.in"
    if [[ -n "$nstlim_opt" ]]; then
        mdin="md_continue.in"
        sed "s/nstlim=${nstlim_first}/nstlim=${nstlim_opt}/" "${AMBER_DIR}/md.in" > "$mdin"
    fi
    pmemd.cuda -O -i "$mdin" -o "${prefix}.out" -p $PRMTOP -c "$inpcrd" -r "${prefix}.rst" -x "${prefix}.nc" -ref "$ref" -e "${prefix}.mden"
}

get_latest_segment() {
    local n=0
    for f in md_*.rst; do
        [[ -f "$f" ]] || continue
        local i="${f#md_}"
        i="${i%.rst}"
        [[ "$i" =~ ^[0-9]+$ ]] && [[ "$i" -gt "$n" ]] && n="$i"
    done
    echo "$n"
}

get_last_time_ps() {
    local f="$1"
    local last=""
    [[ -f "$f" ]] || { echo ""; return; }
    while IFS= read -r line; do
        if [[ "$line" =~ TIME[[:space:]]*\([[:space:]]*PS[[:space:]]*\)[[:space:]]*=[[:space:]]*([0-9]+\.?[0-9]*) ]]; then
            last="${BASH_REMATCH[1]}"
        fi
    done < "$f"
    echo "$last"
}

get_max_nstep() {
    local f="$1"
    local max=""
    [[ -f "$f" ]] || { echo ""; return; }
    while IFS= read -r line; do
        if [[ "$line" =~ NSTEP[[:space:]]*=[[:space:]]*([0-9]+) ]]; then
            local nstep="${BASH_REMATCH[1]}"
            if [[ -z "$max" ]] || [[ "$nstep" -gt "$max" ]]; then
                max="$nstep"
            fi
        fi
    done < "$f"
    echo "$max"
}

# ---------- Resume 逻辑 ----------
SKIP_MD=0
if [[ $resume -eq 1 ]]; then
    latest=$(get_latest_segment)
    if [[ "$latest" -ge 1 ]]; then
        last_ps=$(get_last_time_ps "md_${latest}.out")
        max_nstep=$(get_max_nstep "md_${latest}.out")
        if [[ -n "$last_ps" ]]; then
            target_nstep=$(awk -v ps="$TARGET_PS" -v dt="$DT" 'BEGIN { printf "%.0f", ps / dt * 0.98 }')
            time_ok=$(awk -v t="$last_ps" -v g="$TARGET_PS" 'BEGIN { exit !(t >= g) }' 2>/dev/null && echo "1" || echo "0")
            step_ok="0"
            if [[ -n "$max_nstep" ]] && [[ "$max_nstep" -ge "$target_nstep" ]]; then
                step_ok="1"
            fi
            if [[ "$time_ok" == "1" ]] && [[ "$step_ok" == "1" ]]; then
                echo "[RESUME] Production 已完成（${last_ps} ps, ${max_nstep} steps >= ${TARGET_PS} ps, ${target_nstep} steps），跳过 MD。"
                SKIP_MD=1
            else
                if [[ -n "$max_nstep" ]]; then
                    target_nstep_full=$(awk -v ps="$TARGET_PS" -v dt="$DT" 'BEGIN { printf "%.0f", ps / dt }')
                    nstlim_remain=$((target_nstep_full - max_nstep))
                    if [[ "$nstlim_remain" -le 0 ]]; then
                        remain_ps=$(awk -v a="$TARGET_PS" -v b="$last_ps" 'BEGIN { printf "%.0f", a - b }')
                        if [[ "$remain_ps" -le 0 ]]; then
                            echo "[RESUME] 已超过目标（${last_ps} ps >= ${TARGET_PS} ps），但步数检查异常，跳过 MD。"
                            SKIP_MD=1
                        else
                            nstlim_remain=$(awk -v r="$remain_ps" -v dt="$DT" 'BEGIN { printf "%.0f", r / dt }')
                        fi
                    fi
                else
                    remain_ps=$(awk -v a="$TARGET_PS" -v b="$last_ps" 'BEGIN { printf "%.0f", a - b }')
                    nstlim_remain=$(awk -v r="$remain_ps" -v dt="$DT" 'BEGIN { printf "%.0f", r / dt }')
                fi

                if [[ "$SKIP_MD" == "0" ]]; then
                    if [[ -z "$nstlim_remain" || "$nstlim_remain" -le 0 ]]; then
                        echo "[RESUME] 剩余步数无效（${nstlim_remain}），跳过 MD。"
                        SKIP_MD=1
                    else
                        remain_ps_display=$(awk -v steps="$nstlim_remain" -v dt="$DT" 'BEGIN { printf "%.1f", steps * dt }')
                        if [[ "$time_ok" == "1" ]] && [[ "$step_ok" == "0" ]]; then
                            echo "[RESUME] 时间已达 ${last_ps} ps，但步数不足（${max_nstep} < ${target_nstep}），续跑剩余 ${nstlim_remain} 步（约 ${remain_ps_display} ps）..."
                        else
                            echo "[RESUME] 从 md_${latest}.rst 续跑（当前 ${last_ps} ps，剩余约 ${remain_ps_display} ps，${nstlim_remain} 步）..."
                        fi
                        next=$((latest + 1))
                        run_md "md_${latest}.rst" "md_${latest}.rst" "md_${next}" "$nstlim_remain"
                        echo "[INFO] AMBER MD 续跑完成（段 md_${next}）: $output_dir"
                        SKIP_MD=1
                    fi
                fi
            fi
        fi
    fi

    if [[ "$SKIP_MD" == "0" ]]; then
        if [[ -f equil.rst && ! -f md_1.rst ]]; then
            echo "[RESUME] 从 equil.rst 开始 Production (md_1)..."
            run_md equil.rst equil.rst md_1 "$nstlim_first"
            echo "[INFO] AMBER MD 续跑完成: $output_dir"
            SKIP_MD=1
        elif [[ -f pres.rst && ! -f equil.rst ]]; then
            echo "[RESUME] 从 pres.rst 续跑 equil + md..."
            run_equil
            run_md equil.rst equil.rst md_1 "$nstlim_first"
            echo "[INFO] AMBER MD 续跑完成: $output_dir"
            SKIP_MD=1
        elif [[ -f heat.rst && ! -f pres.rst ]]; then
            echo "[RESUME] 从 heat.rst 续跑 pressure + equil + md..."
            run_pressure
            run_equil
            run_md equil.rst equil.rst md_1 "$nstlim_first"
            echo "[INFO] AMBER MD 续跑完成: $output_dir"
            SKIP_MD=1
        elif [[ -f min2.rst && ! -f heat.rst ]]; then
            echo "[RESUME] 从 min2.rst 续跑 heat + pressure + equil + md..."
            run_heat
            run_pressure
            run_equil
            run_md equil.rst equil.rst md_1 "$nstlim_first"
            echo "[INFO] AMBER MD 续跑完成: $output_dir"
            SKIP_MD=1
        elif [[ -f min1.rst && ! -f min2.rst ]]; then
            echo "[RESUME] 从 min1.rst 续跑 min2 + heat + pressure + equil + md..."
            run_min2
            run_heat
            run_pressure
            run_equil
            run_md equil.rst equil.rst md_1 "$nstlim_first"
            echo "[INFO] AMBER MD 续跑完成: $output_dir"
            SKIP_MD=1
        elif [[ -f system.prmtop && ! -f min1.rst ]]; then
            echo "[RESUME] 从头开始完整 MD..."
            # 不设 SKIP_MD，让后面的完整流程执行
            :
        fi
    fi
fi

# ---------- 完整 MD 流程 ----------
if [[ "$SKIP_MD" == "0" ]]; then
    echo "[INFO] min1..."
    run_min1
    echo "[INFO] min2..."
    run_min2
    echo "[INFO] heat..."
    run_heat
    echo "[INFO] pressure..."
    run_pressure
    echo "[INFO] equil..."
    run_equil
    echo "[INFO] md (${production_ns} ns) -> md_1.*..."
    run_md equil.rst equil.rst md_1 "$nstlim_first"
    echo "[INFO] AMBER MD 完成: $output_dir"
fi

# ---------- 后处理：cpptraj 去水 + RMSD（始终执行）----------
# 注意：cpptraj 必须跑（生成 md_protein_nowat.nc），即使跳过 MMPBSA
if [[ -x "${AMBER_DIR}/postprocess_single_cpptraj.sh" ]]; then
    echo "[INFO] Running cpptraj post-processing..."
    bash "${AMBER_DIR}/postprocess_single_cpptraj.sh" "$output_dir" || {
        echo "[WARN] cpptraj post-processing failed, will try MMPBSA with raw trajectory"
    }
else
    echo "[WARN] postprocess_single_cpptraj.sh not found"
fi

# ---------- MMPBSA ----------
if [[ "$skip_mmpbsa" == "1" ]]; then
    echo "[INFO] --skip_mmpbsa: cpptraj done, 跳过 MMPBSA 计算，立即返回。"
    exit 0
fi

MMPBSA_DIR="$(cd "$AMBER_DIR/../AMBER_MMPBSA" && pwd)"

# 1) MM/GBSA
if [[ -x "${MMPBSA_DIR}/run_mmpbsa_single_gb.sh" ]]; then
    echo "[INFO] Running MM/GBSA (interval=${interval})..."
    export MMPBSA_INTERVAL="${interval}"
    export MMPBSA_STARTFRAME=1
    export MMPBSA_ENDFRAME=9999999
    if [[ -n "$receptor_mask" && -n "$ligand_mask" ]]; then
        echo "[INFO] Using provided masks: receptor=${receptor_mask} ligand=${ligand_mask}"
        bash "${MMPBSA_DIR}/run_mmpbsa_single_gb.sh" \
            --amber_dir "$output_dir" \
            --receptor_mask "$receptor_mask" \
            --ligand_mask "$ligand_mask"
    else
        bash "${MMPBSA_DIR}/run_mmpbsa_single_gb.sh" --amber_dir "$output_dir"
    fi
else
    echo "[WARN] MMPBSA script not found: ${MMPBSA_DIR}/run_mmpbsa_single_gb.sh"
fi

# 2) 生成 mmgbsa_summary.csv（供 Python 外层读取）
echo "[INFO] Generating mmgbsa_summary.csv..."
python3 - <<'PY'
import os, csv, re

mmgbsa_path = "FINAL_RESULTS_MMGBSA_BINDING.dat"
dG = ""
if os.path.isfile(mmgbsa_path):
    with open(mmgbsa_path) as f:
        for line in f:
            m = re.match(r"DELTA TOTAL\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line.strip())
            if m:
                dG = m.group(1)
                break

rmsd_path = "rmsd_bb.dat"
rmsds = []
if os.path.isfile(rmsd_path):
    with open(rmsd_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    rmsds.append(float(parts[1]))
                except ValueError:
                    pass

mean_nm = ""
max_nm = ""
if rmsds:
    mean_nm = f"{sum(rmsds)/len(rmsds)/10.0:.4f}"
    max_nm = f"{max(rmsds)/10.0:.4f}"

with open("mmgbsa_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["mmgbsa_dG_kcal_mol", "rmsd_mean_nm", "rmsd_max_nm", "error"])
    w.writeheader()
    w.writerow({"mmgbsa_dG_kcal_mol": dG, "rmsd_mean_nm": mean_nm, "rmsd_max_nm": max_nm, "error": ""})
PY

echo "[INFO] AMBER 全流程完成: $output_dir"
