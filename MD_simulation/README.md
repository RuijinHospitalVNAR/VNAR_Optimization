# MD 模拟协议（SH3_VNAR / HCG_VNAR 体系）

本文档记录 SH3 与 HCG 两大抗原-VNAR 复合体体系（含全部突变体）的分子动力学模拟完整协议与实际运行命令。

## 体系构建（tleap）

- 力场：`ff14SB`（蛋白质），水模型：`TIP3P`
- 溶剂化：八面体盒子（`solvateoct`），缓冲距离 8.0 Å
- 电荷中和：`addions2` 添加 Na⁺/Cl⁻ 至净电荷为 0

```bash
tleap -f tleap.in
```

`tleap.in`（SH3 与 HCG 相同）：

```
source leaprc.protein.ff14SB
source leaprc.water.tip3p
COM = loadpdb system.pdb
solvateoct COM TIP3PBOX 8.0
charge COM
addions2 COM Na+ 0
addions2 COM Cl- 0
saveAmberParm COM system.prmtop system.inpcrd
savepdb COM system.pdb
quit
```

## 阶段一：预平衡 + 100 ns 生产（NPT）

实际执行命令（`run_single.sh` / `md.sh` 执行链，SH3 与 HCG 一致）：

```bash
# 1. 两阶段能量最小化
pmemd.cuda -O -i min1.in     -o min1.out     -p system.prmtop -c system.inpcrd -r min1.rst -ref system.inpcrd
pmemd.cuda -O -i min2.in     -o min2.out     -p system.prmtop -c min1.rst     -r min2.rst
# 2. NVT 升温
pmemd.cuda -O -i heat.in     -o heat.out     -p system.prmtop -c min2.rst     -r heat.rst -x heat.mdcrd -ref min2.rst -e heat.mden
# 3. NPT 加压预平衡（约束保留）
pmemd.cuda -O -i pressure.in -o pressure.out -p system.prmtop -c heat.rst     -r pres.rst -x pres.mdcrd -ref heat.rst -e pres.mden
# 4. NPT 无约束平衡
pmemd.cuda -O -i equil.in    -o equil.out    -p system.prmtop -c pres.rst     -r equil.rst -x equil.mdcrd -ref pres.rst -e equil.mden
# 5. 100 ns NPT 生产
pmemd.cuda -O -i md.in       -o md_1.out     -p system.prmtop -c equil.rst    -r md_1.rst  -x md_1.nc    -ref equil.rst -e md_1.mden
```

### 各阶段参数（对应 `mdin/` 目录文件）

| 阶段 | 文件 | 关键参数 |
|------|------|----------|
| min1 | `min1.in` | `imaxcyc=15000`（`ncyc=5000`：5000 步最陡下降 + 10000 步共轭梯度），溶质重原子位置约束 `restraint_wt=200.0` kcal·mol⁻¹·Å⁻²，`ntb=1` |
| min2 | `min2.in` | `maxcyc=15000`（5000 SD + 10000 CG），无约束 |
| heat | `heat.in` | 100 ps NVT（`ntb=1`），0→100→200→300 K 分段升温（前 60 ps），后 40 ps 恒温 300 K；Langevin 恒温 `ntt=3, gamma_ln=2.0`；溶质重原子约束 `restraint_wt=2.0` kcal·mol⁻¹·Å⁻² |
| pressure | `pressure.in` | 100 ps NPT（`ntb=2, ntp=1, taup=2.0`），300 K，约束保留 `restraint_wt=2.0` |
| equil | `equil.in` | 100 ps NPT（`ntb=2, ntp=1, taup=2.0`），300 K，**无约束** |
| 生产 | `md.in` | 100 ns NPT，`nstlim=50000000, dt=0.002`，`ntc=2, ntf=2`（SHAKE），`cut=8.0`，Langevin `gamma_ln=2.0`，`iwrap=1` |

## 阶段二：400 ns 延伸（4 × 100 ns，NVT）

100 ns 生产完成后，用 `06_run_md_extend.py` 调度追加 4 段 × 100 ns。

> **系综说明**：延伸段采用 **NVT（`ntb=1`）** 而非 NPT——初始 NPT 延伸在部分体系出现盒子膨胀/不稳定（见 `tasks_manifest.json` 中 box-guard 阈值与 failed 任务），故统一改用 NVT 常温常体积延伸。预平衡已完成充分的 NPT 弛豫，NVT 延伸是稳定可行的替代方案。

```bash
# 调度器（按 tasks_manifest.json 分配 GPU，段间自动衔接）
python3 06_run_md_extend.py                  # 运行全部任务
python3 06_run_md_extend.py --status         # 查看进度
python3 06_run_md_extend.py --gpu 2 --task 0 # 指定 GPU 运行指定任务
```

每段实际命令（第 K 段，前一段 restart 衔接）：

```bash
pmemd.cuda -O -i mdin_nvt_100ns.in -o md_segK.out -p system.prmtop \
  -c md_seg(K-1).rst -r md_segK.rst -x md_segK.nc
```

`mdin_nvt_100ns.in`：100 ns NVT，`nstlim=50000000, dt=0.002`，`ntb=1`，`ntt=3, gamma_ln=2.0`，`ntc=2, ntf=2`，`cut=8.0`，`iwrap=1`。

段 1 从原 100 ns 生产末帧 `md_1.rst`（或经 box 校验的备份 restart）restart；段间通过 `md_segK.rst` 无缝续接，共 5 × 100 ns = 500 ns。

## 任务清单

- **SH3_VNAR**：WT_original、S85I_G91D、S85T、S85T_H95Y、S85I_G91D_H95Y、S85M_E96S、S86G_G91S、S86R_Y92D、Y84A_S86G、E96V 等
- **HCG_VNAR**：I93V_S92Y、I93V_M96S_S85T、I93V_M96S_S85V、I93V_M96S_S92Y、I93V_M96S_H98D 等

完整清单与分层（high/medium）见 `tasks_manifest.json`。

## 文件说明

```
MD_simulation/
├── README.md               ← 本文档
├── run_single.sh           ← 阶段一全流程执行脚本（tleap → min → heat → pressure → equil → 生产）
├── 06_run_md_extend.py     ← 阶段二 500 ns 延伸调度器（GPU 分配 + 段衔接 + box-guard）
└── mdin/                   ← 全部 AMBER 输入文件（实际运行所用版本）
    ├── min1.in / min2.in   ← 两阶段最小化
    ├── heat.in             ← NVT 升温
    ├── pressure.in         ← NPT 加压（约束保留）
    ├── equil.in            ← NPT 无约束平衡
    ├── md.in               ← 100 ns NPT 生产
    └── mdin_nvt_100ns.in   ← 延伸段 100 ns NVT
```

## 结合自由能

生产轨迹经 cpptraj 去水处理后，用 MMPBSA.py 计算 MM/GBSA（`igb=5, saltcon=0.154`）；500 ns 收官后以统一协议（末 20 ns × 100 帧）重算全部体系并添加 per-residue decomposition。
