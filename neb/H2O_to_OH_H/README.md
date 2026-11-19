# NEB: H2O* + * -> OH* + H*  (Rh(111) 3x3, WGS step 1)

Follows VASP Tutorial 2, slides 25-36.

| dir  | what it is |
|------|------------|
| `full/neb/` | regular NEB over the full path, IS -> OH* fcc + H* (7 images) |
| `full/cineb/` | climbing-image NEB over the full path (superseded by the split) |
| `split/neb_IS_to_bridge/` | IS -> OH* bridge + H* (the step with the TS): fresh 7-image `nebinit.py` band, regular NEB then CI run 1 |
| `split/ci_IS_to_bridge_damped/` | CI restart of `neb_IS_to_bridge/` with the slides optimizer (damped MD, POTIM 0.01) |
| `split/ci_IS_to_bridge_fire/` | the same CI restart with VTST FIRE (IOPT 7, POTIM 0) |
| `split/neb_bridge_to_fcc/` | OH* bridge + H* -> OH* fcc + H* (3 images); converged, TS vib in `vib/TS_OH_bridge_fcc` |

The IS is `ads/H2O/top` (E0 = -196.598 eV). `IS/` and `cineb_bridge/` (the failed 4-image run) were removed; they are in git history.

**Final state:** OH* fcc + H* fcc co-adsorbed. It's relaxed as its own thermo job in `ads/OH_H`, so run that first.
Its atom order matches IS (Rh x27, O, H that stays on O, H that leaves), which the interpolation needs.

## 0. Scripts (once, from the repo root)
`nebmake.pl` can't be used: it needs the VTST helper `Vasp.pm`, which isn't on Euclid.
Use `nebinit.py` instead, the alternative from slide 27. It uses ASE and does IDPP interpolation, which also avoids the
squeezed O-H bond (0.74 A) that plain linear interpolation gives here.
`scripts/nebinit.py` has a fixed first line and import; copy it over the `~/bin` version:
```bash
cp scripts/nebinit.py ~/bin/
ls /home/vempatti/.conda/envs/ase/bin/python    # the Python it (and kie.py) runs with; must exist
```
`nebef.pl` and `run_neb.sh` work as-is. `nebef.pl`'s column cuts were checked against these VASP 6.4 OUTCARs.
`run_neb.sh` uses 112 cores = 16 per image, which NCORE = 4 divides evenly.

## 1. Generate images and submit (slides 26-31)
`neb/` already holds everything except POTCAR, which git ignores:
`POSCAR_start` (IS CONTCAR), `POSCAR_end` (converged `ads/OH_H` CONTCAR, E0 = -196.3728 eV, max force 0.019),
`OUTCAR_start`/`OUTCAR_end` (`nebinit.py` copies these into 00/ and 08/ for `nebef.pl`), `INCAR`, `KPOINTS` and `run_neb.sh`.
```bash
cd neb/H2O_to_OH_H/neb
cp ../../../ads/H2O/top/POTCAR .
nebinit.py -n 7                        # makes 00-08
ase gui 0*/POSCAR                      # check the path (slide 28's "ag" is an alias for this)
sbatch run_neb.sh
```
Check that `00/OUTCAR` and `08/OUTCAR` exist after `nebinit.py`.

## 2. Restart / climbing image (slides 29, 31)
Restart (the slide's version has curly quotes and also loops over 00/08):
```bash
for d in 0[1-7]; do cp $d/CONTCAR $d/POSCAR; done
sbatch run_neb.sh
```
Once it reaches `reached required accuracy` at EDIFFG = -0.2, set `EDIFFG = -0.02` and `LCLIMB = .TRUE.` in `INCAR`, then restart the same way.

## 3. Analyze (slides 32-36)
```bash
nebef.pl
```
Watch for a local minimum along the path (slide 36). Then run a vib job (`files/INCAR_vib`) on the highest image's CONTCAR. A true TS has exactly one imaginary mode.

## 4. Climbing-image NEB: `cineb/`
`neb/` ran twice (jobs 49988 and 50028, 100 steps each). After the restart, images 01-06 are below the 0.2 eV/A target.
Image 07 isn't: it oscillates between 0.3 and 0.7 eV/A (0.69 at the end).

| image | 00 | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 |
|---|---|---|---|---|---|---|---|---|---|
| E - E(00), eV | 0.00 | 0.03 | 0.71 | **0.85** | 0.61 | 0.29 | 0.29 | 0.29 | 0.22 |
| O-H2, A | 0.98 | 0.99 | 1.24 | 2.15 | 2.51 | 2.86 | 3.13 | 2.93 | 2.91 |
| max force, eV/A | - | 0.06 | 0.09 | 0.14 | 0.10 | 0.02 | 0.16 | 0.69 | - |

The O-H bond breaks between 02 and 03, so 03 is the climbing image. The barrier is about 0.85 eV before CI.
By image 04, H2 is already in its final fcc site. Images 05-07 are only the OH sliding across a flat stretch (all at 0.29 eV) into its fcc site.
Image 07's leftover force is in that flat tail, not at the TS, so the run moved on to CI-NEB rather than restarting again.
The earlier dip at 05 is gone: 05-07 are flat, not a local minimum.

Why 07 won't converge: the band bends sharply at 06-07 (angle about 81 degrees, against 110-140 degrees elsewhere).
The OH overshoots the bridge site at 06 and doubles back toward fcc at 07, and the OH tilt flips between them.
07's force is mostly spring force (0.73 eV/A along the path, from uneven spacing: 0.75 A to 06, 0.60 A to 08). At the corner, the optimizer oscillates instead of settling.

**Check first: `ads/OH_bridge_H`.** It relaxes image 05 (OH* bridge + H* fcc) as an ordinary ads job (`INCAR`/`KPOINTS`/`run.sh` from `ads/OH_H`).
```bash
cd ads/OH_bridge_H
cp ../H2O/top/POTCAR .
sbatch run.sh
```
- If OH stays on the bridge, it's a local minimum (slide 36). Split the NEB: run CI-NEB from IS to this state (that part has the TS), and treat the bridge-to-fcc hop separately.
- If OH rolls into fcc, there's no intermediate. Run `cineb/` as below, but with `EDIFFG = -0.05`, because 07 won't reach 0.02.

`cineb/` starts from `neb/`'s CONTCARs after the restart (01-07). 00/ and 08/ have their OUTCARs.
Its INCAR is `neb/`'s with `EDIFFG = -0.02` and `LCLIMB = .TRUE.`.
```bash
cd neb/H2O_to_OH_H/cineb
cp ../../../ads/H2O/top/POTCAR .
sbatch run_neb.sh
```
Restart after the 24 h limit, or after NSW runs out, with the same loop as section 2 (`for d in 0[1-7]; ...`) and `sbatch`.
When `out` shows `reached required accuracy`, run `nebef.pl`. The highest image is the TS; verify it with a vib job (section 3).

**Result (job 50062, 100 steps, not converged):** climbing image 03 at 0.83 eV, max force 0.36. Image 07 is still the worst, at 0.73 eV/A.
`cineb/` is superseded by `cineb_bridge/` (section 5).

## 5. Split NEB: `cineb_bridge/` (IS -> OH* bridge + H*)
`ads/OH_bridge_H` (job 50064) converged in 56 steps (max force 0.019 eV/A) with **OH still on the bridge**.
Its O moved to fractional (0.322, 0.512), the bridge between top Rh atoms at (1/3, 1/3) and (1/3, 2/3).
So OH* bridge + H* fcc is a local minimum (slide 36). The full path is two steps:

| state | E0, eV | E - E(IS), eV |
|---|---|---|
| IS: H2O* top | -196.5977 | 0.00 |
| TS (`neb_IS_to_bridge` image 04, before CI) | -195.6418 | ~0.96 |
| OH* bridge + H* fcc (`ads/OH_bridge_H`) | -196.3475 | 0.25 |
| FS: OH* fcc + H* fcc (`ads/OH_H`) | -196.3728 | 0.22 |

The O-H scission (IS -> bridge) holds the TS. The bridge -> fcc OH hop is the flat tail (images 05-07 at about 0.28 eV), so its barrier is only a few hundredths of an eV.
That tail caused 07's kinked, oscillating band.

`cineb_bridge/` runs only the scission step:
- `00/`: IS (POSCAR + OUTCAR).
- `01-04/`: `cineb/`'s last CONTCARs 01-04, nearly converged. The climbing image is 03.
- `05/`: `ads/OH_bridge_H` CONTCAR + OUTCAR.
- `INCAR`: `cineb/`'s with `IMAGES = 4` (still `LCLIMB = .TRUE.`, `EDIFFG = -0.02`).
- `run_neb.sh`: 4 nodes x 16 = 64 cores, 16 per image.
```bash
cd neb/H2O_to_OH_H/cineb_bridge
cp ../../../ads/H2O/top/POTCAR .
sbatch run_neb.sh
```
Restart with `for d in 0[1-4]; do cp $d/CONTCAR $d/POSCAR; done` and `sbatch`.
When it converges, the highest image is the TS. Verify it with a vib job (exactly one imaginary mode).

**Result (job 50076, 100 steps, failed):** the band came apart. Energies above IS: 01 0.11, 02 1.00, 03 0.53, 04 0.74 eV, so the energy is no longer monotonic on either side of the peak.
Image 03 hit 3.3 eV/A at step ~80 and ended at 1.12. The reused cineb images were spaced unevenly once the endpoint changed. Abandoned in favor of `neb_IS_to_bridge/`.

## 6. Fresh 7-image NEB: `neb_IS_to_bridge/` (slide 36, done by the book)
This is the same step as `cineb_bridge/`, but with a new band built the tutorial's way instead of reusing images. Use whichever finishes cleanly; both should give the same TS.
- `POSCAR_start`/`OUTCAR_start`: IS (same as `neb/`).
- `POSCAR_end`/`OUTCAR_end`: relaxed `ads/OH_bridge_H` (E0 = -196.3475 eV).
- `INCAR`: `files/INCAR_neb`, the same as `neb/` (7 images, climbing off, EDIFFG -0.2).
- `KPOINTS` and `run_neb.sh` (112 cores) are the same as `neb/`.
```bash
cd neb/H2O_to_OH_H/neb_IS_to_bridge
cp ../../../ads/H2O/top/POTCAR .
nebinit.py -n 7                        # makes 00-08, copies OUTCARs into 00/ and 08/
ase gui 0*/POSCAR                      # check the path
sbatch run_neb.sh
```
Then follow section 2: restart with `for d in 0[1-7]; do cp $d/CONTCAR $d/POSCAR; done` until `reached required accuracy`.
Then set `EDIFFG = -0.02` and `LCLIMB = .TRUE.` and restart for the climbing image.

**Result (job 50077, regular NEB, converged at EDIFFG -0.2):** a clean band with one peak.

| image | 00 | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 |
|---|---|---|---|---|---|---|---|---|---|
| E - E(IS), eV | 0.00 | 0.02 | 0.14 | 0.74 | **0.96** | 0.79 | 0.58 | 0.33 | 0.25 |
| O-H (breaking), A | 0.98 | 0.99 | 1.01 | 1.24 | 1.82 | 2.19 | 2.50 | 2.78 | 2.99 |
| max force, eV/A | - | 0.09 | 0.16 | 0.20 | 0.20 | 0.11 | 0.10 | 0.07 | - |

Spacing is even (0.43-0.47 A), and the angles are 125-157 degrees, so there are no kinks. Image 04 will be the climbing image.
The old full-path estimate (0.83 eV) was low because its climbing image hadn't converged. 0.96 eV agrees with `cineb_bridge`'s highest image (1.00 eV).

**CI run 1 (job 50112, 100 steps, not converged):** climbing image 04 ended at 0.90 eV with max force 0.27 eV/A, rising. Its O-H went from 1.82 to 1.94 A, and the band angle at 04 is 100 degrees.
It was restarted with the slides' settings (same INCAR, CONTCAR -> POSCAR).

## 7. Second half of the split: `neb_bridge_to_fcc/` (OH* bridge + H* -> OH* fcc + H*)
The OH hops from the bridge into its fcc site, and H* stays put. The path is short (O moves about 0.9 A), so it uses 3 images.
- `POSCAR_start`/`OUTCAR_start`: relaxed `ads/OH_bridge_H` (E0 = -196.3475 eV).
- `POSCAR_end`/`OUTCAR_end`: relaxed `ads/OH_H` (E0 = -196.3728 eV), the FS of the original NEB.
- `INCAR`: `files/INCAR_neb` with `IMAGES = 3` (climbing off, EDIFFG -0.2).
- `run_neb.sh`: 3 nodes x 16 = 48 cores, 16 per image.

The OH tips over during the hop, so plain linear interpolation squeezes O-H to 0.86 A. `nebinit.py` uses IDPP by default, which avoids that. Check it in `ase gui`.
```bash
cd neb/H2O_to_OH_H/neb_bridge_to_fcc
cp ../../../ads/H2O/top/POTCAR .
nebinit.py -n 3                        # makes 00-04, copies OUTCARs into 00/ and 04/
ase gui 0*/POSCAR                      # O-H should stay near 0.98 A in every image
sbatch run_neb.sh
```
Restart with `for d in 0[1-3]; do cp $d/CONTCAR $d/POSCAR; done`.
After `reached required accuracy`, set `EDIFFG = -0.02` and `LCLIMB = .TRUE.` and restart for the climbing image.
The barrier is expected to be only a few hundredths of an eV, so it doesn't control the rate. It completes the energy diagram from H2O* to OH* fcc + H*.

**Results:** the regular NEB (job 50078) converged in 2 steps. CI (job 50079) ran 100 steps without reaching 0.02. The final max forces are 01 0.032, 02 0.061 and 03 0.055 eV/A, still falling slowly.
Energies above the bridge state: 01 0.03, 02 (climbing) **0.09**, 03 0.03 eV; FS -0.03 eV.
After two more restarts (jobs 50114, 50135), CI converged: 02 = **0.088 eV**, all forces < 0.02. TS vib job `vib/TS_OH_bridge_fcc` (job 50162): exactly one imaginary mode, 264i cm-1, so image 02 is a true first-order saddle point.

## 8. Two CI restarts of `neb_IS_to_bridge/`: slides' optimizer vs FIRE
In job 50112, image 04's real force along the path fell from 0.28 to 0.12 eV/A by step 50, then rose to 0.20 by step 100.
It overshot the peak, and the slides' damped MD (`IBRION = 3`, `POTIM = 0.01`) is slow to turn back. Both folders start from 50112's CONTCARs (01-07). 00/ and 08/ have their POSCAR and OUTCAR.

| folder | optimizer | INCAR change |
|---|---|---|
| `ci_IS_to_bridge_damped/` | damped MD (the slides) | none |
| `ci_IS_to_bridge_fire/` | FIRE (VTST) | `IOPT = 7`, `POTIM = 0` (`IBRION = 3` stays) |

Both use `EDIFFG = -0.02`, `LCLIMB = .TRUE.` and 7 nodes (112 cores). Submit both:
```bash
cd neb/H2O_to_OH_H/split/ci_IS_to_bridge_damped
cp ../../../ads/H2O/top/POTCAR .
sbatch run_neb.sh
cd ../ci_IS_to_bridge_fire
cp ../../../ads/H2O/top/POTCAR .
sbatch run_neb.sh
```
Restart either one with `for d in 0[1-7]; do cp $d/CONTCAR $d/POSCAR; done` and `sbatch`.
Track image 04's force along the path: `grep "projections on to tangent" 04/OUTCAR | awk 'NR%10==0{print NR,$NF}'`.
Both should reach the same TS. Report the slides' result if it converges; FIRE is the fallback.

**First attempt (jobs 50163, 50164): crashed after 1 step.** The folders were moved into `split/` while the jobs ran,
so VASP couldn't find its own output files (`forrtl: severe (29): file not found ... XDATCAR`). It was not a convergence problem.
`NSW` is now 300 in both folders. Resubmit from `split/` without restaging; the POSCARs are still the intended start.
Don't move or rename a folder while its job is running. At about 3.5 min per step, 300 steps take about 18 h, which fits within the 24 h limit.

**Second attempt (jobs 50181 damped, 50182 FIRE; ran with NSW 100 because the force-push dropped the 300).**

| image | 01 | 02 | 03 | 04 | 05 | 06 | 07 |
|---|---|---|---|---|---|---|---|
| damped: E - E(IS), eV | 0.02 | 0.16 | **0.96** | 0.69 | 0.77 | 0.74 | 0.45 |
| damped: max force | 0.25 | 0.69 | 0.73 | 0.93 | 0.24 | 0.32 | 0.06 |
| FIRE: E - E(IS), eV | 0.02 | 0.08 | 0.29 | **0.94** | 0.43 | 0.69 | 0.69 |
| FIRE: max force | 0.15 | 0.06 | 0.20 | **0.07** | 0.11 | 0.23 | 0.44 |
| FIRE: breaking O-H, A | 0.98 | 0.99 | 1.05 | 1.56 | 3.31 | 3.35 | 2.85 |

- **damped:** the band fell apart. Forces grew in every middle image, the peak jumped from 04 to 03, the energy is no longer smooth (04 < 05 > 06), and the angle at 04 is 81 degrees. Don't continue it.
- **FIRE:** the climbing image 04 is almost a stationary point. Its max force is 0.066 eV/A, its force along the path is 0.01-0.03, it sits at 0.94 eV, and the breaking O-H is 1.56 A.
  The rest of the band is tangled: in 05-06 the H went past its final site (O-H 3.3 A vs 2.99 in the FS), and 06-07 sit at 0.69 eV, above the FS (0.25).
- TS check: `vib/TS_H2O_to_OH_bridge` is a vib job on FIRE's image 04 CONTCAR. Exactly one imaginary mode (an O-H stretch) confirms the TS at about 0.94 eV, even with the rest of the band untidy.
