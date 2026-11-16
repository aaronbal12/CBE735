# NEB: H2O* + * -> OH* + H*  (Rh(111) 3x3, WGS step 1)

Follows VASP Tutorial 2, slides 25-36.

| dir  | what it is |
|------|------------|
| `IS/`  | H2O* on top: copy of `ads/H2O/top`. Lowest energy (E0 = -196.598 eV) and the only H2O site that converged |
| `neb/` | `INCAR` = `files/INCAR_neb` (7 images, climbing off, EDIFFG -0.2) + `KPOINTS` |
| `cineb/` | climbing-image NEB over the full path (superseded) |
| `cineb_bridge/` | climbing-image NEB, IS -> OH* bridge + H* (the step with the TS) |
| `neb_IS_to_bridge/` | same step as `cineb_bridge/`, done the slides' way: fresh 7-image `nebinit.py` band, regular NEB then CI |

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
| TS (CI image 03, not converged yet) | -195.7699 | ~0.83 |
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
