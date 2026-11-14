# NEB: H2O* + * -> OH* + H*  (Rh(111) 3x3, WGS step 1)

Follows VASP Tutorial 2, slides 25-36.

| dir  | what it is |
|------|------------|
| `IS/`  | H2O* on top: copy of `ads/H2O/top`. Lowest energy (E0 = -196.598 eV) and the only H2O site that converged |
| `neb/` | `INCAR` = `files/INCAR_neb` (7 images, climbing off, EDIFFG -0.2) + `KPOINTS` |

**Final state:** OH* fcc + H* fcc co-adsorbed. It's relaxed as its own thermo job in `ads/OH_H`, so run that first.
Its atom order matches IS (Rh x27, O, H that stays on O, H that leaves), which `nebmake.pl` needs.

## 1. Set up and submit the NEB (slides 26-31), after `ads/OH_H` converges
```bash
cd neb/H2O_to_OH_H/neb
cp ../IS/CONTCAR POSCAR_start
cp ../../../ads/OH_H/CONTCAR POSCAR_end
cp ../../../ads/H2O/top/POTCAR .
nebmake.pl POSCAR_start POSCAR_end 7
cp ../IS/OUTCAR 00/                    # nebef.pl reads the endpoint energies from these
cp ../../../ads/OH_H/OUTCAR 08/
ag 0*/POSCAR                           # check the path
cp ~/CBE735/files/run_neb.sh .         # slide 31 says ~/files/, but the file is in ~/CBE735/files
sbatch run_neb.sh
```
* Linear interpolation squeezes the O-H bond that stays to about 0.74 A in images 03-05 (normal is 0.98 A).
  No atoms overlap, so it should run. If it misbehaves, rebuild the images with `nebinit.py` (slide 27).
* VASP needs the total core count in `run_neb.sh` to be a multiple of IMAGES = 7, and the cores per image to be divisible by NCORE = 4,
  so use 28 or 56 cores. Check with `grep -E "ntasks|-n " run_neb.sh`.

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
