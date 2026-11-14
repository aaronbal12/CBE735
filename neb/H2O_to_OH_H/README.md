# NEB: H2O* + * -> OH* + H*  (Rh(111) 3x3, WGS step 1)

Follows VASP Tutorial 2, slides 25-36.

| dir  | what it is |
|------|------------|
| `IS/`  | H2O* on top: copy of `ads/H2O/top`. Lowest energy (E0 = -196.598 eV) and the only H2O site that converged |
| `neb/` | `INCAR` = `files/INCAR_neb` (7 images, climbing off, EDIFFG -0.2) + `KPOINTS` |

**Final state:** OH* fcc + H* fcc co-adsorbed. It's relaxed as its own thermo job in `ads/OH_H`, so run that first.
Its atom order matches IS (Rh x27, O, H that stays on O, H that leaves), which `nebmake.pl` needs.

## 0. Script fixes (once, from the repo root)
`scripts/` holds copies of the `~/bin` scripts. Two were fixed here; copy them back over the `~/bin` versions:
* `nebmake.pl`: two input checks never fired because of typos (`$totatom1` instead of `$totatoms1`, and the cell compared with itself).
* `nebinit.py`: its first line was `#!/.conda/...`, which doesn't exist. It now points to the same Python as `kie.py` and `order_poscar_elements.py`,
  and it imports NEB from where newer ASE versions keep it.
```bash
cp scripts/nebmake.pl scripts/nebinit.py ~/bin/
ls ~/bin/Vasp.pm      # nebmake.pl needs this helper next to it; it wasn't in the push
```
`nebef.pl` and `run_neb.sh` work as-is. `nebef.pl`'s column cuts were checked against these VASP 6.4 OUTCARs.
`run_neb.sh` uses 112 cores = 16 per image, which NCORE = 4 divides evenly.

## 1. Set up and submit the NEB (slides 26-31), after `ads/OH_H` converges
```bash
cd neb/H2O_to_OH_H/neb
cp ../IS/CONTCAR POSCAR_start
cp ../../../ads/OH_H/CONTCAR POSCAR_end
sed -i '1s/.*/Rh O H/' POSCAR_start POSCAR_end   # nebmake copies line 1 into every image; ag can read the elements from it
cp ../../../ads/H2O/top/POTCAR .
nebmake.pl POSCAR_start POSCAR_end 7
cp ../IS/OUTCAR 00/                    # nebef.pl reads the endpoint energies from these
cp ../../../ads/OH_H/OUTCAR 08/
ag 0*/POSCAR                           # check the path
cp ../../../scripts/run_neb.sh .
sbatch run_neb.sh
```
* Linear interpolation squeezes the O-H bond that stays to about 0.74 A in images 03-05 (normal is 0.98 A).
  No atoms overlap, so it should run. If it misbehaves, rebuild the images with `nebinit.py` instead of `nebmake.pl`.
  It reads `POSCAR_start`/`POSCAR_end` in the same folder and does a smoother (IDPP) interpolation.

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
