# NEB: H2O* + * -> OH* + H*  (Rh(111) 3x3, WGS step 1)

Follows VASP Tutorial 2, slides 25-36.

| dir  | what it is |
|------|------------|
| `IS/`  | H2O* on top: copy of `ads/H2O/top`. Lowest energy (E0 = -196.598 eV) and the only H2O site that converged |
| `neb/` | `INCAR` = `files/INCAR_neb` (7 images, climbing off, EDIFFG -0.2) + `KPOINTS` |

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

## 1. Set up and submit the NEB (slides 26-31), after `ads/OH_H` converges
```bash
cd neb/H2O_to_OH_H/neb
cp ../IS/CONTCAR POSCAR_start
cp ../../../ads/OH_H/CONTCAR POSCAR_end
cp ../IS/OUTCAR OUTCAR_start           # nebinit.py copies these into 00/ and 08/ for nebef.pl
cp ../../../ads/OH_H/OUTCAR OUTCAR_end
cp ../../../ads/H2O/top/POTCAR .
nebinit.py -n 7                        # makes 00-08
ag 0*/POSCAR                           # check the path
cp ../../../scripts/run_neb.sh .
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
