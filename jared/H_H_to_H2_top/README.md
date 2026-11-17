# NEB: H* + H* -> H2* (top), Rh(111) 3x3

Corrected setup for `../H_H_NEB_Calculation_Climb`, which didn't converge. That folder is left as it was.

## Why the old run didn't converge
1. **The final state wasn't the one intended.** 9 images were generated (folders 00-10), but the INCAR has `IMAGES = 7`.
   VASP therefore used `08/` as the final state: an unrelaxed interpolated image with H2 stretched to 1.14 A, 3.5 A above the surface.
   The real FS (`10/`) was ignored. 00/ and 08/ also had no OUTCARs.
2. **There was no barrier for CI to climb to.** The intended FS (`POSCAR_end` = `ads/H2/fcc_hollow` CONTCAR) is H2 that desorbed to 4 A above the surface.
   The energy rose steadily, then went flat once H2 left the surface (images 05-07 within 0.015 eV), so the climbing image had no maximum to find.
   That plateau (-188.81 eV) equals H2 gas + clean slab (-188.82).
3. **Images were too far apart and the band folded.** Images were 1.8 A apart (about 0.7 A is typical), and the band bent back on itself at image 02 (24 degrees).
   The H atoms in 01-03 had wandered to other hollow sites. 01 and 02 were stuck at 0.44-0.54 eV/A, mostly spring force.
4. `LCHARG = .FAST.` isn't valid; VASP reads it as `.FALSE.`. Harmless.

## The new band
- **IS (`POSCAR_start`):** the same relaxed 2H* structure as before. Both H are in fcc hollows next to the top Rh at (1/3, 1/3), H-H 2.69 A.
  The only change: the first H's flags are now `T T T` instead of `F F F`. The geometry is unchanged.
  With `F F F`, the constraint can carry over into the images, and that H could never move.
- **FS (`POSCAR_end`, `OUTCAR_end`):** `ads/H2/top` (H2 flat on top of that same Rh, H-H 0.94 A, E0 = -189.1795 eV). The atom order matches the IS.
  - It's a very shallow minimum held mostly by symmetry. Off-symmetry, H2 splits or desorbs. `vib/H2_top` shows no imaginary mode, but its softest mode is 65 cm-1.
  - It's still a valid NEB endpoint, because VASP never moves endpoints.
  - It keeps the band on the surface, which the gas-phase FS didn't.
- **7 images** (the course default): the path is about 1.85 A, so the spacing is about 0.23 A. `INCAR` = `files/INCAR_neb` (`IMAGES = 7`, climbing off, EDIFFG -0.2).
- **run_neb.sh:** 7 nodes x 16 = 112 cores.

## Run
`OUTCAR_start` isn't in the repo. Copy it from wherever the IS relaxation ran.
```bash
cd jared/H_H_to_H2_top/neb
cp <IS relaxation folder>/OUTCAR OUTCAR_start
cp ../../../ads/H2/top/POTCAR .        # Rh H order
nebinit.py -n 7                        # makes 00-08, copies OUTCARs into 00/ and 08/
ls -d 0*                               # must be exactly 00-08; IMAGES = 7 needs 7 middle folders
ase gui 0*/POSCAR                      # both H should move from the fcc sites toward the top Rh
sbatch run_neb.sh
```
Restart until `out` says `reached required accuracy`:
```bash
for d in 0[1-7]; do cp $d/CONTCAR $d/POSCAR; done
sbatch run_neb.sh
```
Then switch on the climbing image in the same folder and restart:
```bash
sed -i 's/^EDIFFG.*/EDIFFG   = -0.02/; s/^LCLIMB.*/LCLIMB   = .TRUE./' INCAR
for d in 0[1-7]; do cp $d/CONTCAR $d/POSCAR; done
sbatch run_neb.sh
```
Check with `nebef.pl`: column 2 (max force) must be below the EDIFFG value for images 01-07, and column 4 is energy relative to 00.

## What to expect
- **Most likely:** a steady rise from 2H* to H2* top (about +0.44 eV, from the old band), with no peak or only a tiny one near the end.
  Then H-H formation has no barrier beyond its reaction energy, and the climbing image ends up next to the FS.
  The H2* -> H2(g) desorption that follows adds +0.36 eV, so the barrier for the whole H* + H* -> H2(g) step is about +0.8 eV, equal to its reaction energy.
- **If there's a real maximum:** it's the TS. Verify it with a vib job (exactly one imaginary mode).

Reference E0 (eV): clean Rh -182.0572 (`clean/Rh`), H2 gas -6.7595 (`gas/H2/relax`), H2* top -189.1795 (`ads/H2/top`), H* fcc -185.8832 (`ads/H/fcc_hollow`).
