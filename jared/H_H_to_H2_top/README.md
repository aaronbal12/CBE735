# NEB: H* + H* -> H2* (top), Rh(111) 3x3

Corrected setup for `../H_H_NEB_Calculation_Climb`, which didn't converge. That folder is left as it was.

## Why the old run didn't converge
1. **The final state wasn't the one intended.** 9 images were generated (folders 00-10), but the INCAR has `IMAGES = 7`.
   VASP therefore used `08/` as the final state: an unrelaxed interpolated image with H2 stretched to 1.14 A, 3.5 A above the surface.
   The real FS (`10/`) was ignored. 00/ and 08/ also had no OUTCARs.
2. **There was no barrier for CI to climb to.** The FS had H2 desorbed 4 A into the vacuum.
   The energy rose steadily, then went flat once H2 left the surface (images 05-07 within 0.015 eV), so the climbing image had no maximum to find.
3. **Images were too far apart and the band folded.** Images were 1.8 A apart (about 0.7 A is typical), and the band bent back on itself at image 02 (24 degrees).
   The H atoms in 01-03 had wandered to other hollow sites. 01 and 02 were stuck at 0.44-0.54 eV/A, mostly spring force.
4. The IS had one H fixed (`F F F`), so the IS wasn't fully relaxed. `LCHARG = .FAST.` isn't valid; VASP reads it as `.FALSE.`.

## What changed
Image 04 of the old band was H2* on top: E = -189.1786 eV, H-H 0.92 A, 1.6 A up. That matches the converged `ads/H2/top` (-189.1795 eV).

**Caveat: H2* top is at best a very shallow minimum.** It only holds in the symmetric geometry (H2 centred on the Rh atom, both H at the same height).
Started off-symmetry, H2 tends to split or desorb. `vib/H2_top` shows no imaginary mode, but the softest mode is 65 cm-1 (8 meV).
In the old band, image 04 sat on the rising slope (03 -> 04 -> 05 all uphill), not in a dip.
So treat it as a convenient end point, not a real intermediate. Expect the NEB below to rise steadily with at most a tiny maximum near the end.
The step splits in two:

| step | from -> to | how |
|---|---|---|
| A | 2H* fcc -> H2* top | **NEB** (this folder). Contains the barrier for H-H bond formation |
| B | H2* top -> H2(g) + * | no NEB: desorption with no barrier beyond its own energy. Use dE = E(H2 gas) + E(clean) - E(H2* top) |

Reference energies (E0, eV): clean Rh -182.0572 (`clean/Rh`), H2 gas -6.7595 (`gas/H2/relax`), H2* top -189.1795 (`ads/H2/top`), H* fcc -185.8832 (`ads/H/fcc_hollow`).
- Step B: dE = -6.7595 - 182.0572 + 189.1795 = **+0.36 eV**.
- The old band's plateau (-188.81) equals H2(g) + clean slab (-188.82). That's consistent: those images were just H2 in the vacuum.

**IS:** the same sites as before, both H in fcc hollows next to the top Rh at (1/3, 1/3): H1 at (1/9, 4/9) and H2 at (4/9, 4/9), H-H 2.69 A.
H2* sits on exactly that Rh atom, so each H moves only about 1.3 A.
**FS:** `ads/H2/top` CONTCAR (H2 flat on top, H-H 0.94 A). The atom order matches the IS: Rh x27, then the H nearer (1/9, 4/9), then the other H.

## 1. Relax the IS (`IS_relax/`, 1 node)
`POSCAR` is the old `POSCAR_start` with the fixed H freed (`T T T`). `INCAR` and `KPOINTS` are the same as `ads/H2/top`, so the energies compare.
```bash
cd jared/H_H_to_H2_top/IS_relax
cp ../../../ads/H2/top/POTCAR .        # Rh H order
sbatch run.sh
```
Check that `out` says `reached required accuracy`.

## 2. NEB (`neb/`, 5 images, 5 nodes x 16 = 80 cores)
`INCAR` = `files/INCAR_neb` with `IMAGES = 5` (climbing off, EDIFFG -0.2).
The path is only about 1.85 A long, so 5 images give about 0.3 A spacing. `POSCAR_end`/`OUTCAR_end` = `ads/H2/top`.
```bash
cd ../neb
cp ../IS_relax/CONTCAR POSCAR_start
cp ../IS_relax/OUTCAR OUTCAR_start
cp ../../../ads/H2/top/POTCAR .
nebinit.py -n 5                        # makes 00-06, copies OUTCARs into 00/ and 06/
ls 0*                                  # must be exactly 00-06; IMAGES = 5 needs 5 middle folders
ase gui 0*/POSCAR                      # H should move from the fcc sites toward the top Rh
sbatch run_neb.sh
```
Restart until `out` says `reached required accuracy`:
```bash
for d in 0[1-5]; do cp $d/CONTCAR $d/POSCAR; done
sbatch run_neb.sh
```
Then switch on the climbing image in the same folder and restart:
```bash
sed -i 's/^EDIFFG.*/EDIFFG   = -0.02/; s/^LCLIMB.*/LCLIMB   = .TRUE./' INCAR
for d in 0[1-5]; do cp $d/CONTCAR $d/POSCAR; done
sbatch run_neb.sh
```
Check with `nebef.pl`. Column 2 is max force, which must be below the EDIFFG value for images 01-05. Column 4 is energy relative to 00.
Verify the highest image with a vib job: exactly one imaginary mode.

**If the energy only rises from 2H* to H2* with no maximum** (likely, given the caveat above), H2 formation has no barrier beyond its reaction energy.
The barrier for the whole H* + H* -> H2(g) step is then just its reaction energy: E(H2 gas) + E(clean) - E(IS relaxed).
That's about +0.8 eV, judging by the old band (01-03 at about -189.62, plateau at -188.81).
If that's acceptable for the project, ask the instructor; the IS relaxation alone then gives the number, and the NEB is just a check.
