"""Set up frequency, H2 and gas-phase jobs from converged relaxations.

  vib/<ads>_<site>/   IBRION=5 on the relaxed CONTCAR; slab fixed, adsorbate free
  ads/H2/top/         H2 molecule lying flat over a top Rh (relaxation)
  gas/<mol>/relax/    isolated H2, CO, CO2 in a 15 A cubic box (Gamma point)
  gas/<mol>/vib/      frequencies of the isolated molecule (all atoms free)

INCARs come from files/INCAR_relax and files/INCAR_vib with ALGO = Normal.
"""
import os
import re
import shutil

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
N_SLAB = 27
BOX = 15.0

# vib dir name -> converged relaxation whose CONTCAR is the starting geometry
VIB = {
    "H_fcc": "ads/H/fcc_hollow",        # H started on bridge relaxes into fcc
    "O_fcc": "ads/O/fcc_hollow",
    "OH_fcc": "ads/OH/fcc_hollow",
    "H2O_top": "ads/H2O/top",
    "CO_hcp": "ads/CO/hcp_hollow",
}

# Isolated molecules: PBE equilibrium bond lengths (A), centred in the box
GAS = {
    "H2": [("H", (0, 0, -0.3755)), ("H", (0, 0, 0.3755))],
    "CO": [("C", (0, 0, -0.5715)), ("O", (0, 0, 0.5715))],
    "CO2": [("C", (0, 0, 0)), ("O", (0, 0, -1.176)), ("O", (0, 0, 1.176))],
}


def incar(template, **changes):
    """Copy a files/ INCAR template, replacing (or dropping, if None) tags."""
    out = []
    for line in open(os.path.join(ROOT, "files", template)).read().splitlines():
        m = re.match(r"\s*([A-Z]+)\s*=", line)
        if m and m.group(1) in changes:
            val = changes.pop(m.group(1))
            if val is None:
                continue
            comment = line[line.index("#"):] if "#" in line else ""
            line = f"{m.group(1):<8} = {val} {comment}".rstrip()
        out.append(line)
    assert not changes, f"tags not in {template}: {changes}"
    return "\n".join(out) + "\n"


def read_poscar(path):
    L = open(path).read().splitlines()
    cell = [[float(x) for x in l.split()[:3]] for l in L[2:5]]
    els, cnt = L[5].split(), [int(x) for x in L[6].split()]
    i = 8 if L[7].strip()[0] in "sS" else 7
    direct = L[i].strip()[0] in "dD"
    xyz = []
    for l in L[i + 1:i + 1 + sum(cnt)]:
        f = [float(x) for x in l.split()[:3]]
        if direct:
            f = [sum(f[k] * cell[k][j] for k in range(3)) for j in range(3)]
        xyz.append(f)
    return cell, els, cnt, xyz


def write_poscar(path, title, cell, els, cnt, xyz, free):
    with open(path, "w") as f:
        f.write(f"{title}\n 1.0000000000000000\n")
        for v in cell:
            f.write("  " + "  ".join(f"{x:20.16f}" for x in v) + "\n")
        f.write("  " + "  ".join(els) + "\n  " + "  ".join(map(str, cnt)) + "\n")
        f.write("Selective dynamics\nCartesian\n")
        for p, fr in zip(xyz, free):
            flag = "T   T   T" if fr else "F   F   F"
            f.write("  " + "  ".join(f"{x:20.16f}" for x in p) + f"   {flag}\n")


def job_dir(rel, incar_text, kpoints):
    d = os.path.join(ROOT, rel)
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "INCAR"), "w").write(incar_text)
    open(os.path.join(d, "KPOINTS"), "w").write(kpoints)
    shutil.copy(os.path.join(ROOT, "ads", "CO", "bridge", "run.sh"), d)
    return d


def main():
    slab_kp = open(os.path.join(ROOT, "files", "KPOINTS")).read()
    gamma_kp = "Gamma point only\n0\nGamma\n1 1 1\n0 0 0\n"
    vib_slab = incar("INCAR_vib", ALGO="Normal")
    relax_slab = incar("INCAR_relax", ALGO="Normal", ISMEAR="0")
    # Isolated molecules: no dipole correction, small smearing
    gas_mods = dict(ALGO="Normal", SIGMA="0.01", IDIPOL=None, LDIPOL=None)
    relax_gas = incar("INCAR_relax", ISMEAR="0", **gas_mods)
    vib_gas = incar("INCAR_vib", **gas_mods)

    # 1) adsorbate frequencies from converged CONTCARs
    for name, src in VIB.items():
        cell, els, cnt, xyz = read_poscar(os.path.join(ROOT, src, "CONTCAR"))
        d = job_dir(f"vib/{name}", vib_slab, slab_kp)
        free = [i >= N_SLAB for i in range(len(xyz))]
        write_poscar(os.path.join(d, "POSCAR"), f"Rh(111) 3x3 + {name} vib (from {src})",
                     cell, els, cnt, xyz, free)

    # 2) H2 "adsorption": molecule flat over a top Rh, 2.0 A above it (H-Rh ~2.03 A)
    cell, els, cnt, xyz = read_poscar(os.path.join(ROOT, "clean", "Rh", "CONTCAR"))
    z_top = max(p[2] for p in xyz)
    x0, y0 = 4.0305086527633209, 2.3270152556440182
    h2 = [[x0 - 0.3755, y0, z_top + 2.0], [x0 + 0.3755, y0, z_top + 2.0]]
    d = job_dir("ads/H2/top", relax_slab, slab_kp)
    write_poscar(os.path.join(d, "POSCAR"), "Rh(111) 3x3 + H2 top (flat)", cell,
                 ["Rh", "H"], [N_SLAB, 2], xyz + h2, [False] * N_SLAB + [True] * 2)

    # 3) isolated molecules in a box: relaxation + frequencies
    box = [[BOX, 0, 0], [0, BOX, 0], [0, 0, BOX]]
    for mol, atoms in GAS.items():
        els = list(dict.fromkeys(a[0] for a in atoms))
        cnt = [sum(a[0] == e for a in atoms) for e in els]
        xyz = [[BOX / 2 + c for c in a[1]] for e in els for a in atoms if a[0] == e]
        for job, text in (("relax", relax_gas), ("vib", vib_gas)):
            d = job_dir(f"gas/{mol}/{job}", text, gamma_kp)
            write_poscar(os.path.join(d, "POSCAR"), f"{mol} gas {job}", box, els, cnt,
                         xyz, [True] * len(xyz))


if __name__ == "__main__":
    main()
