"""Build adsorbate POSCARs on the 3x3 Rh(111) slab (clean/Rh/POSCAR_Rh).

Sites (Cartesian x, y) are chosen near the cell centre on the top layer
(z = 4.388 A). Stacking of the slab is ...C B A (top):
  top    : above top-layer Rh at (4.0305, 2.3270)
  bridge : midpoint of two top-layer Rh, (4.0305, 2.3270)-(5.3740, 4.6540)
  hcp    : above a 2nd-layer Rh,   (4.0305, 3.8784)
  fcc    : above a 3rd-layer Rh,   (5.3740, 3.1027)
Heights come from target metal-adsorbate bond lengths (typical DFT-PBE values
for Rh(111)). The slab keeps its constraints; adsorbate atoms are free (T T T).
"""
import math
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
A_NN = 2.6870057685088806          # in-plane Rh-Rh distance
R_BRIDGE = A_NN / 2                # lateral distance bridge site -> Rh
R_HOLLOW = A_NN / math.sqrt(3)     # lateral distance hollow site -> Rh

SITES = {
    "top": (4.0305086527633209, 2.3270152556440182, 0.0),
    "bridge": (4.7022615, 3.4905229, R_BRIDGE),
    "hcp_hollow": (4.0305086527633209, 3.8783587594066966, R_HOLLOW),
    "fcc_hollow": (5.3740115370177612, 3.1026870075253576, R_HOLLOW),
}


def height(bond, lateral):
    """Height above the Rh plane giving the requested Rh-X bond length."""
    return math.sqrt(bond ** 2 - lateral ** 2)


def site_kind(site):
    return "top" if site == "top" else ("bridge" if site == "bridge" else "hollow")


# Target Rh-X bond lengths (A) per site kind for the surface-bound atom.
BONDS = {
    "H": {"top": 1.57, "bridge": 1.75, "hollow": 1.86},        # Rh-H
    "O": {"top": 1.76, "bridge": 1.95, "hollow": 2.03},        # Rh-O
    "OH": {"top": 2.03, "bridge": 2.12, "hollow": 2.20},       # Rh-O(H)
    "CO": {"top": 1.83, "bridge": 2.02, "hollow": 2.10},       # Rh-C
    "H2O": {"top": 2.30, "bridge": 2.60, "hollow": 2.60},      # Rh-O (flat water)
}
D_CO = {"top": 1.160, "bridge": 1.180, "hollow": 1.190}        # C-O
D_OH = 0.975                                                  # O-H
HOH = math.radians(104.5)
D_CO2 = 1.176                                                 # C=O in CO2
H_CO2 = 3.10                                                  # physisorption height


def adsorbate(name, site):
    """Return list of (element, (dx, dy, dz)) relative to (site_x, site_y, z_top)."""
    kind = site_kind(site)
    lat = SITES[site][2]
    if name in ("H", "O"):
        return [(name, (0.0, 0.0, height(BONDS[name][kind], lat)))]
    if name == "CO":
        h = height(BONDS["CO"][kind], lat)
        return [("C", (0.0, 0.0, h)), ("O", (0.0, 0.0, h + D_CO[kind]))]
    if name == "OH":
        h = height(BONDS["OH"][kind], lat)
        if kind == "top":
            # On-top OH tilts: Rh-O-H ~ 105 deg (H tipped toward a hollow, +y)
            tilt = math.radians(180 - 105)
            return [("O", (0.0, 0.0, h)),
                    ("H", (0.0, D_OH * math.sin(tilt), h + D_OH * math.cos(tilt)))]
        return [("O", (0.0, 0.0, h)), ("H", (0.0, 0.0, h + D_OH))]
    if name == "H2O":
        # Near-flat water, O over site, molecular plane tilted ~10 deg so the
        # H atoms sit slightly above O; H atoms point along +x away from site.
        h = height(BONDS["H2O"][kind], lat)
        tilt = math.radians(10)
        half = HOH / 2
        along = D_OH * math.cos(half)          # along bisector
        perp = D_OH * math.sin(half)           # across bisector (in plane)
        dx, dz = along * math.cos(tilt), along * math.sin(tilt)
        return [("O", (0.0, 0.0, h)),
                ("H", (dx, +perp, h + dz)),
                ("H", (dx, -perp, h + dz))]
    if name == "CO2":
        # Linear CO2 physisorbed flat (O=C=O parallel to surface), C over site.
        return [("C", (0.0, 0.0, H_CO2)),
                ("O", (+D_CO2, 0.0, H_CO2)),
                ("O", (-D_CO2, 0.0, H_CO2))]
    raise ValueError(name)


def read_slab(path):
    with open(path) as f:
        lines = f.read().splitlines()
    cell = lines[2:5]
    n = int(lines[6].split()[0])
    atoms = lines[9:9 + n]
    return cell, atoms


def write(path, title, cell, slab, ads):
    order = []
    for el, _ in ads:
        if el not in order:
            order.append(el)
    # Rh first, then adsorbate species grouped in order of appearance
    grouped = [a for el in order for a in ads if a[0] == el]
    counts = [len(slab)] + [sum(1 for a in ads if a[0] == el) for el in order]
    with open(path, "w") as f:
        f.write(title + "\n 1.0000000000000000\n")
        f.write("\n".join(cell) + "\n")
        f.write(" " + "  ".join(["Rh"] + order) + "\n")
        f.write(" " + "  ".join(str(c) for c in counts) + "\n")
        f.write("Selective dynamics\nCartesian\n")
        f.write("\n".join(slab) + "\n")
        for _, (x, y, z) in grouped:
            f.write(f"  {x:.16f}  {y:.16f}  {z:.16f}   T   T   T\n")


def main():
    cell, slab = read_slab(os.path.join(ROOT, "clean", "Rh", "POSCAR_Rh"))
    z_top = max(float(l.split()[2]) for l in slab)
    for name in ("H", "O", "OH", "H2O", "CO", "CO2"):
        for site, (sx, sy, _) in SITES.items():
            ads = [(el, (sx + dx, sy + dy, z_top + dz))
                   for el, (dx, dy, dz) in adsorbate(name, site)]
            path = os.path.join(ROOT, "ads", name, site, "POSCAR")
            write(path, f"Rh(111) 3x3 + {name} {site}", cell, slab, ads)
            print(path)


if __name__ == "__main__":
    main()
