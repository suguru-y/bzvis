#!/usr/bin/env python3
"""
bzvis.py: Creating a CIF file to visualize Brillouin zone with VESTA
=============================================
@author: Suguru Yoshida
Cleaned up and improved, thanks to ChatGPT's suggestions.
version 0.9 (July 5, 2025)

How to use:
    $ python bzvis. -i POSCAR -o myBZ
    or
    $ python bzvis. -i hoge.cif -o myBZ
    then,
    myBZ.cif will be created.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Tuple

import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from scipy.spatial import ConvexHull, Voronoi

# -----------------------------------------------------------------------------
# Command Line Interface
# -----------------------------------------------------------------------------

def _get_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate a CIF with the first Brillouin zone (POSCAR, CIF, …)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", default="POSCAR", help="Structure file (POSCAR, CIF, …)")
    p.add_argument("-o", "--output", default="bz", help="Output basename (writes <name>.cif)")
    p.add_argument("-s", "--symmetrize", action="store_true", help="Reduce to primitive cell by symmetry")
    p.add_argument("-n", "--order", type=int, default=1, metavar="N", help="Neighbourhood order (1→3³, 2→5³, …)")
    return p

# -----------------------------------------------------------------------------
# Centring tables
# -----------------------------------------------------------------------------

_REAL2RECIP = {
    "P": "P", "I": "F", "F": "I", "A": "A", "B": "B", "C": "C", "R": "R", "H": "H",
}
_CENTERING_OFFSETS = {
    "P": [],
    "I": [[0.5, 0.5, 0.5]],
    "F": [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]],
    "A": [[0.0, 0.5, 0.5]],
    "B": [[0.5, 0.0, 0.5]],
    "C": [[0.5, 0.5, 0.0]],
    "R": [],
    "H": [],
}
HALF_VECTORS_ABC = {
    "A": np.array([0.0, 0.5, 0.5]),
    "B": np.array([0.5, 0.0, 0.5]),
    "C": np.array([0.5, 0.5, 0.0]),
}
PIVOT_FRAC = np.array([0.5, 0.5, 0.5])
TOL = 1e-3

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def _centres_structure(struct: Structure, vec: np.ndarray, tol: float = TOL) -> bool:
    f = struct.frac_coords % 1.0
    shifted = (f + vec) % 1.0
    sorter = lambda a: a[np.lexsort(a.T)]
    return np.allclose(sorter(f), sorter(shifted), atol=tol)


def _make_sampling_points(order: int, extra: List[List[float]]) -> List[Tuple[float, float, float]]:
    rng = range(-order, order + 1)
    pts = [(i, j, k) for i in rng for j in rng for k in rng]
    for off in extra:
        pts.extend([(i + off[0], j + off[1], k + off[2]) for i in rng for j in rng for k in rng])
    return pts

# -----------------------------------------------------------------------------
# Main routine
# -----------------------------------------------------------------------------

def main(ns: argparse.Namespace) -> None:
    path = Path(ns.input).expanduser().resolve()
    if not path.exists():
        sys.exit(f"[bzvis] Input not found: {path}")

    # 1. Get Real-Space Structure ---------------------------------------------
    struct = Structure.from_file(str(path))
    if ns.symmetrize:
        struct = struct.get_primitive_structure()

    # 2. Get Real-Space Centering ---------------------------------------------
    try:
        real_letter = SpacegroupAnalyzer(struct, symprec=1e-3, angle_tolerance=5).get_space_group_symbol()[0]
    except Exception:
        real_letter = "P"

    if real_letter == "C":
        for letter, vec in HALF_VECTORS_ABC.items():
            if _centres_structure(struct, vec):
                real_letter = letter
                break

    recip_letter = _REAL2RECIP.get(real_letter, "P")
    extra = _CENTERING_OFFSETS[recip_letter]

    # 3. Create Reciprocal lattice -------------------------------------------
    lat_R = struct.lattice.reciprocal_lattice
    mat_R = lat_R.matrix
    mat_R_inv = np.linalg.inv(mat_R)

    # 4. k-point sampling ----------------------------------------------------
    frac_pts = _make_sampling_points(ns.order, extra)
    cart_pts = np.dot(np.asarray(frac_pts), mat_R)

    # 5. Create Voronoi cell -------------------------------------------------
    vor = Voronoi(cart_pts)
    centre_idx = frac_pts.index((0.0, 0.0, 0.0))
    verts_cart = np.array([vor.vertices[i] for i in vor.regions[vor.point_region[centre_idx]]])
    verts_cart = ConvexHull(verts_cart).points

    # 6. Shift so (0 0 0) => (1/2 1/2 1/2) -----------------------------------
    verts_frac_shift = np.dot(verts_cart, mat_R_inv) + PIVOT_FRAC
    verts_cart_shift = np.dot(verts_frac_shift, mat_R)

    # 7. Compute max radius ---------------------------------------------------
    pivot_cart = np.dot(PIVOT_FRAC, mat_R)
    max_radius = float(np.linalg.norm(verts_cart_shift - pivot_cart, axis=1).max())

    # 8. Build CIF -----------------------------------------------------------
    bz = Structure(lat_R, ["C"], [PIVOT_FRAC])
    for v in verts_cart_shift:
        bz.append("X", v, coords_are_cartesian=True)

    outfile = Path(f"{ns.output}.cif")
    cif_str = bz.to(fmt="cif")  # returns CIF as string

    with open(outfile, "w", encoding="utf-8") as fh:
        fh.write(f"# max_BZ_radius  {max_radius:.6f}  Angstrom^-1\n")
        fh.write(cif_str)

    print(
        f"[bzvis] {outfile} written – {len(verts_cart_shift)} vertices; max R = {max_radius:.6f} Å⁻¹; "
        f"real {real_letter} → reciprocal {recip_letter}."
    )


if __name__ == "__main__":
    main(_get_parser().parse_args())
