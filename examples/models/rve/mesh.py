#!/usr/bin/env python3
"""
Periodic 2D microstructure mesh (square RVE with circular fibres) using Gmsh.

- Domain: [0, L] x [0, L]
- Fibres: circles (disks) inside the unit cell; can be generated randomly.
- Geometry: matrix = square minus fibre disks (holes); fibre boundaries are kept as separate physical group.
- Periodicity: left<->right and bottom<->top via gmsh.model.mesh.setPeriodic
- Output: .msh (and optionally .geo_unrolled)

Tested with Gmsh Python API (4.x).
"""

import math
import random
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
import gmsh


@dataclass
class Fibre:
    x: float
    y: float
    r: float


def dist2(a: Tuple[float, float], b: Tuple[float, float]) -> float:
    return (a[0]-b[0])**2 + (a[1]-b[1])**2


def place_random_fibres(
    L: float,
    r: float,
    n: int,
    min_gap: float = 0.02,
    seed: int = 1,
    periodic: bool = True,
    max_tries: int = 200000,
) -> List[Fibre]:
    """
    Place n identical radius-r fibres.
    If periodic=True, enforce non-overlap under periodic distance.
    """
    rng = random.Random(seed)
    fibres: List[Fibre] = []

    def periodic_d2(p, q):
        # minimum-image convention
        dx = p[0] - q[0]
        dy = p[1] - q[1]
        dx -= L * round(dx / L)
        dy -= L * round(dy / L)
        return dx*dx + dy*dy

    min_center_dist = 2.0 * r + min_gap

    tries = 0
    while len(fibres) < n and tries < max_tries:
        tries += 1
        x = rng.random() * L
        y = rng.random() * L
        ok = True
        for f in fibres:
            d2 = periodic_d2((x, y), (f.x, f.y)) if periodic else dist2((x, y), (f.x, f.y))
            if d2 < (min_center_dist ** 2):
                ok = False
                break
        if ok:
            fibres.append(Fibre(x, y, r))

    if len(fibres) < n:
        raise RuntimeError(
            f"Could only place {len(fibres)}/{n} fibres. "
            f"Try smaller r, fewer fibres, or smaller min_gap."
        )

    return fibres


def add_disk_occurrences_for_periodicity(L: float, fibres: List[Fibre]) -> List[Fibre]:
    """
    If a disk intersects a boundary, we add periodic image(s) shifted by +/-L
    so that the boolean cut respects periodic geometry.

    We keep only those image disks whose *center* lies within an expanded window
    [-r, L+r] so that they can affect the unit cell cut.
    """
    out = []
    for f in fibres:
        shifts = [(0.0, 0.0)]
        # If disk crosses left/right boundary, include images shifted in x
        if f.x - f.r < 0.0:
            shifts.append((+L, 0.0))
        if f.x + f.r > L:
            shifts.append((-L, 0.0))
        # If disk crosses bottom/top, include images shifted in y
        if f.y - f.r < 0.0:
            shifts.append((0.0, +L))
        if f.y + f.r > L:
            shifts.append((0.0, -L))
        # Corners: if crosses in both x and y, include combined shifts
        extra = []
        for sx, sy in shifts:
            for tx, ty in shifts:
                extra.append((sx + tx, sy + ty))
        # unique
        uniq = list({(sx, sy) for (sx, sy) in extra})

        for sx, sy in uniq:
            cx = f.x + sx
            cy = f.y + sy
            # keep if relevant for cutting the unit cell
            if (-f.r <= cx <= L + f.r) and (-f.r <= cy <= L + f.r):
                out.append(Fibre(cx, cy, f.r))
    return out


def build_periodic_rve(
    L: float,
    fibres: List[Fibre],
    mesh_size_matrix: float,
    mesh_size_fibre: float,
    algo: int = 6,  # 6 = Frontal-Delaunay for 2D (often good)
    recombine_quads: bool = False,
    save_geo: bool = False,
    msh_filename: str = "rve_periodic.msh",
):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("periodic_rve")

    occ = gmsh.model.occ

    # --- Base square
    square = occ.addRectangle(0.0, 0.0, 0.0, L, L)

    # --- Add disks (including periodic image occurrences when needed)
    disks = []
    fibre_occ = add_disk_occurrences_for_periodicity(L, fibres)
    for f in fibre_occ:
        disks.append(occ.addDisk(f.x, f.y, 0.0, f.r, f.r))

    # --- Boolean fragments: keep both matrix and fibre disks as separate surfaces
    # This will allow us to assign different physical groups (materials)
    if disks:
        fragments = occ.fragment([(2, square)], [(2, d) for d in disks])
        all_surfs = fragments[0]  # list of (dim, tag)
    else:
        all_surfs = [(2, square)]

    occ.synchronize()

    # --- Separate matrix and fibre surfaces, and remove any outside the square
    matrix_surfs = []
    fibre_surfs = []
    outside_surfs = []
    for (dim, tag) in all_surfs:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        # Matrix surface covers the whole domain
        if abs(xmin - 0.0) < 1e-6 and abs(ymin - 0.0) < 1e-6 and abs(xmax - L) < 1e-6 and abs(ymax - L) < 1e-6:
            matrix_surfs.append((dim, tag))
        # Fibre surfaces must be fully inside the square
        elif (xmin >= 0.0 - 1e-6 and ymin >= 0.0 - 1e-6 and xmax <= L + 1e-6 and ymax <= L + 1e-6):
            fibre_surfs.append((dim, tag))
        else:
            outside_surfs.append((dim, tag))

    # Remove any surfaces outside the square
    if outside_surfs:
        gmsh.model.removeEntities([s for s in outside_surfs], recursive=True)

    # --- Identify fibre boundary curves (circles) for physical group
    fibre_boundary_curves = set()
    for (dim, tag) in fibre_surfs:
        bnd = gmsh.model.getBoundary([(dim, tag)], oriented=False, recursive=False)
        for (bdim, btag) in bnd:
            if bdim == 1:
                fibre_boundary_curves.add(btag)

    # --- Outer boundary curves of the matrix surface(s)
    all_matrix_boundaries = set()
    for s in matrix_surfs:
        bnd = gmsh.model.getBoundary([s], oriented=False, recursive=False)
        for (dim, tag) in bnd:
            if dim == 1:
                all_matrix_boundaries.add(tag)

    def curve_bbox(curve_tag: int):
        return gmsh.model.getBoundingBox(1, curve_tag)

    tol = 1e-8
    left, right, bottom, top = [], [], [], []
    for c in all_matrix_boundaries:
        xmin, ymin, zmin, xmax, ymax, zmax = curve_bbox(c)
        # classify curves lying on the box boundaries
        if abs(xmin - 0.0) < tol and abs(xmax - 0.0) < tol:
            left.append(c)
        elif abs(xmin - L) < tol and abs(xmax - L) < tol:
            right.append(c)
        elif abs(ymin - 0.0) < tol and abs(ymax - 0.0) < tol:
            bottom.append(c)
        elif abs(ymin - L) < tol and abs(ymax - L) < tol:
            top.append(c)

    # --- Physical groups
    # Matrix surfaces:
    mat_tags = [tag for (dim, tag) in matrix_surfs if dim == 2]
    pg_matrix = gmsh.model.addPhysicalGroup(2, mat_tags)
    gmsh.model.setPhysicalName(2, pg_matrix, "matrix")

    # Fibre surfaces (2D):
    # Assign all circular inclusions to a single physical group
    fibre_tags = [tag for (dim, tag) in fibre_surfs if dim == 2]
    if fibre_tags:
        pg_fibre = gmsh.model.addPhysicalGroup(2, fibre_tags)
        gmsh.model.setPhysicalName(2, pg_fibre, "fibre")

    # Fibre boundaries (1D):
    if fibre_boundary_curves:
        pg_fibre_bnd = gmsh.model.addPhysicalGroup(1, sorted(fibre_boundary_curves))
        gmsh.model.setPhysicalName(1, pg_fibre_bnd, "FIBRE_BOUNDARY")

    # Outer boundaries:
    if left:
        pg_left = gmsh.model.addPhysicalGroup(1, left)
        gmsh.model.setPhysicalName(1, pg_left, "left")
    if right:
        pg_right = gmsh.model.addPhysicalGroup(1, right)
        gmsh.model.setPhysicalName(1, pg_right, "right")
    if bottom:
        pg_bottom = gmsh.model.addPhysicalGroup(1, bottom)
        gmsh.model.setPhysicalName(1, pg_bottom, "bottom")
    if top:
        pg_top = gmsh.model.addPhysicalGroup(1, top)
        gmsh.model.setPhysicalName(1, pg_top, "top")


    # --- Uniform mesh size everywhere
    # Use mesh_size_matrix as the global mesh size
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size_matrix)

    # --- Periodic boundary constraints (geometry-based)
    # Right is the "master" and Left is the "slave" with translation +L in x
    # Top is master and Bottom is slave with translation +L in y
    # The affineTransform is a 4x4 matrix in row-major order.
    if left and right:
        T_x = [1, 0, 0, L,
               0, 1, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 1]
        gmsh.model.mesh.setPeriodic(1, left, right, T_x)

    if bottom and top:
        T_y = [1, 0, 0, 0,
               0, 1, 0, L,
               0, 0, 1, 0,
               0, 0, 0, 1]
        gmsh.model.mesh.setPeriodic(1, bottom, top, T_y)

    # --- Meshing
    gmsh.option.setNumber("Mesh.Algorithm", algo)
    if recombine_quads:
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)  # simple
    gmsh.model.mesh.generate(2)

    if save_geo:
        gmsh.write("rve_periodic.geo_unrolled")

    gmsh.write(msh_filename)
    gmsh.finalize()


def main():
    # --- RVE params
    L = 1.0
    r = 0.10
    n_fibres = 6
    min_gap = 0.02

    fibres = place_random_fibres(L=L, r=r, n=n_fibres, min_gap=min_gap, seed=4, periodic=True)

    build_periodic_rve(
        L=L,
        fibres=fibres,
        mesh_size_matrix=0.03,  # Refined by factor 2
        mesh_size_fibre=0.01,  # Refined by factor 2
        algo=6,
        recombine_quads=False,
        save_geo=False,
        msh_filename="rve_periodic.msh",
    )

    print("Wrote: rve_periodic.msh")


if __name__ == "__main__":
    main()
