#!/usr/bin/env python3
# Copyright 2026 Riccardo Savio
# SPDX-License-Identifier: Apache-2.0  (see LICENSE and NOTICE)
# https://github.com/riccardosavio02/academic-portfolio
"""
generate_inp.py -- Abaqus/Explicit deck generator for the NbC/Al2024 Whipple
shield of Moore & Teixeira-Dias (2025), Target 2.

One script for three sphere discretisations:

    --sphere ogrid      C3D8R O-grid (cubed sphere), exact surface r = R
    --sphere cartesian  C3D8R Cartesian-clip, uniform spacing dx = 2R/n_div
    --sphere particles  PC3D particles from t = 0 on the Cartesian-clip grid

Usage:
    python generate_inp.py CONFIG.yaml [--sphere MESH] [--out JOB.inp]
                           [--set KEY=VALUE ...]

--set overrides one config entry without editing the file, e.g.
    --set material_al2024.damage_initiation=johnson_cook
    --set simulation.time=5.0e-7

Optional keys and their defaults:

    geometry.sphere.mesh                        ogrid | cartesian | particles
    geometry.sphere.section                     own (default) | shared_al: the sphere is
        written inside the Al plate part, in the same elset and *Solid Section
        (AL_CTRL, so the same threshold). Abaqus: particles from different section
        definitions never interact; with a shared section the sphere particles and
        the Al plate particles do. Not available with mesh=particles.
    geometry.nbc.section                        own (default) | shared_al: the NbC plate is
        written inside the Al plate part and shares ONE *Solid Section with it
        (elset ALL_BODIES = ALL_NBC + ALL_AL, AL_CTRL). One section has one material,
        so the NbC plate becomes Al2024 (density included) unless the user material
        handles both (material_al2024.model: vumat_multi); with sphere.section=shared_al
        the three bodies are one section definition
    geometry.front                             nbc (default) | al: which plate faces the
        projectile. al swaps the plates (Al2024 in front, NbC behind);
        thicknesses and meshes stay with their plate
    simulation.boundary                         encastre (default) | free
    simulation.surface_split                    1 (default) | N: each plate surface (NBC_SURF, AL_SURF) written as
        N surfaces NBC_SURF_k / AL_SURF_k on slabs of whole element layers (assembly elsets by label range), each
        with exterior and interior faces, all of them in the general contact. With fine meshes (Al plate
        510 x 510 x 21) Abaqus stops the input processor with "THE SIZE OF SURFACE AL_SURF EXCEEDS THE STORAGE
        LIMITATIONS OF ABAQUS. SUBDIVISION OF THE SURFACE IS REQUIRED". Only with the self inclusions of the
        eroding contact (no sphere pair, no plate pair, no particles, no particle_thickness).
    simulation.eroding_contact                  false (default) | true
    simulation.eroding_surface                  interior (default) | eroding: how the eroding
        contact surfaces are written. interior = exterior faces + static INTERIOR faces;
        eroding = one dynamically evolving ERODING surface per body, which the manual
        prefers ("approximately the same results", less memory, faster; Modeling Guide,
        Element-Based Surface Definition)
    simulation.scale_factor                     none (default) | X: *Dynamic, Explicit, scale factor=X (the stable
        time increment is multiplied by X). Together with the linear bulk viscosity x10 it stabilises
        small benchmark runs with conversion; alone it does not cure a blow-up
    simulation.penalty_scale                    none (default) | X: *Contact Controls Assignment, type=SCALE PENALTY
        with factor X on the whole general contact domain (x10 cut the ETOTAL drift of a reduced model from
        +12.6 % to +2.6 %)
    simulation.particle_thickness               none (default) | {sphere: mm, nbc: mm, al: mm}: writes
        *Surface Property Assignment, property=THICKNESS on the parent surfaces SPH_SURF,
        NBC_SURF, AL_SURF, which sets the contact thickness of the particles generated
        from their elements (manual: with 1 or 2 particles per direction the default is
        "barely touching the closest face of the parent element", so PPD 2 halves the
        contact radius). Needs eroding_contact (and sphere_contact for the sphere)
    simulation.sphere_contact                   false (default) | true: with eroding_contact,
        adds the sphere surface (exterior + interior faces) and the inclusions
        sphere/NbC and sphere/Al. Abaqus 2025, "Finite element conversion to SPH
        particles", Table 4: particles inherit contact only from PAIRS of
        element-based surfaces (and from ALL EXTERIOR, exterior faces only), so
        without these pairs the sphere and its particles never see the faces
        opened inside the plates
    material_al2024.damage_initiation           legacy_ductile (default) | johnson_cook | none: no *Damage
        Initiation / *Damage Evolution at all (Al JC plasticity + Mie-Grueneisen only)
    material_al2024.damage_evolution            {type, value} (default) | none: no *Damage
        Evolution, so no element deletion (it cannot be combined with conversion)
    material_al2024.conversion                  strain (default) | user: Moore's criterion
        for the Al plate and the sphere, D_JC >= 1 OR sigma_1 >= 2.6 GPa, computed by
        Subroutines/vusdfld_conv.f (*Depvar, convert=2 with 4 state variables,
        *User Defined Field, *Initial Conditions, type=SOLUTION with the flag at 1).
        Needs conversion_criterion: USER in section_controls.al/.sphere and the
        job launched with user=vusdfld_conv.f
    material_nbc.delete_on_failure              false (default) | true: *Depvar, delete=8
    material_nbc.convert_on_failure             false (default) | true: *Depvar, convert=8:
        the NbC converts when ABQ_JH2 sets MPSTATUS (SDV8) to 0, i.e. at JH-2 failure
        (Moore: eps_f = 0.1 = FS, 5th value of the 3rd user_material row). Needs
        conversion_criterion: USER in section_controls.nbc
    material_nbc.model                          jh2 (default) | al2024_jc: the NbC plate
        gets the Al2024 Johnson-Cook + Mie-Grueneisen block (NbC density, no
        damage, no ABQ_JH2): diagnostic only, tests whether JH2 particles interact
        | al2024: the NbC section references the Al2024 material itself (Al density
        and damage): control of geometry.nbc.section=shared_al, same material, own section
    section_controls.<part>.element_conversion   true (default) | false: part stays FEM
    section_controls.<part>.particles_per_direction      1 (default)
    section_controls.<part>.kernel                       cubic (default) | quadratic | quintic
    section_controls.<part>.sph_formulation              classical (default) | nsph | xsph
    section_controls.<part>.smoothing_length             constant (default) | variable
    section_controls.<part>.smoothing_length_factor      scales the automatic smoothing length
        (second data line, second field; Abaqus 2025 keyword reference). Default: blank
    section_controls.<part>.smoothing_length_value       user-specified smoothing length in mm
        (second data line, first field) instead of the automatic one. Default: blank
    section_controls.<part>.tensile_instability_control  false (default) | true
        writes SPH TENSILE INSTABILITY CONTROL=YES (Monaghan 2000 artificial
        stress); the manual restricts it to particles that start uniformly spaced

Every *Section Controls carries sph conversion=PER ELEMENT: without it Abaqus
2025 reads element conversion=YES as a background-grid conversion and rejects
the conversion criterion.

Intentional differences from the legacy generators:
  * *Eos Shear, type=ELASTIC is written as *Elastic, type=SHEAR (same data;
    R2019x already flagged the old form for removal).
  * Every number is written exactly as it appears in the YAML. The legacy
    scripts rounded the conversion threshold to 3 decimals (0.0005 -> 0.001)
    and hard-coded the damage-evolution energy.
  * The deck is pure ASCII with LF line endings.
  * The heading states the boundary conditions that are actually applied.
  * Sets and surfaces that nothing references are no longer written.

Units: mm, tonne, s, MPa.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass

import yaml

SPHERE_MESHES = ("ogrid", "cartesian", "particles")
DAMAGE_MODELS = ("legacy_ductile", "johnson_cook", "none")
NBC_MODELS = ("jh2", "al2024_jc", "al2024")
BOUNDARIES = ("encastre", "free")

SPHERE_FIRST_ELEMENT = 2000000

# Solution-dependent state variables exposed by the built-in ABQ_JH2 VUMAT
JH2_STATE_VARIABLES = [
    ("PEEQ", "Equiv plastic strain"),
    ("PEEQ_RATE", "Equiv plastic strain rate"),
    ("DUCTCRT", "Ductile damage criterion"),
    ("DAMAGE", "Damage variable"),
    ("DELTAP", "Pressure increment bulking"),
    ("YIELD", "Yield strength"),
    ("EVOL", "Volumetric strain mu"),
    ("MPSTATUS", "Material point status"),
]


# -- formatting --------------------------------------------------------------

def num(value):
    """Exact text for a real number, always with a '.' or an exponent."""
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"non-finite number in the config: {value!r}")
    mantissa, e, exponent = repr(value).partition("e")
    if "." not in mantissa:
        mantissa += "."
    return mantissa + e + exponent


def nums(values):
    return ", ".join(num(v) for v in values)


def id_lines(ids, per_line=8):
    ids = list(ids)
    for i in range(0, len(ids), per_line):
        yield " " + ", ".join(str(v) for v in ids[i:i + per_line])


def yes_no(flag):
    return "YES" if flag else "NO"


class Deck:
    """Line writer for the .inp: ASCII only, LF line endings.

    The deck is written to a temporary file that replaces the target only once
    it is complete, so a failure never leaves a truncated .inp behind.
    """

    def __init__(self, path):
        self.path = path
        self._tmp_path = path + ".part"
        self._fh = None
        self.n_lines = 0

    def __enter__(self):
        self._fh = open(self._tmp_path, "w", encoding="ascii", newline="\n")
        return self

    def __exit__(self, exc_type, exc, tb):
        self._fh.close()
        if exc_type is None:
            os.replace(self._tmp_path, self.path)
        else:
            os.remove(self._tmp_path)
        return False

    def __call__(self, *lines):
        for line in lines:
            self._fh.write(line)
            self._fh.write("\n")
        self.n_lines += len(lines)


# -- plates ------------------------------------------------------------------

def _grid(lo, hi, n):
    return [lo + (hi - lo) * i / n for i in range(n + 1)]


class BoxMesh:
    """Structured block of nx * ny * nz C3D8R elements.

    Node IDs run x fastest, then y, then z, starting at first_node. Element IDs
    follow the same order starting at first_elem.
    """

    def __init__(self, x_range, y_range, z_range, nx, ny, nz, first_node, first_elem):
        self.nx, self.ny, self.nz = int(nx), int(ny), int(nz)
        self.first_node, self.first_elem = first_node, first_elem
        self.xs = _grid(*x_range, self.nx)
        self.ys = _grid(*y_range, self.ny)
        self.zs = _grid(*z_range, self.nz)
        self.last_node = first_node + (self.nx + 1) * (self.ny + 1) * (self.nz + 1) - 1
        self.last_elem = first_elem + self.nx * self.ny * self.nz - 1

    @property
    def n_nodes(self):
        return self.last_node - self.first_node + 1

    @property
    def n_elems(self):
        return self.last_elem - self.first_elem + 1

    def node(self, i, j, k):
        return self.first_node + i + (self.nx + 1) * (j + (self.ny + 1) * k)

    def layer_range(self, k):
        """First and last node ID of the k-th z level (the IDs are contiguous)."""
        first = self.node(0, 0, k)
        return first, first + (self.nx + 1) * (self.ny + 1) - 1

    def lateral_nodes(self):
        """Nodes on the four lateral faces (x and y at their bounds), ascending."""
        return [self.node(i, j, k)
                for k in range(self.nz + 1)
                for j in range(self.ny + 1)
                for i in range(self.nx + 1)
                if i in (0, self.nx) or j in (0, self.ny)]

    def write_nodes(self, w):
        nid = self.first_node
        for z in self.zs:
            for y in self.ys:
                for x in self.xs:
                    w(f" {nid}, {x:.5f}, {y:.5f}, {z:.5f}")
                    nid += 1

    def write_elements(self, w):
        sx = self.nx + 1
        sxy = sx * (self.ny + 1)
        eid = self.first_elem
        for k in range(self.nz):
            for j in range(self.ny):
                for i in range(self.nx):
                    n1 = self.node(i, j, k)
                    w(f" {eid}, {n1}, {n1 + 1}, {n1 + 1 + sx}, {n1 + sx}, "
                      f"{n1 + sxy}, {n1 + 1 + sxy}, {n1 + 1 + sx + sxy}, {n1 + sx + sxy}")
                    eid += 1


# -- sphere ------------------------------------------------------------------

@dataclass
class SphereMesh:
    nodes: dict          # node ID -> (x, y, z)
    elements: dict       # element ID -> 8-node tuple (C3D8R) or node ID (PC3D)
    element_type: str
    tip_node: int        # lowest node: the first to reach the target
    spacing: float = 0.0 # particle spacing dx, PC3D only


def _lowest_node(nodes):
    return min(nodes, key=lambda nid: nodes[nid][2])


def ogrid_sphere(radius, n, alpha, z_centre, first_node, first_elem):
    """O-grid (cubed sphere): an inner cube [-a, a]^3 with a = alpha * R, wrapped
    by six radial blocks whose outer face lies exactly on r = R.

    Element sizes are not uniform (edge ratio about 6 at n_div = 8), and so
    neither are the SPH smoothing lengths after conversion.
    """
    a = radius * alpha
    node_ids = {}
    nodes = {}

    def node_at(x, y, z):
        key = (round(x, 9), round(y, 9), round(z, 9))
        if key not in node_ids:
            nid = first_node + len(node_ids)
            node_ids[key] = nid
            nodes[nid] = (x, y, z + z_centre)
        return node_ids[key]

    def jacobian(conn):
        p = [nodes[nid] for nid in conn]
        v1 = [p[1][c] - p[0][c] for c in range(3)]
        v2 = [p[3][c] - p[0][c] for c in range(3)]
        v3 = [p[4][c] - p[0][c] for c in range(3)]
        return (v1[0] * (v2[1] * v3[2] - v2[2] * v3[1])
                - v1[1] * (v2[0] * v3[2] - v2[2] * v3[0])
                + v1[2] * (v2[0] * v3[1] - v2[1] * v3[0]))

    def positive(conn):
        """Reorder the connectivity so that the element is not inverted."""
        if jacobian(conn) > 0:
            return conn
        n1, n2, n3, n4, n5, n6, n7, n8 = conn
        for alt in ((n1, n4, n3, n2, n5, n8, n7, n6),
                    (n5, n6, n7, n8, n1, n2, n3, n4),
                    (n5, n8, n7, n6, n1, n4, n3, n2)):
            if jacobian(alt) > 0:
                return alt
        return conn

    def cube_node(ix, iy, iz):
        return node_at(-a + ix * 2 * a / n, -a + iy * 2 * a / n, -a + iz * 2 * a / n)

    elements = {}
    eid = first_elem

    # inner cube
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                elements[eid] = positive((
                    cube_node(ix, iy, iz), cube_node(ix + 1, iy, iz),
                    cube_node(ix + 1, iy + 1, iz), cube_node(ix, iy + 1, iz),
                    cube_node(ix, iy, iz + 1), cube_node(ix + 1, iy, iz + 1),
                    cube_node(ix + 1, iy + 1, iz + 1), cube_node(ix, iy + 1, iz + 1)))
                eid += 1

    # six radial blocks: (normal axis, normal sign, first and second in-face axes)
    for fax, fsgn, sax, tax in ((0, +1, 1, 2), (0, -1, 1, 2), (1, +1, 0, 2),
                                (1, -1, 0, 2), (2, +1, 0, 1), (2, -1, 0, 1)):
        for iu in range(n):
            for iv in range(n):
                s0 = -a + iu * 2 * a / n
                s1 = s0 + 2 * a / n
                t0 = -a + iv * 2 * a / n
                t1 = t0 + 2 * a / n

                def blend(s, t, k, fax=fax, fsgn=fsgn, sax=sax, tax=tax):
                    """Point at fraction k/n from the cube face to the sphere."""
                    p = [0.0, 0.0, 0.0]
                    p[fax] = fsgn * a
                    p[sax] = s
                    p[tax] = t
                    norm = math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
                    sp = [p[c] / norm * radius for c in range(3)]
                    tk = k / n
                    return node_at(p[0] * (1 - tk) + sp[0] * tk,
                                   p[1] * (1 - tk) + sp[1] * tk,
                                   p[2] * (1 - tk) + sp[2] * tk)

                for k in range(n):
                    b00 = blend(s0, t0, k)
                    b10 = blend(s1, t0, k)
                    b11 = blend(s1, t1, k)
                    b01 = blend(s0, t1, k)
                    t00 = blend(s0, t0, k + 1)
                    t10 = blend(s1, t0, k + 1)
                    t11 = blend(s1, t1, k + 1)
                    t01 = blend(s0, t1, k + 1)
                    elements[eid] = positive((b00, b10, b11, b01, t00, t10, t11, t01))
                    eid += 1

    return SphereMesh(nodes, elements, "C3D8R", _lowest_node(nodes))


def _cells_inside(radius, n):
    """Cells of the grid [-R, R]^3 (spacing 2R/n) whose centre lies within r <= R."""
    dx = 2.0 * radius / n
    r2 = radius * radius
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                cx = -radius + (ix + 0.5) * dx
                cy = -radius + (iy + 0.5) * dx
                cz = -radius + (iz + 0.5) * dx
                if cx * cx + cy * cy + cz * cz <= r2:
                    yield ix, iy, iz, cx, cy, cz


def cartesian_sphere(radius, n, z_centre, first_node, first_elem):
    """Cartesian-clip C3D8R sphere: the cells of a uniform grid whose centre is
    inside the sphere. Staircase surface, but dx = dy = dz everywhere, hence
    uniform SPH smoothing lengths after conversion.
    """
    dx = 2.0 * radius / n
    node_ids = {}
    nodes = {}

    def grid_node(ix, iy, iz):
        key = (ix, iy, iz)
        if key not in node_ids:
            nid = first_node + len(node_ids)
            node_ids[key] = nid
            nodes[nid] = (-radius + ix * dx, -radius + iy * dx, -radius + iz * dx + z_centre)
        return node_ids[key]

    elements = {}
    eid = first_elem
    for ix, iy, iz, _, _, _ in _cells_inside(radius, n):
        elements[eid] = (
            grid_node(ix, iy, iz), grid_node(ix + 1, iy, iz),
            grid_node(ix + 1, iy + 1, iz), grid_node(ix, iy + 1, iz),
            grid_node(ix, iy, iz + 1), grid_node(ix + 1, iy, iz + 1),
            grid_node(ix + 1, iy + 1, iz + 1), grid_node(ix, iy + 1, iz + 1))
        eid += 1

    return SphereMesh(nodes, elements, "C3D8R", _lowest_node(nodes))


def particle_sphere(radius, n, z_centre, first_node, first_elem):
    """PC3D sphere, SPH from t = 0: one particle at the centre of every cell of
    the Cartesian-clip grid, one single-node element per particle.
    """
    nodes = {}
    elements = {}
    nid, eid = first_node, first_elem
    for _, _, _, cx, cy, cz in _cells_inside(radius, n):
        nodes[nid] = (cx, cy, cz + z_centre)
        elements[eid] = nid
        nid += 1
        eid += 1

    return SphereMesh(nodes, elements, "PC3D", _lowest_node(nodes),
                      spacing=2.0 * radius / n)


def build_sphere(mesh, g, z_centre, first_node, first_elem):
    radius, n = g["radius"], int(g["n_div"])
    if mesh == "ogrid":
        return ogrid_sphere(radius, n, g["alpha"], z_centre, first_node, first_elem)
    if mesh == "cartesian":
        return cartesian_sphere(radius, n, z_centre, first_node, first_elem)
    return particle_sphere(radius, n, z_centre, first_node, first_elem)


# -- keyword blocks ------------------------------------------------------------

def section_controls(sc):
    if not sc["element_conversion"]:
        # no FEM-to-SPH conversion: the part stays finite elements (failed
        # elements are removed by the material's own deletion, if any)
        return [
            f"*Section Controls, name={sc['name']}, "
            f"distortion control={yes_no(sc['distortion_control'])}, "
            f"length ratio={num(sc['length_ratio'])}",
            " 1., 1., 1.",
        ]
    per_direction = int(sc.get("particles_per_direction", 1))
    # Monaghan (2000) artificial stress; the manual says it is valid only for
    # particles that start from a uniform distribution
    tic = (", sph tensile instability control=YES"
           if sc.get("tensile_instability_control", False) else "")
    sph = ""
    for key, keyword, choices in (("kernel", "kernel", ("cubic", "quadratic", "quintic")),
                                  ("sph_formulation", "sph formulation", ("classical", "nsph", "xsph")),
                                  ("smoothing_length", "sph smoothing length", ("constant", "variable"))):
        if key in sc:
            _check_choice(f"section_controls.{sc['name']}.{key}", str(sc[key]).lower(), choices)
            sph += f", {keyword}={str(sc[key]).upper()}"
    # second data line: smoothing length, smoothing length factor, minimum particles, -, XSPH coefficient
    if "smoothing_length_factor" in sc or "smoothing_length_value" in sc:
        value = num(sc["smoothing_length_value"]) if "smoothing_length_value" in sc else ""
        scale = num(sc["smoothing_length_factor"]) if "smoothing_length_factor" in sc else ""
        factor = f" {value}, {scale}, , , , ,"
    else:
        factor = " , , , , , ,"
    return [
        f"*Section Controls, name={sc['name']}, "
        f"distortion control={yes_no(sc['distortion_control'])}, "
        f"length ratio={num(sc['length_ratio'])}, "
        f"element conversion=YES, "
        f"conversion criterion={sc['conversion_criterion']}, "
        f"sph conversion=PER ELEMENT{tic}{sph}",
        " 1., 1., 1.",
        factor,
        " , , , , , ,",
        # USER: the material decides (*Depvar, convert=n), no threshold on the data line
        f" {per_direction}," if str(sc["conversion_criterion"]).upper() == "USER"
        else f" {per_direction}, {num(sc['sct'])}",
    ]


def vumat_multi_material(name, mat_al, mat_nbc, nsolid, ppd, first_label, nbc_first, nbc_last, transfer, fs=None,
                         sph_first=0, body=False):
    """ONE user material for Al2024 AND NbC: vumat_hvi.f, material name prefix
    HVI. props 1-16 Al2024 (G rho0 c0 s Gamma0 A B n C edot0 D1-D5 sigma_max), 17 Nsolid, 18 PPD, 19 first element
    label, 20-21 first and last label of the NbC elements, 22 energy transfer, 23-24 unused, 25-56 the 32 NbC
    constants in the ABQ_JH2 order (material_nbc.user_material, 8 per row; FS = 5th of the 3rd row, overridden by
    material_nbc.fs: negative = convert at fracture, D = 1). The branch of each point is chosen by the label of its
    element (of its parent for a particle: labels must be contiguous). Density by field variable 1 (0 = Al2024,
    1 = NbC)."""
    jc, rd, dd, eos = (mat_al["johnson_cook_plastic"], mat_al["rate_dependent"], mat_al["damage_ductile"],
                       mat_al["eos_usup"])
    al = [mat_al["eos_shear_G"], mat_al["density"], eos["c0"], eos["S1"], eos["gamma0"],
          jc["A"], jc["B"], jc["n"], rd["C"], rd["edot0"],
          dd["d1"], dd["d2"], dd["d3"], dd["d4"], dd["d5"], mat_al["damage_maxs"]["s11"]]
    rows = [list(r) + [0.0] * (8 - len(r)) for r in mat_nbc["user_material"]["values"]]
    while len(rows) < 4:
        rows.append([0.0] * 8)
    if fs is not None:
        rows[2][4] = fs
    # 23 = first label of the sphere: with body=True the VUMAT writes SDV9 BODY (1 sphere, 2 NbC, 3 Al2024) in every
    # point, particles included, to colour the three bodies of one section apart in the viewer
    # 24 = eps_c: Al2024 also converts when EPSP reaches it; 0 = off (material_al2024.eps_c)
    props = al + [nsolid, ppd, first_label, nbc_first, nbc_last, 1 if transfer else 0, sph_first if body else 0,
                  float(mat_al.get("eps_c", 0.0))] \
        + [v for r in rows[:4] for v in r]
    return ["**",
            f"** ONE user material for Al2024 and NbC (vumat_hvi.f): Al2024 + sphere + NbC plate in one section, "
            f"NbC elements {nbc_first}..{nbc_last}, Nsolid {nsolid}",
            f"*Material, name={name}",
            "*Density, dependencies=1",
            f" {num(mat_al['density'])}, , 0.",
            f" {num(mat_nbc['density'])}, , 1.",
            "** SDV4 CONV: the element converts to SPH particles when the VUMAT sets it to 0",
            "*Depvar, convert=4", " 9," if body else " 8,",
            ' 1, EPSP, "equivalent plastic strain"',
            ' 2, DMG, "damage D"',
            ' 3, EPSV, "accumulated volumetric log strain"',
            ' 4, CONV, "conversion flag, 1 = stays FEM"',
            ' 5, SIG1, "largest principal stress"',
            ' 6, PRESS, "pressure"',
            ' 7, DPBULK, "bulking pressure (NbC)"',
            ' 8, EPRATE, "plastic strain rate (NbC)"'] + \
        ([' 9, BODY, "1 sphere, 2 NbC, 3 Al2024"'] if body else []) + \
        [f"*User Material, constants={len(props)}"] + \
        [" " + nums(props[k:k + 8]) for k in range(0, len(props), 8)]


def vumat_hvi_material(name, mat_al, nsolid, ppd, first_label, transfer):
    """Al2024 as the user material of vumat_hvi.f (verified on single-element tests): Mie-Grueneisen, JC without thermal term, JC damage of Moore 2025 (d1-d5 of damage_ductile), conversion
    when D >= 1 or sigma_1 >= damage_maxs.s11 (2600, He 2020), particles = failed material. Props 1-17 as in the
    README of vumat_hvi; with transfer, 18 PPD and 19 first element label: the particle takes its parent's energy
    and state at birth (parent of particle L = first_label + (L - nsolid - 1) / PPD^3)."""
    jc, rd, dd, eos = (mat_al["johnson_cook_plastic"], mat_al["rate_dependent"], mat_al["damage_ductile"],
                       mat_al["eos_usup"])
    props = [mat_al["eos_shear_G"], mat_al["density"], eos["c0"], eos["S1"], eos["gamma0"],
             jc["A"], jc["B"], jc["n"], rd["C"], rd["edot0"],
             dd["d1"], dd["d2"], dd["d3"], dd["d4"], dd["d5"], mat_al["damage_maxs"]["s11"], nsolid]
    if transfer:
        props += [ppd, first_label]
    return ["**",
            f"** Al2024 user material (vumat_hvi.f): elements {first_label}..{nsolid}, particles above {nsolid}",
            f"*Material, name={name}",
            "*Density", f" {num(mat_al['density'])},",
            "** SDV4 CONV: the element converts to SPH particles when the VUMAT sets it to 0",
            "*Depvar, convert=4", " 6,",
            ' 1, EPSP, "equivalent plastic strain"',
            ' 2, DMG, "JC damage D"',
            ' 3, EPSV, "accumulated volumetric log strain"',
            ' 4, CONV, "conversion flag, 1 = stays FEM"',
            ' 5, SIG1, "largest principal stress"',
            ' 6, PRESS, "pressure"',
            "** G, rho0, c0, s, Gamma0, A, B, n / C, edot0, D1-D5, sigma_max / Nsolid"
            + (", PPD, first element label" if transfer else ""),
            f"*User Material, constants={len(props)}",
            " " + nums(props[0:8]), " " + nums(props[8:16]), " " + nums(props[16:])]


def _check_choice(key, value, choices):
    if value not in choices:
        raise ValueError(f"{key}: {value!r} is not one of {', '.join(choices)}")


# -- deck ----------------------------------------------------------------------

def generate(cfg, sphere_mesh, out_path, source="config"):
    geo, sim = cfg["geometry"], cfg["simulation"]
    mat_nbc, mat_al = cfg["material_nbc"], cfg["material_al2024"]
    ctrl, out_req = cfg["section_controls"], cfg["output_request"]
    g_nbc, g_al, g_sph = geo["nbc"], geo["al"], geo["sphere"]

    boundary = sim.get("boundary", "encastre")
    eroding_contact = bool(sim.get("eroding_contact", False))
    sphere_contact = bool(sim.get("sphere_contact", False))
    eroding_surface = sim.get("eroding_surface", "interior")
    _check_choice("simulation.eroding_surface", eroding_surface, ("interior", "eroding"))
    damage =mat_al.get("damage_initiation", "legacy_ductile")
    _check_choice("simulation.boundary", boundary, BOUNDARIES)
    _check_choice("material_al2024.damage_initiation", damage, DAMAGE_MODELS)
    nbc_model = mat_nbc.get("model", "jh2")
    _check_choice("material_nbc.model", nbc_model, NBC_MODELS)
    nbc_shared = g_nbc.get("section", "own") == "shared_al"
    _check_choice("geometry.nbc.section", g_nbc.get("section", "own"), ("own", "shared_al"))
    if nbc_shared:
        nbc_model = "al2024"   # one *Solid Section has one material
    nbc_name = {"jh2": mat_nbc["name"], "al2024_jc": "NBC_AS_AL2024_JC",
                "al2024": mat_al["name"]}[nbc_model]
    nbc_inst = "AL-1" if nbc_shared else "NBC-1"
    if (g_al["lx"], g_al["ly"]) != (g_nbc["lx"], g_nbc["ly"]):
        raise ValueError("geometry: the Al plate must have the same lx and ly as the NbC plate")
    front = geo.get("front", "nbc")
    _check_choice("geometry.front", front, ("nbc", "al"))
    particles = sphere_mesh == "particles"
    shared = g_sph.get("section", "own") == "shared_al"
    _check_choice("geometry.sphere.section", g_sph.get("section", "own"), ("own", "shared_al"))
    if shared and particles:
        raise ValueError("geometry.sphere.section=shared_al needs a C3D8R sphere (ogrid or cartesian)")
    sph_inst = "AL-1" if shared else "SPH-1"
    if sphere_contact and (not eroding_contact or shared or particles):
        raise ValueError("simulation.sphere_contact needs eroding_contact: true and a C3D8R "
                         "sphere in its own section")
    notes = []
    # material_al2024.model: native (Abaqus JC + EoS + damage cards), vumat_hvi or vumat_multi (user material)
    al_model = str(mat_al.get("model", "native")).lower()
    _check_choice("material_al2024.model", al_model, ("native", "vumat_hvi", "vumat_multi"))
    al_vumat = al_model == "vumat_hvi"
    al_multi = al_model == "vumat_multi"
    if al_multi:
        # F5: one material, one section for the three bodies (labels contiguous: NbC, Al2024, sphere)
        if not (shared and nbc_shared) or particles or sphere_contact:
            raise ValueError("material_al2024.model: vumat_multi needs geometry.nbc.section = shared_al, "
                             "geometry.sphere.section = shared_al, a C3D8R sphere and sphere_contact: false")
        if str(ctrl["al"]["conversion_criterion"]).upper() != "USER":
            raise ValueError("material_al2024.model: vumat_multi needs section_controls.al.conversion_criterion: USER")
        notes.append("ONE user material (HVI_MULTI) for NbC, Al2024 and sphere in ONE section: launch with "
                     "user=vumat_hvi.f")
    egap = 1 if al_multi else 10        # element label gap between the bodies (multi: contiguous labels)
    al_transfer = bool(mat_al.get("energy_transfer", True))
    if al_vumat:
        if shared or particles or nbc_shared:
            raise ValueError("material_al2024.model: vumat_hvi needs the plates and the C3D8R sphere in their "
                             "own sections (one material each, own Nsolid)")
        for p in ("al", "sphere"):
            if str(ctrl[p]["conversion_criterion"]).upper() != "USER":
                raise ValueError(f"material_al2024.model: vumat_hvi needs section_controls.{p}."
                                 "conversion_criterion: USER")
        notes.append("Al2024 plate and sphere = vumat_hvi user material: launch with "
                     "user=vumat_hvi.f"
                     + ("" if al_transfer else " (no energy transfer at birth)"))
    al_mat_plate = "HVI_MULTI" if al_multi else ("AL2024_USER_PLATE" if al_vumat else mat_al["name"])
    al_mat_sph = "AL2024_USER_SPH" if al_vumat else mat_al["name"]
    if shared:
        notes.append("sphere in part PLATE_AL, elset ALL_AL, section controls AL_CTRL "
                     f"(sct {ctrl['al']['sct']}); section_controls.sphere ignored")
    if nbc_shared:
        notes.append("NbC plate in part PLATE_AL, one *Solid Section (elset ALL_BODIES) with "
                     f"material {mat_al['name']} and AL_CTRL (sct {ctrl['al']['sct']}): the NbC "
                     "plate IS Al2024; material_nbc and section_controls.nbc ignored")
    elif nbc_model == "al2024":
        notes.append(f"NbC plate with the {mat_al['name']} material in its own section "
                     "(material_nbc.model: al2024): diagnostic deck")

    # meshes
    hx, hy = g_nbc["lx"] / 2, g_nbc["ly"] / 2
    # the front plate sits in z >= 0, the rear one in z <= 0, the sphere above the front plate
    z_nbc = (0.0, g_nbc["thickness"]) if front == "nbc" else (-g_nbc["thickness"], 0.0)
    z_al = (-g_al["thickness"], 0.0) if front == "nbc" else (0.0, g_al["thickness"])
    nbc = BoxMesh((-hx, hx), (-hy, hy), z_nbc,
                  g_nbc["nx"], g_nbc["ny"], g_nbc["nz"], 1, 1)
    al = BoxMesh((-hx, hx), (-hy, hy), z_al,
                 g_al["nx"], g_al["ny"], g_al["nz"], nbc.last_node + 10, nbc.last_elem + egap)
    z_centre = (g_nbc if front == "nbc" else g_al)["thickness"] + g_sph["gap"] + g_sph["radius"]
    # node layers of the two plates on the interface z = 0, tied together
    nbc_if, al_if = ("NBC_BOT", "AL_TOP") if front == "nbc" else ("NBC_TOP", "AL_BOT")
    if front == "al":
        notes.append("plates swapped: Al2024 in front, NbC behind (geometry.front: al)")
    # own section: element labels are local to the SPHERE part, any range works.
    # shared section: the sphere lives in part PLATE_AL, so its labels must not collide with the
    # plate's (labels 2000000+ would overwrite Al elements of the plate).
    sph_first = al.last_elem + egap if shared else SPHERE_FIRST_ELEMENT
    sphere = build_sphere(sphere_mesh, g_sph, z_centre, al.last_node + 10, sph_first)
    sph_first_node, sph_last_node = min(sphere.nodes), max(sphere.nodes)
    sph_first_elem, sph_last_elem = min(sphere.elements), max(sphere.elements)
    if shared and sph_first_elem <= al.last_elem:
        raise ValueError("shared_al: sphere element labels overlap the Al plate's")

    print(f"NbC    : {nbc.n_nodes} nodes ({nbc.first_node}..{nbc.last_node}), "
          f"{nbc.n_elems} C3D8R")
    print(f"Al2024 : {al.n_nodes} nodes ({al.first_node}..{al.last_node}), "
          f"{al.n_elems} C3D8R")
    spacing = f", dx={sphere.spacing:.4f} mm" if particles else ""
    print(f"Sphere : {len(sphere.nodes)} nodes ({sph_first_node}..{sph_last_node}), "
          f"{len(sphere.elements)} {sphere.element_type} ({sphere_mesh}){spacing}, "
          f"tip node {sphere.tip_node}")

    sphere_text = {
        "ogrid": "C3D8R O-grid (cubed sphere), exact surface r = R, FEM-to-SPH conversion",
        "cartesian": "C3D8R Cartesian-clip, uniform dx = 2R/n_div, FEM-to-SPH conversion",
        "particles": "PC3D particles from t = 0, Cartesian-clip grid dx = 2R/n_div",
    }[sphere_mesh]
    boundary_text = ("lateral plate edges ENCASTRE (Moore 2025 Sec. 4.2 applies none)"
                     if boundary == "encastre" else
                     "none, free plates (as in Moore 2025 Sec. 4.2)")
    source = source.encode("ascii", "replace").decode("ascii")

    with Deck(out_path) as w:
        w("*Heading",
          " Moore & Teixeira-Dias (2025) Target 2: NbC/Al2024 Whipple shield",
          f"** Generated by generate_inp.py from {source}, sphere mesh: {sphere_mesh}",
          "** NbC   : C3D8R, FEM-to-SPH conversion, JH2 (built-in ABQ_JH2 VUMAT)",
          "** Al2024: C3D8R, FEM-to-SPH conversion, native Johnson-Cook + Mie-Grueneisen EoS",
          f"** Sphere: Al2024, {sphere_text}",
          f"** Boundary conditions: {boundary_text}",
          "** Units: mm, tonne, s, MPa",
          "*Preprint, echo=NO, model=NO, history=NO, contact=NO")

        # -- part: NbC plate (or its sets only, inside PLATE_AL with nbc.section=shared_al)
        bot_first, bot_last = nbc.layer_range(0 if front == "nbc" else nbc.nz)
        nbc_sets = ("*Nset, nset=ALL_NBC, generate", f" {nbc.first_node}, {nbc.last_node}, 1",
                    "*Elset, elset=ALL_NBC, generate", f" {nbc.first_elem}, {nbc.last_elem}, 1",
                    f"*Nset, nset={nbc_if}, generate", f" {bot_first}, {bot_last}, 1",
                    f"*Surface, type=NODE, name={nbc_if}_SURF", f" {nbc_if},")
        if not nbc_shared:
            w("**",
              "** NbC bumper plate" if front == "nbc" else "** NbC plate (rear: geometry.front = al)",
              "*Part, name=PLATE_NBC",
              "*Node")
            nbc.write_nodes(w)
            w("*Element, type=C3D8R")
            nbc.write_elements(w)
            w(*nbc_sets)
            w(f"*Solid Section, elset=ALL_NBC, material={nbc_name}, "
              f"controls={ctrl['nbc']['name']}",
              " 1.0,",
              "*End Part")

        # -- part: Al2024 plate
        w("**",
          "** Al2024 rear plate" if front == "nbc" else "** Al2024 plate (front: geometry.front = al)",
          "*Part, name=PLATE_AL",
          "*Node")
        al.write_nodes(w)
        if nbc_shared:
            w("** NbC plate nodes: same part and section as the Al2024 plate", "*Node")
            nbc.write_nodes(w)
        if shared:
            w(f"** Sphere nodes ({sphere_text}), same part and section as the plate", "*Node")
            for nid, (x, y, z) in sphere.nodes.items():
                w(f" {nid}, {x:.5f}, {y:.5f}, {z:.5f}")
        w("*Element, type=C3D8R")
        al.write_elements(w)
        if nbc_shared:
            w("** NbC plate elements", "*Element, type=C3D8R")
            nbc.write_elements(w)
        if shared:
            w("** Sphere elements", "*Element, type=C3D8R")
            for eid, conn in sphere.elements.items():
                w(f" {eid}, " + ", ".join(str(nid) for nid in conn))
        top_first, top_last = al.layer_range(al.nz if front == "nbc" else 0)
        w("*Nset, nset=ALL_AL, generate", f" {al.first_node}, {al.last_node}, 1")
        if shared:
            w(f" {sph_first_node}, {sph_last_node}, 1")
        w("*Elset, elset=ALL_AL, generate", f" {al.first_elem}, {al.last_elem}, 1")
        if shared:
            w(f" {sph_first_elem}, {sph_last_elem}, 1",
              "*Elset, elset=ALL_SPH, generate", f" {sph_first_elem}, {sph_last_elem}, 1")
        w(f"*Nset, nset={al_if}, generate", f" {top_first}, {top_last}, 1",
          f"*Surface, type=NODE, name={al_if}_SURF", f" {al_if},")
        section_elset = "ALL_AL"
        if nbc_shared:
            w(*nbc_sets)
            w("** ONE section definition for every body in this part: particles of different",
              "** section definitions never interact (Abaqus 2025, SPH limitations)",
              "*Elset, elset=ALL_BODIES", " ALL_NBC, ALL_AL")
            section_elset = "ALL_BODIES"
        w(f"*Solid Section, elset={section_elset}, material={al_mat_plate}, "
          f"controls={ctrl['al']['name']}",
          " 1.0,",
          "*End Part")

        # -- part: sphere (own section)
        if not shared:
          w("**",
          f"** Sphere: {sphere_text}",
          "*Part, name=SPHERE",
          "*Node")
          for nid, (x, y, z) in sphere.nodes.items():
              w(f" {nid}, {x:.5f}, {y:.5f}, {z:.5f}")
          w(f"*Element, type={sphere.element_type}")
          if particles:
              for eid, nid in sphere.elements.items():
                  w(f" {eid}, {nid}")
          else:
              for eid, conn in sphere.elements.items():
                  w(f" {eid}, " + ", ".join(str(nid) for nid in conn))
          w("*Elset, elset=ALL_SPH, generate", f" {sph_first_elem}, {sph_last_elem}, 1")
          if particles:
              w("** PC3D section data: half of the particle spacing",
                f"*Solid Section, elset=ALL_SPH, material={mat_al['name']}",
                f" {num(sphere.spacing / 2)},")
          else:
              w(f"*Solid Section, elset=ALL_SPH, material={al_mat_sph}, "
                f"controls={ctrl['sphere']['name']}",
                " 1.0,")
          w("*End Part")

        # -- assembly
        w("**",
          "*Assembly, name=Assembly",
          *(() if nbc_shared else ("*Instance, name=NBC-1, part=PLATE_NBC", "*End Instance")),
          "*Instance, name=AL-1, part=PLATE_AL", "*End Instance")
        if not shared:
            w("*Instance, name=SPH-1, part=SPHERE", "*End Instance")
        w("**",
          f"*Nset, nset=TIP_NODE, instance={sph_inst}", f" {sphere.tip_node},",
          f"*Nset, nset=ALL_SPH_N, instance={sph_inst}, generate",
          f" {sph_first_node}, {sph_last_node}, 1")
        def body_surface(name, elset):
            if eroding_surface == "eroding":
                return (f"*Surface, name={name}, type=ELEMENT", f" {elset}, ERODING")
            return (f"*Surface, name={name}, type=ELEMENT", f" {elset},", f" {elset}, INTERIOR")

        split = int(sim.get("surface_split", 1))
        plate_surfs = []
        if split > 1 and (particles or sphere_contact or sim.get("plate_pair", True)
                          or sim.get("particle_thickness")):
            raise ValueError("simulation.surface_split > 1 works only with the eroding-contact self inclusions: "
                             "no particles, sphere_contact, plate_pair or particle_thickness")
        if eroding_contact or particles:
            if eroding_surface == "eroding":
                w("** Plate surfaces as dynamically evolving ERODING surfaces: they follow the",
                  "** currently exposed faces of the elements not yet failed or converted")
            else:
                w("** Plate surfaces with exterior AND interior faces: faces exposed when",
                  "** elements fail or convert can then take part in contact")
            if split > 1:
                w(f"** surfaces split into {split} slabs of whole element layers (simulation.surface_split)")
                for base, inst, plate in (("NBC_SURF", nbc_inst, nbc), ("AL_SURF", "AL-1", al)):
                    per = plate.nx * plate.ny
                    cuts = [round(k * plate.nz / split) for k in range(split + 1)]
                    for k in range(split):
                        nm = f"{base}_{k + 1}"
                        w(f"*Elset, elset={nm}_E, instance={inst}, generate",
                          f" {plate.first_elem + cuts[k] * per}, {plate.first_elem + cuts[k + 1] * per - 1}, 1",
                          *body_surface(nm, f"{nm}_E"))
                        plate_surfs.append(nm)
            else:
                w(*body_surface("NBC_SURF", f"{nbc_inst}.ALL_NBC"),
                  *body_surface("AL_SURF", "AL-1.ALL_AL"))
                plate_surfs += ["NBC_SURF", "AL_SURF"]
        if sphere_contact:
            w(*body_surface("SPH_SURF", "SPH-1.ALL_SPH"))
        if particles:
            w("*Surface, name=SPH_NODE_SURF, type=NODE", " ALL_SPH_N,")
        w("**",
          "** Tie: NbC bottom face (secondary) to Al2024 top face (main)" if front == "nbc" else
          "** Tie: NbC top face (secondary) to Al2024 bottom face (main)",
          "*Tie, name=NBC_AL_TIE, adjust=YES",
          f" {nbc_inst}.{nbc_if}_SURF, AL-1.{al_if}_SURF")
        if boundary == "encastre":
            w(f"*Nset, nset=NBC_EDGES_ASM, instance={nbc_inst}")
            w(*id_lines(nbc.lateral_nodes()))
            w("*Nset, nset=AL_EDGES_ASM, instance=AL-1")
            w(*id_lines(al.lateral_nodes()))
        w("*End Assembly")

        # -- material: NbC
        jc, rd = mat_al["johnson_cook_plastic"], mat_al["rate_dependent"]
        eos = mat_al["eos_usup"]
        if nbc_model == "al2024_jc":
            notes.append("NbC plate with the Al2024 JC + EoS material (NbC density), "
                         "no JH2 and no damage: diagnostic deck")
            w("**",
              "** NbC DIAGNOSTIC: Al2024 Johnson-Cook + Mie-Grueneisen with the NbC density,",
              "** no damage (material_nbc.model: al2024_jc). Not a physical NbC model.",
              f"*Material, name={nbc_name}",
              "*Density", f" {num(mat_nbc['density'])},",
              "*Plastic, hardening=JOHNSON COOK",
              " " + nums([jc["A"], jc["B"], jc["n"], jc["m"], jc["Tmelt"], jc["Troom"]]),
              "*Rate Dependent, type=JOHNSON COOK",
              " " + nums([rd["C"], rd["edot0"]]),
              "*Eos, type=USUP",
              " " + nums([eos["c0"], eos["S1"], eos["gamma0"]]),
              "*Elastic, type=SHEAR",
              f" {num(mat_al['eos_shear_G'])},")
        um = mat_nbc["user_material"]
        if nbc_model == "jh2":
            n_constants = int(um["constants"])
            rows = um["values"]
            if len(rows) > math.ceil(n_constants / 8) or any(len(row) > 8 for row in rows):
                raise ValueError("material_nbc.user_material: at most 8 values per row and "
                                 "ceil(constants / 8) rows")
            nbc_delete = bool(mat_nbc.get("delete_on_failure", False))
            nbc_convert = bool(mat_nbc.get("convert_on_failure", False))
            if nbc_delete and nbc_convert:
                raise ValueError("material_nbc: delete_on_failure and convert_on_failure are "
                                 "mutually exclusive (*Depvar DELETE and CONVERT)")
            if nbc_convert and str(ctrl["nbc"]["conversion_criterion"]).upper() != "USER":
                raise ValueError("material_nbc.convert_on_failure needs "
                                 "section_controls.nbc.conversion_criterion: USER")
            nbc_depvar = ("*Depvar, delete=8" if nbc_delete else
                          "*Depvar, convert=8" if nbc_convert else "*Depvar")
            if int(mat_nbc["depvar"]) != len(JH2_STATE_VARIABLES):
                raise ValueError(f"material_nbc.depvar: ABQ_JH2 has {len(JH2_STATE_VARIABLES)} "
                                 "state variables")
            w("**",
              "** NbC: Johnson-Holmquist 2 through the built-in ABQ_JH2 VUMAT, selected by the",
              "** ABQ_JH2 prefix of the material name (no user subroutine to compile)",
              f"*Material, name={nbc_name}",
              "*Density", f" {num(mat_nbc['density'])},",
              # delete_on_failure: remove the element when ABQ_JH2 sets MPSTATUS (SDV8) to 0;
              # convert_on_failure: convert it to particles instead (CONVERSION CRITERION=USER)
              nbc_depvar,
              f" {len(JH2_STATE_VARIABLES)},")
            w(*(f' {i}, {key}, "{text}"' for i, (key, text) in enumerate(JH2_STATE_VARIABLES, 1)))
            w("** One YAML row per data line. Abaqus reads 8 constants per line, so the row",
              "** layout decides which constant each value is: do not merge the rows.",
              f"*User Material, constants={n_constants}")
            w(*(" " + nums(row) for row in rows))

        # -- material: Al2024
        jc, rd = mat_al["johnson_cook_plastic"], mat_al["rate_dependent"]
        dd, de, eos = mat_al["damage_ductile"], mat_al["damage_evolution"], mat_al["eos_usup"]
        d1_to_d5 = [dd[k] for k in ("d1", "d2", "d3", "d4", "d5")]
        no_evolution = str(de).lower() == "none"
        al_conv = str(mat_al.get("conversion", "strain")).lower()
        _check_choice("material_al2024.conversion", al_conv, ("strain", "user"))
        al_user = al_conv == "user"
        if al_user:
            parts = ["al"] + ([] if (shared or particles) else ["sphere"])
            for p in parts:
                if str(ctrl[p]["conversion_criterion"]).upper() != "USER":
                    raise ValueError(f"material_al2024.conversion: user needs "
                                     f"section_controls.{p}.conversion_criterion: USER")
            if not no_evolution:
                notes.append("material_al2024.conversion: user with *Damage Evolution: element "
                             "deletion and conversion cannot be combined; set damage_evolution: none")
            notes.append("Al/sphere conversion by vusdfld_conv.f: launch with user=vusdfld_conv.f")
        w("**",
          f"*Material, name={mat_al['name']}",
          "*Density", f" {num(mat_al['density'])},",
          "** Johnson-Cook plasticity: A, B, n, m, T_melt, T_transition",
          "*Plastic, hardening=JOHNSON COOK",
          " " + nums([jc["A"], jc["B"], jc["n"], jc["m"], jc["Tmelt"], jc["Troom"]]),
          "*Rate Dependent, type=JOHNSON COOK",
          " " + nums([rd["C"], rd["edot0"]]))
        if damage == "none":
            w("** No damage model on the Al2024 (material_al2024.damage_initiation: none): plasticity and EoS only")
        elif damage == "johnson_cook":
            if "damage_maxs" in mat_al:
                notes.append("damage_maxs is not written: MAXS is a traction-separation "
                             "criterion (cohesive/XFEM), inert on C3D8R")
            w("** Johnson-Cook damage initiation: d1..d5, T_melt, T_transition, reference strain rate",
              "*Damage Initiation, criterion=JOHNSON COOK",
              " " + nums(d1_to_d5 + [jc["Tmelt"], jc["Troom"], rd["edot0"]]))
            if no_evolution:
                w("** No *Damage Evolution: initiation only (JCCRT as an indicator), no element",
                  "** deletion, which cannot be combined with conversion")
            else:
                w(f"*Damage Evolution, type={de['type']}",
                  f" {num(de['value'])},")
        else:
            dm = mat_al["damage_maxs"]
            w("** LEGACY damage block (damage_initiation: legacy_ductile), kept for older decks.",
              "** Known issues:",
              "**  - criterion=DUCTILE reads d1..d5 as ONE tabular point of the fracture",
              "**    strain (eps_f = d1 at triaxiality d2, strain rate d3, ...)",
              "**  - *Damage Evolution belongs to the initiation just above it, i.e. MAXS,",
              "**    a traction-separation criterion that is inert on C3D8R elements",
              "*Damage Initiation, criterion=DUCTILE",
              " " + nums(d1_to_d5),
              "*Damage Initiation, criterion=MAXS",
              " " + nums([dm["s11"], dm["s22"], dm["s33"]]))
            if not no_evolution:
                w(f"*Damage Evolution, type={de['type']}",
                  f" {num(de['value'])},")
        if al_user:
            w("** Moore's conversion for the Al plate and the sphere, computed by vusdfld_conv.f:",
              "** convert when D_JC >= 1 or sigma_1 >= 2.6 GPa (flag SDV2: 1 = FEM, 0 = convert)",
              "*Depvar, convert=2",
              " 4,",
              ' 1, JCDMG, "JC damage accumulation D"',
              ' 2, CONVFLG, "FEM to SPH conversion flag"',
              ' 3, SIGMAX1, "Max principal stress"',
              ' 4, PEEQOLD, "PEEQ previous increment"',
              "*User Defined Field")
        w("** Mie-Grueneisen EoS, linear Us-Up Hugoniot: c0, s, Gamma0",
          "*Eos, type=USUP",
          " " + nums([eos["c0"], eos["S1"], eos["gamma0"]]),
          "** Deviatoric response of the EoS material: shear modulus",
          "*Elastic, type=SHEAR",
          f" {num(mat_al['eos_shear_G'])},")
        if "specific_heat" in mat_al or "inelastic_heat_fraction" in mat_al:
            w("** Thermal cards: Abaqus ignores them unless the step is adiabatic.",
              "** Thermal effects are excluded on purpose (plain *Dynamic, Explicit).")
            if "specific_heat" in mat_al:
                w("*Specific Heat", f" {num(mat_al['specific_heat'])},")
            if "inelastic_heat_fraction" in mat_al:
                w("*Inelastic Heat Fraction", f" {num(mat_al['inelastic_heat_fraction'])},")
        if al_multi:
            # the native Al card above stays in the deck unused (the single section refers to HVI_MULTI)
            w(*vumat_multi_material("HVI_MULTI", mat_al, mat_nbc, sph_last_elem,
                                    int(ctrl["al"].get("particles_per_direction", 1)), nbc.first_elem,
                                    nbc.first_elem, nbc.last_elem, al_transfer, mat_nbc.get("fs"),
                                    sph_first_elem, bool(mat_al.get("body_sdv", False))))
        if al_vumat:
            # the native card above stays in the deck unused (no section refers to it)
            ppd_al = int(ctrl["al"].get("particles_per_direction", 1))
            ppd_sph = int(ctrl["sphere"].get("particles_per_direction", 1))
            w(*vumat_hvi_material(al_mat_plate, mat_al, al.last_elem, ppd_al, al.first_elem, al_transfer))
            w(*vumat_hvi_material(al_mat_sph, mat_al, sph_last_elem, ppd_sph, sph_first_elem, al_transfer))

        # -- FEM-to-SPH conversion
        w("**",
          "** FEM-to-SPH conversion. Run with cpus=1 and no domains=: Abaqus 2025",
          "** (no Feature Pack) crashes in the packager with >= 2 domains when the",
          "** converting sections differ.",
          "** Last data line: particles per isoparametric direction, conversion threshold")
        controlled = [] if nbc_shared else [("nbc", "NbC plate")]
        controlled.append(("al", "Al2024 plate"))
        if not particles and not shared:
            controlled.append(("sphere", "sphere"))
        for key, label in controlled:
            w(f"** {label}")
            w(*section_controls(ctrl[key]))

        # -- interaction and initial state
        w("**",
          "*Surface Interaction, name=FRIC",
          "*Friction", f" {num(sim['friction'])},",
          "**",
          "*Initial Conditions, type=VELOCITY",
          f" ALL_SPH_N, 3, {num(sim['velocity'])}")
        if al_user:
            # SDVs start at 0: without this every Al element would be flagged "convert" at t = 0
            # Abaqus reads, for EVERY element line, as many values as the largest number of state variables
            # of any material in the model (Keywords Guide, *INITIAL CONDITIONS, TYPE=SOLUTION): 8 with the
            # ABQ_JH2 NbC plate, 4 otherwise. First line: set label + up to 7 values; then up to 8 per line.
            n_sdv = max(4, len(JH2_STATE_VARIABLES) if nbc_model == "jh2" else 4)
            sdv0 = [0.0, 1.0, 0.0, 0.0] + [0.0] * (n_sdv - 4)

            def sdv_lines(label):
                out = [f" {label}, " + ", ".join(num(v) for v in sdv0[:7])]
                rest = sdv0[7:]
                while rest:
                    out.append(" " + ", ".join(num(v) for v in rest[:8]))
                    rest = rest[8:]
                return out

            w("** State variables of vusdfld_conv.f: JCDMG, CONVFLG (1 = stays FEM), SIGMAX1, PEEQOLD,",
              f"** then zeros up to the {n_sdv} values Abaqus expects per element line",
              "*Initial Conditions, type=SOLUTION",
              *sdv_lines("AL-1.ALL_AL"))
            if not (shared or particles):
                w(*sdv_lines("SPH-1.ALL_SPH"))
        if al_multi:
            # one material: CONV (SDV4) = 1 on every element (8 values per element line) and the density label
            # (field variable 1: 1 on the NbC nodes, 0 on the others; the plates have separate nodes, tied)
            w("** HVI_MULTI state variables: EPSP, DMG, EPSV, CONV (= 1), SIG1, PRESS, DPBULK, EPRATE",
              "*Initial Conditions, type=SOLUTION",
              " AL-1.ALL_BODIES, 0., 0., 0., 1., 0., 0., 0.", " 0., 0." if mat_al.get("body_sdv", False) else " 0.",
              "** material label for the density: 1 on the NbC nodes, 0 on the Al2024 plate and the sphere",
              "*Initial Conditions, type=FIELD, variable=1",
              " AL-1.ALL_AL, 0.", " AL-1.ALL_NBC, 1.")
        if al_vumat:
            # vumat_hvi: SDV4 CONV = 1 at t = 0, otherwise every element converts at the first increment. Abaqus
            # reads as many values per element line as the largest number of SDVs in the model (8 with ABQ_JH2)
            n_sdv = max(6, len(JH2_STATE_VARIABLES) if nbc_model == "jh2" else 6)
            sdv0 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0] + [0.0] * (n_sdv - 6)

            def sdv_lines_um(label):
                out = [f" {label}, " + ", ".join(num(v) for v in sdv0[:7])]
                rest = sdv0[7:]
                while rest:
                    out.append(" " + ", ".join(num(v) for v in rest[:8]))
                    rest = rest[8:]
                return out

            w("** vumat_hvi state variables: EPSP, DMG, EPSV, CONV (= 1, stays FEM), SIG1, PRESS,",
              f"** then zeros up to the {n_sdv} values Abaqus expects per element line",
              "*Initial Conditions, type=SOLUTION",
              *sdv_lines_um("AL-1.ALL_AL"), *sdv_lines_um("SPH-1.ALL_SPH"))

        # -- step
        b1, b2 = sim["bulk_viscosity"]
        w("**",
          "*Step, name=Impact, nlgeom=YES",
          "*Dynamic, Explicit" + (f", scale factor={num(float(sim['scale_factor']))}"
                                  if sim.get("scale_factor") else ""), f" , {num(sim['time'])}",
          *(["** Variable mass scaling (simulation.mass_scaling_dt): mass is added ONLY to the elements and particles whose",
             "** stable increment falls below the threshold (otherwise single particles drive dt to ~4e-11 s)",
             f"*Variable Mass Scaling, dt={num(float(sim['mass_scaling_dt']))}, type=below min, frequency=1"]
            if sim.get("mass_scaling_dt") else []),
          "*Bulk Viscosity", f" {num(b1)}, {num(b2)}")
        if boundary == "encastre":
            w("*Boundary", " NBC_EDGES_ASM, ENCASTRE", " AL_EDGES_ASM, ENCASTRE")
        w("*Contact",
          "*Contact Inclusions, ALL EXTERIOR")
        if particles:
            w("** PC3D particles are not part of ALL EXTERIOR: explicit pairs",
              "** (self-contact is not allowed for node-based surfaces)",
              "*Contact Inclusions", " SPH_NODE_SURF, NBC_SURF",
              "*Contact Inclusions", " SPH_NODE_SURF, AL_SURF")
        if eroding_contact:
            w("** Eroding contact: ALL EXTERIOR only covers the initial skin of the mesh.",
              "** With the interior faces, crater walls opened by failure or conversion",
              "** take part in contact (without them the kinetic energy plateaus)",
              *[x for nm in plate_surfs for x in ("*Contact Inclusions", f" {nm},")])
            if sim.get("plate_pair", True):
                w("*Contact Inclusions", " NBC_SURF, AL_SURF")
            else:
                w("** NO NBC/Al pair: one section, the SPH kernel couples the bodies (a pair on top "
                  "blows up); the plates stay tied")
        if sphere_contact:
            w("** Sphere against the plates, interior faces included on both sides: the pairs also",
              "** generate the particle (node-based) inclusions, see the manual's Table 4",
              "*Contact Inclusions", " SPH_SURF, NBC_SURF",
              "*Contact Inclusions", " SPH_SURF, AL_SURF")
        w("*Contact Property Assignment", " ,  , FRIC",
          "*Contact Controls Assignment, nodal erosion=YES")
        if sim.get("penalty_scale"):
            # k x10 cut the ETOTAL drift of a reduced model from +12.6 % to
            # +2.6 % at no cost; x0.3 blew up. Data line: surface 1, surface 2 (blank = whole domain), factor
            w("** Penalty stiffness of the whole general contact domain scaled (simulation.penalty_scale)",
              "*Contact Controls Assignment, type=SCALE PENALTY",
              f" , , {num(float(sim['penalty_scale']))}")
        thick = sim.get("particle_thickness", {}) or {}
        if thick:
            surf_of = {"sphere": "SPH_SURF", "nbc": "NBC_SURF", "al": "AL_SURF"}
            lines = []
            for body, value in thick.items():
                if body not in surf_of:
                    raise ValueError(f"simulation.particle_thickness.{body}: use sphere, nbc or al")
                if body == "sphere" and not sphere_contact:
                    raise ValueError("simulation.particle_thickness.sphere needs sphere_contact: true "
                                     "(the sphere needs its element-based surface SPH_SURF)")
                if not (eroding_contact or particles):
                    raise ValueError("simulation.particle_thickness needs eroding_contact: true "
                                     "(NBC_SURF and AL_SURF)")
                lines.append(f" {surf_of[body]}, {num(float(value))},")
            w("** Contact thickness of the generated particles (Abaqus 2025, FE conversion, Interactions:",
              "** with 1 or 2 particles per direction it is set to barely touch the parent faces, i.e. it",
              "** shrinks with the particle spacing; this assignment on the parent surface overrides it)",
              "*Surface Property Assignment, property=THICKNESS", *lines)

        # -- output
        def field_only(requested, where):
            # Abaqus 2025 Output Variable Identifiers: JCCRT is History-only ("Field: no"); asking for it
            # in the field output aborts the input processor ("OUTPUT REQUEST JCCRT IS NOT AVAILABLE").
            # DMICRT (all active damage initiation criteria) is the field variable.
            keep = [v for v in requested if v != "JCCRT"]
            if len(keep) != len(requested):
                notes.append(f"output_request.{where}: JCCRT dropped (history-only variable); "
                             "use DMICRT for the field output")
            return keep
        field_al = field_only(list(out_req["field_al"]), "field_al")
        out_req = dict(out_req, field_sphere=field_only(list(out_req["field_sphere"]), "field_sphere"),
                       field_nbc=field_only(list(out_req["field_nbc"]), "field_nbc"))
        if al_vumat or al_multi:
            # user material: no PEEQ / TEMP from Abaqus, the state lives in the SDVs
            def user_fields(req):
                keep = [v for v in req if v not in ("PEEQ", "TEMP", "DMICRT")]
                return keep + ([] if "SDV" in keep else ["SDV"])
            field_al = user_fields(field_al)
            out_req = dict(out_req, field_sphere=user_fields(out_req["field_sphere"]))
        # SDV exists only with ABQ_JH2: an Al2024 NbC plate gets the Al output
        field_nbc = out_req["field_nbc"] if nbc_model == "jh2" else field_al
        w("**",
          f"*Output, field, number interval={int(cfg['output']['field_intervals'])}",
          f"*Element Output, elset={nbc_inst}.ALL_NBC", " " + ", ".join(field_nbc),
          "*Element Output, elset=AL-1.ALL_AL", " " + ", ".join(field_al),
          f"*Element Output, elset={sph_inst}.ALL_SPH", " " + ", ".join(out_req["field_sphere"]),
          "*Node Output", " " + ", ".join(out_req["node"]),
          f"*Output, history, frequency={int(out_req['history_freq'])}",
          "*Node Output, nset=TIP_NODE", " U3, V3",
          "*Energy Output", " " + ", ".join(out_req["energy"]),
          "**",
          "*End Step")

    for note in notes:
        print(f"note   : {note}")
    print(f"Written: {out_path}  ({w.n_lines} lines)")


def apply_override(cfg, assignment):
    """Set one config entry from a KEY=VALUE string; KEY is dotted, VALUE is YAML."""
    key, sep, raw = assignment.partition("=")
    if not sep or not key:
        raise ValueError(f"--set expects KEY=VALUE, got {assignment!r}")
    *parents, leaf = key.split(".")
    node = cfg
    for part in parents:
        node = node.setdefault(part, {})
    node[leaf] = yaml.safe_load(raw)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate the Abaqus/Explicit .inp of the NbC/Al2024 Whipple shield.")
    parser.add_argument("config", help="YAML configuration file")
    parser.add_argument("--sphere", choices=SPHERE_MESHES,
                        help="sphere discretisation (default: geometry.sphere.mesh)")
    parser.add_argument("--out", help="output deck (default: <output.job_name>.inp)")
    parser.add_argument("--set", dest="overrides", action="append", default=[],
                        metavar="KEY=VALUE",
                        help="override a config entry, e.g. simulation.time=5.0e-7")
    args = parser.parse_args(argv)

    with open(args.config, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    try:
        for assignment in args.overrides:
            apply_override(cfg, assignment)
        sphere_mesh = args.sphere or cfg["geometry"]["sphere"].get("mesh")
        if sphere_mesh not in SPHERE_MESHES:
            parser.error("choose the sphere mesh with --sphere {"
                         + ",".join(SPHERE_MESHES) + "} or geometry.sphere.mesh")
        out_path = args.out or cfg["output"]["job_name"] + ".inp"
        generate(cfg, sphere_mesh, out_path, source=os.path.basename(args.config))
    except KeyError as exc:
        sys.exit(f"error: missing config key {exc}")
    except ValueError as exc:
        sys.exit(f"error: {exc}")


if __name__ == "__main__":
    main()
