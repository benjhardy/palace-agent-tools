"""Gmsh (OpenCASCADE) helpers for RF models that Palace can read, plus a mesh report.

Conventions (the ones Palace and the lint rules expect):

* Every volume region is a 3-D physical group (a material); every surface that carries a boundary condition,
  port or lumped element is a 2-D physical group. Volume tags start at 1, surface tags at 101.
* Conductors much thicker than the skin depth are removed from the mesh; their surfaces become one group
  (PEC or Conductivity boundary). Zero-thickness traces are sheets (PEC/Conductivity on an internal face).
* Outer faces are grouped per domain side (``-x``, ``+x``, ``-y``, ``+y``, ``-z``, ``+z``), named by the
  ``walls`` map you pass to ``build``. A wall you leave out of the config is PMC (Palace's natural condition).
* The mesh is written as Gmsh MSH 2.2 with only physical groups saved (``Mesh.SaveAll = 0``), next to a
  ``<name>.tags.json`` with group names -> tags and the geometry the lint rules use (bounding boxes, areas,
  normals, patch axes, which faces are outer, domain and structure boxes) and the mesh report (tets, quality,
  exact Nedelec unknowns per order).

``gmsh`` is imported lazily, so the rest of the package works without it.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .estimate import nedelec_unknowns
from .rf import UNITS, wavelength

SIDES = ("-x", "+x", "-y", "+y", "-z", "+z")


def _gmsh():
    try:
        import gmsh
    except ImportError as e:  # pragma: no cover
        raise ImportError("palace_tools.mesh needs gmsh: pip install gmsh") from e
    return gmsh


@dataclass
class Item:
    name: str
    kind: str  # volume | air | conductor | sheet | port
    dimtags: list
    eps_r: float = 1.0
    group: str | None = None
    priority: int = 2
    structure: bool = True


@dataclass
class RFModel:
    """Build an RF model with named parts, then ``build()``, ``mesh()`` and ``write()``.

    Lengths are in ``units`` (written to tags.json; use the same units in the Palace config: Model.L0)."""
    name: str = "model"
    units: str = "mm"
    verbosity: int = 2
    items: list = field(default_factory=list)

    def __post_init__(self):
        g = _gmsh()
        if not g.isInitialized():
            g.initialize()
        g.option.setNumber("General.Terminal", 1 if self.verbosity > 0 else 0)
        g.option.setNumber("General.Verbosity", self.verbosity)
        g.model.add(self.name)
        self.g = g
        self.occ = g.model.occ
        self.groups: dict = {}
        self.built = False

    # ------------------------------------------------------------------ parts
    def _add(self, it: Item) -> str:
        if any(i.name == it.name for i in self.items):
            raise ValueError(f"duplicate part name {it.name!r}")
        self.items.append(it)
        return it.name

    def box(self, name, x0, y0, z0, dx, dy, dz, eps_r=1.0, kind="volume", group=None):
        """Dielectric (or air) box. kind: volume | air | conductor."""
        t = self.occ.addBox(x0, y0, z0, dx, dy, dz)
        pr = {"air": 1, "volume": 2, "conductor": 3}[kind]
        return self._add(Item(name, kind, [(3, t)], eps_r, group or name, pr, structure=kind != "air"))

    def substrate(self, name, x0, y0, x1, y1, z0, h, eps_r):
        return self.box(name, x0, y0, z0, x1 - x0, y1 - y0, h, eps_r)

    def cylinder(self, name, x, y, z0, r, height, eps_r=1.0, kind="volume", axis="z", group=None):
        d = {"x": (height, 0, 0), "y": (0, height, 0), "z": (0, 0, height)}[axis]
        t = self.occ.addCylinder(x, y, z0, *d, r)
        pr = {"air": 1, "volume": 2, "conductor": 3}[kind]
        return self._add(Item(name, kind, [(3, t)], eps_r, group or name, pr, structure=kind != "air"))

    def via(self, name, x, y, r, z0, z1, group="conductor"):
        """Plated via as a removed conductor; its wall joins the ``group`` surface set."""
        return self.cylinder(name, x, y, z0, r, z1 - z0, kind="conductor", group=group)

    def trace(self, name, x0, y0, x1, y1, z, t=0.0, group=None):
        """Rectangular trace in the plane z. t = 0: zero-thickness sheet; t > 0: removed conductor volume."""
        if t > 0:
            v = self.occ.addBox(x0, y0, z, x1 - x0, y1 - y0, t)
            return self._add(Item(name, "conductor", [(3, v)], 1.0, group or name, 3))
        s = self.occ.addRectangle(x0, y0, z, x1 - x0, y1 - y0)
        return self._add(Item(name, "sheet", [(2, s)], 1.0, group or name, 4))

    def polygon_trace(self, name, pts, z, t=0.0, group=None):
        """Polygon trace (list of (x, y)) at height z; holes: cut with booleans, never multi-loop surfaces
        (a multi-loop addPlaneSurface produced an invalid solid with the holes counted as material)."""
        p = [self.occ.addPoint(x, y, z) for x, y in pts]
        ln = [self.occ.addLine(p[i], p[(i + 1) % len(p)]) for i in range(len(p))]
        s = self.occ.addPlaneSurface([self.occ.addCurveLoop(ln)])
        if t > 0:
            ext = self.occ.extrude([(2, s)], 0, 0, t)
            vols = [dt for dt in ext if dt[0] == 3]
            return self._add(Item(name, "conductor", vols, 1.0, group or name, 3))
        return self._add(Item(name, "sheet", [(2, s)], 1.0, group or name, 4))

    def rect(self, name, plane: str, at: float, u: tuple, v: tuple, kind="port", group=None):
        """Axis-aligned rectangle in plane 'x', 'y' or 'z' = at, spanning u and v (the other two axes in xyz
        order). kind 'port' (lumped/wave port, lumped element) or 'sheet' (trace, resistive sheet)."""
        (u0, u1), (v0, v1) = u, v
        if plane == "z":
            c = [(u0, v0, at), (u1, v0, at), (u1, v1, at), (u0, v1, at)]
        elif plane == "y":
            c = [(u0, at, v0), (u1, at, v0), (u1, at, v1), (u0, at, v1)]
        elif plane == "x":
            c = [(at, u0, v0), (at, u1, v0), (at, u1, v1), (at, u0, v1)]
        else:
            raise ValueError("plane must be x, y or z")
        p = [self.occ.addPoint(*q) for q in c]
        ln = [self.occ.addLine(p[i], p[(i + 1) % 4]) for i in range(4)]
        s = self.occ.addPlaneSurface([self.occ.addCurveLoop(ln)])
        return self._add(Item(name, kind, [(2, s)], 1.0, group or name, 5))

    def air_box(self, name="air", pad=None, pad_lambda=None, f_ghz=None, z_min=None, eps_r=1.0, bbox=None):
        """Air (domain) box around everything added so far. pad: (x, y, z_top) or per-side dict
        {'-x':..,'+x':..,'-y':..,'+y':..,'-z':..,'+z':..}; pad_lambda: padding in free-space wavelengths at f_ghz
        (absorbing walls: >= 0.3-0.5 lambda0 was needed to stop fake wall loss). z_min clamps the bottom (e.g. a
        ground plane at z = 0)."""
        self.occ.synchronize()
        if bbox is None:
            bb = None
            for it in self.items:
                for dt in it.dimtags:
                    b = self.occ.getBoundingBox(*dt)
                    bb = list(b) if bb is None else [min(bb[i], b[i]) for i in range(3)] + [max(bb[i], b[i]) for i in range(3, 6)]
            if bb is None:
                raise ValueError("add parts before the air box, or pass bbox")
            if pad_lambda is not None:
                lam = wavelength(f_ghz) / UNITS[self.units]
                pad = {s: pad_lambda * lam for s in SIDES}
            if isinstance(pad, (int, float)):
                pad = {s: float(pad) for s in SIDES}
            elif isinstance(pad, (tuple, list)):
                pad = {"-x": pad[0], "+x": pad[0], "-y": pad[1], "+y": pad[1], "-z": pad[2], "+z": pad[2]}
            pad = {s: (pad or {}).get(s, 0.0) for s in SIDES}
            x0, y0, z0 = bb[0] - pad["-x"], bb[1] - pad["-y"], bb[2] - pad["-z"]
            x1, y1, z1 = bb[3] + pad["+x"], bb[4] + pad["+y"], bb[5] + pad["+z"]
            if z_min is not None:
                z0 = max(z0, z_min)
        else:
            x0, y0, z0, x1, y1, z1 = bbox
        return self.box(name, x0, y0, z0, x1 - x0, y1 - y0, z1 - z0, eps_r, kind="air")

    # ------------------------------------------------------------------ build
    def build(self, walls: dict | None = None) -> dict:
        """Fragment, remove conductors, classify surfaces and create physical groups.

        walls: side -> group name for outer faces, e.g. {'-z': 'ground', '+z': 'top', '-x': 'end', '+x': 'end'};
        unmapped sides get 'wall-x' style names. Returns the group table."""
        g, occ = self.g, self.occ
        walls = walls or {}
        vol_items = [i for i in self.items if i.dimtags[0][0] == 3]
        srf_items = [i for i in self.items if i.dimtags[0][0] == 2]
        objs = [dt for i in vol_items for dt in i.dimtags]
        tools = [dt for i in srf_items for dt in i.dimtags]
        if len(objs) + len(tools) > 1:
            _, omap = occ.fragment(objs, tools)
        else:
            omap = [[dt] for dt in objs + tools]
        occ.synchronize()
        k = 0
        owner: dict[int, Item] = {}
        for it in vol_items:
            for _ in it.dimtags:
                for d, t in omap[k]:
                    if d == 3 and (t not in owner or it.priority >= owner[t].priority):
                        owner[t] = it
                k += 1
        srf_parts: dict[str, set] = {}
        for it in srf_items:
            for _ in it.dimtags:
                srf_parts.setdefault(it.name, set()).update(t for d, t in omap[k] if d == 2)
                k += 1
        cond_vols = {t for t, it in owner.items() if it.kind == "conductor"}
        kept = {t: it for t, it in owner.items() if it.kind != "conductor"}
        cond_srf: dict[str, set] = {}
        for t in cond_vols:
            for _, s in g.model.getBoundary([(3, t)], oriented=False):
                cond_srf.setdefault(owner[t].group, set()).add(abs(s))
        if cond_vols:
            occ.remove([(3, t) for t in cond_vols])
            occ.synchronize()
        self._kept = kept
        dom = None
        for t in kept:
            b = g.model.getBoundingBox(3, t)
            dom = list(b) if dom is None else [min(dom[i], b[i]) for i in range(3)] + [max(dom[i], b[i]) for i in range(3, 6)]
        self.domain_bbox = dom
        span = max(dom[3] - dom[0], dom[4] - dom[1], dom[5] - dom[2])
        tol = 1e-7 * span

        def n_kept_adj(s):
            up, _ = g.model.getAdjacencies(2, s)
            return sum(1 for u in up if u in kept)

        def side_of(s):
            b = g.model.getBoundingBox(2, s)
            for ax in range(3):
                if abs(b[ax + 3] - b[ax]) < tol:
                    if abs(b[ax] - dom[ax]) < tol:
                        return "-" + "xyz"[ax]
                    if abs(b[ax] - dom[ax + 3]) < tol:
                        return "+" + "xyz"[ax]
            return None

        assigned: dict[int, str] = {}
        groups2: dict[str, list] = {}
        self._wall_groups = set()
        # conductors first, then ports / sheets, then the walls
        for grp, ss in cond_srf.items():
            for s in sorted(ss):
                if n_kept_adj(s) >= 1 and s not in assigned:
                    assigned[s] = grp
                    groups2.setdefault(grp, []).append(s)
        for it in srf_items:
            for s in sorted(srf_parts.get(it.name, ())):
                if n_kept_adj(s) >= 1:
                    if s in assigned and assigned[s] != it.group:
                        # a port/sheet drawn on a conductor face: the port wins (it is the more specific part)
                        groups2[assigned[s]].remove(s)
                    assigned[s] = it.group
                    groups2.setdefault(it.group, []).append(s)
        outer_other = []
        for t in kept:
            for _, s in g.model.getBoundary([(3, t)], oriented=False):
                s = abs(s)
                if s in assigned or n_kept_adj(s) != 1:
                    continue
                sd = side_of(s)
                name = walls.get(sd, f"wall{sd}") if sd else walls.get("other", "wall-other")
                self._wall_groups.add(name)
                assigned[s] = name
                groups2.setdefault(name, []).append(s)
                if not sd:
                    outer_other.append(s)
        # physical groups
        gid = 1
        vgroups: dict[str, list] = {}
        for t, it in sorted(kept.items()):
            vgroups.setdefault(it.group, []).append(t)
        self.groups = {"volumes": {}, "surfaces": {}}
        eps = {it.group: it.eps_r for it in vol_items}
        for name, ts in vgroups.items():
            g.model.addPhysicalGroup(3, sorted(ts), gid, name=name)
            self.groups["volumes"][name] = {"tag": gid, "entities": sorted(ts), "eps_r": eps.get(name, 1.0),
                                            "structure": any(kept[t].structure for t in ts)}
            gid += 1
        sid = 101
        for name, ss in groups2.items():
            ss = sorted(set(ss))
            if not ss:
                continue
            g.model.addPhysicalGroup(2, ss, sid, name=name)
            self.groups["surfaces"][name] = {"tag": sid, "entities": ss}
            sid += 1
        self._outer = {s for s in assigned if n_kept_adj(s) == 1}
        self._side_of = side_of
        self.built = True
        return {k: {n: v["tag"] for n, v in d.items()} for k, d in self.groups.items()}

    # ------------------------------------------------------------------ mesh
    def mesh(self, h_max: float, f_ghz: float | None = None, ppw: float = 10.0, refine: dict | None = None,
             order: int = 1, algorithm3d: int = 10, optimize: bool = True, threads: int = 4) -> dict:
        """Size fields: per volume min(h_max, lambda_material / ppw at f_ghz); near the curves of the surface
        groups in refine={'groups': [...], 'h_min':..., 'dist_min':..., 'dist_max':...} a Threshold field from
        h_min to the background. order 2 makes curved (second-order) geometry, useful for coaxial/curved parts."""
        if not self.built:
            raise RuntimeError("call build() first")
        g = self.g
        fld = g.model.mesh.field
        fields = []
        for name, v in self.groups["volumes"].items():
            h = h_max
            if f_ghz:
                h = min(h, wavelength(f_ghz, v["eps_r"]) / UNITS[self.units] / ppw)
            c = fld.add("Constant")
            fld.setNumber(c, "VIn", h)
            fld.setNumber(c, "VOut", 1e22)
            fld.setNumbers(c, "VolumesList", v["entities"])
            fld.setNumber(c, "IncludeBoundary", 1)
            fields.append(c)
        if refine:
            curves = sorted({abs(c) for n in refine["groups"] for s in self.groups["surfaces"][n]["entities"]
                             for _, c in g.model.getBoundary([(2, s)], oriented=False)})
            dist = fld.add("Distance")
            fld.setNumbers(dist, "CurvesList", curves)
            fld.setNumber(dist, "Sampling", 40)
            thr = fld.add("Threshold")
            fld.setNumber(thr, "InField", dist)
            fld.setNumber(thr, "SizeMin", refine["h_min"])
            fld.setNumber(thr, "SizeMax", h_max)
            fld.setNumber(thr, "DistMin", refine.get("dist_min", refine["h_min"]))
            fld.setNumber(thr, "DistMax", refine.get("dist_max", 10 * refine["h_min"]))
            fields.append(thr)
        mn = fld.add("Min")
        fld.setNumbers(mn, "FieldsList", fields)
        fld.setAsBackgroundMesh(mn)
        for k, v in (("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0),
                     ("Mesh.MeshSizeFromCurvature", 12), ("Mesh.MeshSizeMax", h_max), ("Mesh.Algorithm", 6),
                     ("Mesh.Algorithm3D", algorithm3d), ("Mesh.MshFileVersion", 2.2), ("Mesh.SaveAll", 0),
                     ("General.NumThreads", threads), ("Mesh.Optimize", 1 if optimize else 0)):
            g.option.setNumber(k, v)
        g.model.mesh.generate(3)
        if order > 1:
            g.model.mesh.setOrder(order)
        return self.report()

    def report(self) -> dict:
        return mesh_report(self.g)

    # ------------------------------------------------------------------ geometry for the linter
    def geometry(self) -> dict:
        g = self.g
        out = {"units": self.units, "L0": UNITS[self.units], "domain_bbox": self.domain_bbox,
               "volumes": {}, "surfaces": {}}
        sb = None
        for name, v in self.groups["volumes"].items():
            bb = _union_bbox(g, 3, v["entities"])
            out["volumes"][name] = {"tag": v["tag"], "bbox": bb, "eps_r": v["eps_r"],
                                    "volume": float(sum(g.model.occ.getMass(3, t) for t in v["entities"]))}
            dom = self.domain_bbox
            tol = 1e-6 * max(dom[3] - dom[0], dom[4] - dom[1], dom[5] - dom[2])
            inside = all(bb[k] > dom[k] + tol and bb[k + 3] < dom[k + 3] - tol for k in range(3))
            if v["structure"] and inside:  # a dielectric body inside the domain counts as structure
                sb = bb if sb is None else _merge(sb, bb)
        for name, s in self.groups["surfaces"].items():
            ents = s["entities"]
            bb = _union_bbox(g, 2, ents)
            info = {"tag": s["tag"], "bbox": bb, "area": float(sum(g.model.occ.getMass(2, t) for t in ents)),
                    "outer": all(t in self._outer for t in ents),
                    "sides": sorted({self._side_of(t) for t in ents if t in self._outer and self._side_of(t)})}
            info.update(_patch_frame(g, ents, bb))
            out["surfaces"][name] = info
            info["wall"] = name in self._wall_groups
            if not info["wall"]:  # conductors, sheets and ports: what radiates / carries the near field
                sb = bb if sb is None else _merge(sb, bb)
        out["structure_bbox"] = sb
        return out

    def write(self, path, extra: dict | None = None) -> dict:
        """Write <path>.msh (MSH 2.2), <path>.brep and <path>.tags.json. Returns the tags/geometry/report dict."""
        path = Path(path).with_suffix("")
        path.parent.mkdir(parents=True, exist_ok=True)
        self.g.write(str(path.with_suffix(".msh")))
        try:
            self.g.write(str(path.with_suffix(".brep")))
        except Exception:  # noqa: BLE001  (brep export is a convenience)
            pass
        info = self.geometry()
        info["mesh_file"] = path.with_suffix(".msh").name
        info["mesh"] = self.report()
        if extra:
            info.update(extra)
        path.with_suffix(".tags.json").write_text(json.dumps(info, indent=1), encoding="utf-8")
        return info

    def close(self):
        if self.g.isInitialized():
            self.g.finalize()


def _union_bbox(g, dim, ents):
    bb = None
    for t in ents:
        b = list(g.model.getBoundingBox(dim, t))
        bb = b if bb is None else _merge(bb, b)
    return [float(x) for x in bb]


def _merge(a, b):
    return [min(a[i], b[i]) for i in range(3)] + [max(a[i], b[i]) for i in range(3, 6)]


def _patch_frame(g, ents, bb) -> dict:
    """Normal, planarity and in-plane axes/lengths of a surface group (Palace measures lumped elements on such a
    box). Axis-aligned flat groups get the coordinate axes; others a PCA frame of the mesh nodes or BRep points."""
    ext = [bb[k + 3] - bb[k] for k in range(3)]
    span = max(ext) or 1.0
    normals = []
    for t in ents:
        try:
            lo, hi = g.model.getParametrizationBounds(2, t)
            uv = [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2]
            n = g.model.getNormal(t, uv)
            normals.append(np.asarray(n[:3], float))
        except Exception:  # noqa: BLE001
            pass
    planar = True
    if normals:
        n0 = normals[0] / (np.linalg.norm(normals[0]) or 1)
        planar = all(abs(abs(float(np.dot(n0, n / (np.linalg.norm(n) or 1)))) - 1) < 1e-6 for n in normals)
    flat_axes = [k for k in range(3) if ext[k] < 1e-6 * span + 1e-6]
    if flat_axes:
        k = flat_axes[0]
        axes = np.eye(3)
        lengths = [0.0 if j == k else float(e) for j, e in enumerate(ext)]
        normal = axes[k].tolist()
        return {"planar": planar, "normal": normal, "axes": axes.tolist(), "lengths": lengths}
    pts = []
    for t in ents:
        for _, p in g.model.getBoundary([(2, t)], recursive=True, oriented=False):
            pts.append(g.model.getValue(0, abs(p), []))
    pts = np.array(pts, float)
    if len(pts) < 3:
        return {"planar": planar, "normal": normals[0].tolist() if normals else None}
    c = pts - pts.mean(axis=0)
    _, sv, vt = np.linalg.svd(c, full_matrices=False)
    proj = c @ vt.T
    lengths = (proj.max(axis=0) - proj.min(axis=0)).tolist()
    planar = planar and lengths[2] < 1e-6 * span
    return {"planar": bool(planar), "normal": vt[2].tolist(), "axes": vt.tolist(), "lengths": lengths}


def mesh_report(g=None) -> dict:
    """Tet count, quality (minSICN, gamma) and exact Nedelec unknowns for orders 1-3 of the current Gmsh mesh."""
    g = g or _gmsh()
    tets, conn = [], []
    for etype in (4, 11, 29, 30):  # linear and curved tetrahedra
        tags, nodes = g.model.mesh.getElementsByType(etype)
        if len(tags):
            nper = len(nodes) // len(tags)
            tets.append(tags)
            conn.append(np.asarray(nodes, np.int64).reshape(-1, nper)[:, :4])
    if not tets:
        return {"tets": 0}
    tags = np.concatenate(tets)
    c = np.concatenate(conn)
    c.sort(axis=1)
    e = np.concatenate([c[:, [i, j]] for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))])
    f = np.concatenate([c[:, list(ix)] for ix in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))])
    n_e = len(np.unique(e, axis=0))
    n_f = len(np.unique(f, axis=0))
    n_t = len(c)
    q = {}
    for name in ("minSICN", "gamma"):
        try:
            v = np.asarray(g.model.mesh.getElementQualities(tags, name))
            q[name] = {"min": float(v.min()), "p01": float(np.quantile(v, 0.01)), "mean": float(v.mean()),
                       "below_0.1": int((v < 0.1).sum())}
        except Exception:  # noqa: BLE001
            pass
    verdict = "ok"
    if q.get("minSICN", {}).get("min", 1) <= 0:
        verdict = "invalid elements (minSICN <= 0)"
    elif q.get("minSICN", {}).get("min", 1) < 0.05 or q.get("gamma", {}).get("min", 1) < 0.02:
        verdict = "poor elements: check thin layers/slivers (Palace solves, but conditioning and accuracy suffer)"
    return {"tets": int(n_t), "edges": int(n_e), "faces": int(n_f),
            "unknowns": {str(p): int(nedelec_unknowns(n_e, n_f, n_t, p)) for p in (1, 2, 3)},
            "quality": q, "verdict": verdict}


# ============================================================================================ templates

def rectangular_cavity(out: Path, a: float, b: float, d: float, h: float, units: str = "mm", verbosity: int = 1):
    """Air-filled PEC box a x b x d (x, y, z). All walls are one group 'walls'. Returns tags/geometry/report."""
    m = RFModel(Path(out).stem, units, verbosity)
    try:
        m.box("cavity", 0, 0, 0, a, b, d, kind="volume")
        m.build(walls={s: "walls" for s in SIDES})
        m.mesh(h_max=h)
        return m.write(out, extra={"template": "rectangular_cavity", "a": a, "b": b, "d": d})
    finally:
        m.close()


def parallel_plate_line(out: Path, w: float, h: float, length: float, eps_r: float, h_mesh: float,
                        units: str = "mm", verbosity: int = 1):
    """TEM parallel-plate line along x: PEC plates at z = 0 and z = h ('plates'), PMC side walls ('sides',
    left out of the config), lumped ports 'port1'/'port2' covering the whole end faces (x = 0, x = length).
    With PMC sides the TEM mode is exact: Z0 = eta0 / sqrt(eps_r) * h / w, and a full-face lumped port matches its
    uniform field exactly."""
    m = RFModel(Path(out).stem, units, verbosity)
    try:
        m.box("dielectric", 0, 0, 0, length, w, h, eps_r=eps_r)
        m.rect("port1", "x", 0.0, (0, w), (0, h))
        m.rect("port2", "x", length, (0, w), (0, h))
        m.build(walls={"-z": "plates", "+z": "plates", "-y": "sides", "+y": "sides"})
        m.mesh(h_max=h_mesh)
        return m.write(out, extra={"template": "parallel_plate_line", "w": w, "h": h, "length": length,
                                   "eps_r": eps_r})
    finally:
        m.close()


def microstrip_line(out: Path, w: float, h: float, length: float, eps_r: float, side: float, top: float,
                    h_mesh: float, h_edge: float, ports: str = "lumped", f_ghz: float | None = None,
                    units: str = "mm", verbosity: int = 1):
    """Microstrip thru line along x on a substrate of height h (ground plane = z = 0, group 'ground').

    Strip: zero-thickness sheet 'strip' at z = h, width w, full length. Air to `top` above the board, board and
    air extend `side` beyond the strip edges on both sides (groups 'sides' and 'top').
    ports='lumped': vertical sheets 'port1'/'port2' from the ground to the strip at x = 0 and x = length; the end
    walls ('ends', leave them PMC) do not touch the strip, the way a PEC end wall would.
    ports='wave': the whole end faces are wave ports 'port1'/'port2' (a port that spans the domain face holds the
    whole quasi-TEM field)."""
    m = RFModel(Path(out).stem, units, verbosity)
    try:
        y0, y1 = -w / 2 - side, w / 2 + side
        m.box("substrate", 0, y0, 0, length, y1 - y0, h, eps_r=eps_r)
        m.box("air", 0, y0, h, length, y1 - y0, top, kind="air")
        m.trace("strip", 0, -w / 2, length, w / 2, h)
        walls = {"-z": "ground", "+z": "top", "-y": "sides", "+y": "sides"}
        if ports == "lumped":
            m.rect("port1", "x", 0.0, (-w / 2, w / 2), (0, h))
            m.rect("port2", "x", length, (-w / 2, w / 2), (0, h))
            walls.update({"-x": "ends", "+x": "ends"})
        elif ports == "wave":
            walls.update({"-x": "port1", "+x": "port2"})
        else:
            raise ValueError("ports must be 'lumped' or 'wave'")
        m.build(walls=walls)
        groups = ["strip"] + (["port1", "port2"] if ports == "lumped" else [])
        m.mesh(h_max=h_mesh, f_ghz=f_ghz, refine={"groups": groups, "h_min": h_edge, "dist_min": h_edge,
                                                  "dist_max": 4 * h})
        return m.write(out, extra={"template": "microstrip_line", "w": w, "h": h, "length": length,
                                   "eps_r": eps_r, "ports": ports})
    finally:
        m.close()
