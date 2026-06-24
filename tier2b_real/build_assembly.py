"""Tier-2(b) assembly: real BONE + real DENTIN(tooth24) + conforming PDL layer + root-form IMPLANT
(tooth23, Ti), coupled by *TIE at the socket walls.  Writes tier2b_real.inp.

Parts (cached by prep_meshes.py):
  BONE     real mandible crop (real alveolar sockets are voids)
  DENTIN   real tooth-24 solid
  PDL      0.2 mm conforming offset layer on the tooth-24 surface (inner nodes SHARED with dentin
           -> tooth/PDL interface is conforming; standard uniform-PDL idealisation)
  TI       real tooth-23 solid, titanium (the natural tooth is "extracted"; a root-analog implant
           fills the socket and is OSSEOINTEGRATED -> tied directly to bone, no PDL)
  BIOFILM  thin dysbiotic-growth collar carved from the bone crest around each neck
Couplings (*TIE, ADJUST=YES):  PDL-outer <-> bone socket(24);  implant <-> bone socket(23).
Steps: (1) dysbiotic biofilm growth;  (2) occlusal load on both crowns -> transmitted through bone.
Units mm, MPa.
"""
import json
import os
import sys
import numpy as np

FEM = "/home/nishioka/IKM_Hiwi/FEM"
sys.path.insert(0, FEM)
from biofilm_conformal_tet import compute_vertex_normals  # noqa: E402

OUT = f"{FEM}/tier2b_real"
T24_AXIS = np.array([-63.9, -41.2])
T23_AXIS = np.array([-69.4, -41.0])
CREST_Z = 29.0
PDL_THICK = 0.25
EPS_GROWTH = 0.19
# crop-box (must match prep_meshes.py)
X0, X1, Y0, Y1, Z0, Z1 = -76.0, -59.0, -47.0, -37.5, 15.0, 31.0
MATS = {"BONE": (13700., 0.30), "CORTICAL": (13700., 0.30), "CANCELLOUS": (1000., 0.30),
        "TI": (110000., 0.34), "DENTIN": (18000., 0.31), "PDL": (50., 0.45), "BIOFILM": (1.0, 0.45),
        "CROWN": (95000., 0.30),   # monolithic lithium-disilicate-class ceramic crown
        "GINGIVA": (3.0, 0.45),    # peri-implant mucosa (conformal cuff, mechanically negligible)
        "ENAMEL": (84000., 0.30),  # enamel cap on the natural neighbour-tooth clinical crown
        "CEMENTUM": (15000., 0.31),   # thin layer on the tooth root (over the dentin, under the PDL)
        "PULP": (2.0, 0.45),          # tooth pulp (root canal + chamber), embedded in dentin
        "ABUTSCREW": (110000., 0.34)} # abutment screw down the implant axis, embedded in TI
GUM_Z1 = 31.9                 # gingival margin height (gum cuff top); biofilm SULC_Z1 = 31.3
ENAM_Z0 = 31.5                # enamel covers the supragingival clinical crown (CEJ -> occlusal)
CORTICAL_THICK = 1.8          # cortical shell / lamina-dura thickness (mm) for the two-layer bone
# C3D4 face -> local node triple (Abaqus 1-based)
FACE_NODES = {1: (0, 1, 2), 2: (0, 3, 1), 3: (1, 3, 2), 4: (2, 3, 0)}
# anatomical peri-implant/peri-tooth biofilm sleeve (supracrestal, in the gingival sulcus, ON the
# hard-tissue surface): a thin annulus from the bone crest up into the sulcus, tied to the implant /
# tooth neck.  Replaces the old subcrestal bone-carved collar for the coupled (crown) assembly.
SULC_Z0, SULC_Z1 = CREST_Z + 0.0, CREST_Z + 2.3    # crest -> ~2.3 mm up the sulcus
# 6-tet (Kuhn) split of a hex with corners 0..7 (0-3 bottom ring, 4-7 top ring)
_HEX6 = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]


def biofilm_sleeve(ax, r_in, r_out, z0=SULC_Z0, z1=SULC_Z1, nth=56, nz=8):
    """Thin annular biofilm sleeve (body of revolution) hugging a neck in the sulcus. r_in / r_out may
    be scalars OR per-z-level arrays (length nz+1) so the sleeve FOLLOWS a flaring profile (e.g. a tooth
    crown widening coronally) and stays in contact. Inner shell = r_in (tied to the neck)."""
    th = np.linspace(0, 2 * np.pi, nth, endpoint=False)
    zs = np.linspace(z0, z1, nz + 1)
    rin = np.full(nz + 1, r_in, float) if np.isscalar(r_in) else np.asarray(r_in, float)
    rout = np.full(nz + 1, r_out, float) if np.isscalar(r_out) else np.asarray(r_out, float)
    nd = np.empty(((nz + 1) * nth * 2, 3))

    def nid(it, jz, sh):
        return (it * (nz + 1) + jz) * 2 + sh
    for it in range(nth):
        for jz in range(nz + 1):
            for sh, rr in enumerate((rin[jz], rout[jz])):
                nd[nid(it, jz, sh)] = (ax[0] + rr * np.cos(th[it]), ax[1] + rr * np.sin(th[it]), zs[jz])
    tets = []
    for it in range(nth):
        j2 = (it + 1) % nth
        for jz in range(nz):
            n = [nid(it, jz, 0), nid(j2, jz, 0), nid(j2, jz, 1), nid(it, jz, 1),
                 nid(it, jz + 1, 0), nid(j2, jz + 1, 0), nid(j2, jz + 1, 1), nid(it, jz + 1, 1)]
            for a, b, c, d in _HEX6:
                tets.append([n[a], n[b], n[c], n[d]])
    return nd, np.array(tets, dtype=np.int64)


def offset_shell(bn, bt, axis, z0=SULC_Z0, z1=SULC_Z1, r0=0.0, r1=0.15):
    """CONFORMAL layer on a neck: the body's lateral surface in [z0,z1], radially offset to [r0,r1], as
    a thin prism shell that HUGS the (irregular) surface exactly -- the same construction as the PDL, so
    an irregular natural tooth gets a perfectly-fitting layer (no circular over/under-shoot). r0=0 puts
    the inner shell ON the surface (biofilm); r0>0 stacks an outer layer (gingiva over the biofilm).
    Returns (nodes, tets, n_inner); nodes[:n_inner] = inner shell."""
    ff = free_faces(bt, np.arange(len(bt)))     # free_faces is module-level (resolved at call time)
    faces = []
    for key in ff:
        tri = list(key); p = bn[tri]; c = p.mean(0)
        if not (z0 <= c[2] <= z1):
            continue
        n = np.cross(p[1] - p[0], p[2] - p[0]); nn = np.linalg.norm(n)
        if nn == 0 or abs(n[2] / nn) > 0.6:        # keep the lateral wall, drop top/bottom caps
            continue
        faces.append(tri)
    if not faces:
        return np.zeros((0, 3)), np.zeros((0, 4), np.int64), 0
    used = np.unique(np.array(faces).ravel())
    remap = -np.ones(len(bn), np.int64); remap[used] = np.arange(len(used))
    surf = bn[used].astype(float)
    rxy = surf[:, :2] - axis; rl = np.hypot(rxy[:, 0], rxy[:, 1]); rl[rl == 0] = 1.0
    u = rxy / rl[:, None]
    inner = surf.copy(); inner[:, 0] += r0 * u[:, 0]; inner[:, 1] += r0 * u[:, 1]
    outer = surf.copy(); outer[:, 0] += r1 * u[:, 0]; outer[:, 1] += r1 * u[:, 1]
    N = len(used); nodes = np.vstack([inner, outer]); tets = []
    for tri in faces:
        a, b, c = (int(remap[t]) for t in tri); a2, b2, c2 = a + N, b + N, c + N
        tets += [[a, b, c, c2], [a, b, c2, b2], [a, b2, c2, a2]]   # triangular prism -> 3 tets
    return nodes, np.array(tets, np.int64), N


def solid_revolve(zs, radii, centers, nth=20):
    """Solid body of revolution: at each z-level a centre + a ring of radius `radii[k]` about
    `centers[k]` (a 2-D point, so the axis may bend to follow a tooth). Cone/cylinder tet fill.
    Used for the pulp (along the dentin centre-line), the abutment screw, and the papilla."""
    th = np.linspace(0, 2 * np.pi, nth, endpoint=False); nz = len(zs)
    nd = []
    def cid(k): return k * (nth + 1)
    def rid(k, i): return k * (nth + 1) + 1 + i
    for k in range(nz):
        nd.append([centers[k][0], centers[k][1], zs[k]])
        for i in range(nth):
            nd.append([centers[k][0] + radii[k] * np.cos(th[i]),
                       centers[k][1] + radii[k] * np.sin(th[i]), zs[k]])
    tets = []
    for k in range(nz - 1):
        for i in range(nth):
            j = (i + 1) % nth
            c0, c1 = cid(k), cid(k + 1); a, b = rid(k, i), rid(k, j); a2, b2 = rid(k + 1, i), rid(k + 1, j)
            tets += [[c0, a, b, c1], [a, b, c1, b2], [a, c1, b2, a2]]
    return np.array(nd, float), np.array(tets, np.int64)


def dentin_centerline(npts, tets, zs):
    """Centre (x,y) of a body at each z-level -> follows a tilted natural-tooth axis."""
    ec = npts[tets].mean(1); out = []
    for z in zs:
        m = np.abs(ec[:, 2] - z) <= (zs[1] - zs[0]) * 0.9
        out.append(ec[m, :2].mean(0) if m.any() else (out[-1] if out else npts[:, :2].mean(0)))
    return out


def clean_tets(conn, nodes, vol_min=2.5e-4):
    """Flip negative-volume tets (winding) and drop sliver tets (volume < vol_min mm^3, which
    Abaqus rejects as zero/small). conn 0-based into nodes."""
    p = nodes[conn]
    v6 = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0])
    out = conn.copy()
    neg = v6 < 0
    out[neg] = out[neg][:, [0, 1, 3, 2]]          # swap last two -> flip orientation
    keep = np.abs(v6) / 6.0 > vol_min
    return out[keep], int((~keep).sum()), int(neg.sum())


def boundary_face_centroids(tets, nodes, exclude_planes):
    """Centroids of the free (boundary) faces of a tet block, excluding faces lying on the artificial
    crop planes (so 'outer surface' = real mandible cortex + socket lamina dura, not the cut faces)."""
    cnt = {}
    for t in tets:
        for (a, b, c) in FACE_NODES.values():
            key = tuple(sorted((int(t[a]), int(t[b]), int(t[c]))))
            cnt[key] = cnt.get(key, 0) + 1
    cents = []
    for key, n in cnt.items():
        if n != 1:
            continue
        fc = nodes[list(key)].mean(axis=0)
        on_plane = any(abs(fc[ax] - val) < 0.8 for ax, val in exclude_planes)
        if not on_plane:
            cents.append(fc)
    return np.array(cents)


def free_faces(tets, elem_gids):
    """Return dict sorted-face -> (gid, face_id) for faces that appear exactly once."""
    seen = {}
    for loc, (t, g) in enumerate(zip(tets, elem_gids)):
        for fid, (a, b, c) in FACE_NODES.items():
            key = tuple(sorted((int(t[a]), int(t[b]), int(t[c]))))
            if key in seen:
                seen[key] = None            # shared -> not free
            else:
                seen[key] = (g, fid)
    return {k: v for k, v in seen.items() if v is not None}


def main():
    # optional argv: <implant_cache.npz> <job_basename>  (defaults = original root-analog model)
    imp_cache = sys.argv[1] if len(sys.argv) > 1 else "cache_implant.npz"
    job = sys.argv[2] if len(sys.argv) > 2 else "tier2b_real"
    WITH_CROWN = "crown" in job          # ceramic crown seated on the abutment -> occlusal moment arm
    ANAT_BIO = job != "tier2b_real"       # anatomical biofilm: supracrestal sulcular sleeve on the neck
    #                                       (all leveled-up coupled jobs); only original tier2b_real
    #                                       keeps the legacy subcrestal bone-carved collar.
    if len(sys.argv) > 3:                 # optional crown Young's modulus override (MPa), design sweep
        MATS["CROWN"] = (float(sys.argv[3]), 0.30)
    # MATS_OVERRIDE env var: JSON dict {mat_name: [E_MPa, nu]} for sensitivity sweeps. Keeps the headline
    # workflow unchanged when unset.
    _mats_override = os.environ.get("MATS_OVERRIDE")
    if _mats_override:
        for k, v in json.loads(_mats_override).items():
            MATS[k] = (float(v[0]), float(v[1]))
        print("[build_assembly] MATS_OVERRIDE applied:", _mats_override)
    bone = np.load(f"{OUT}/cache_bone.npz")
    dent = np.load(f"{OUT}/cache_dentin.npz")
    imp = np.load(f"{OUT}/{imp_cache}")
    bn, bt = bone["nodes"], bone["tets"]
    dn, dt = dent["nodes"], dent["tets"]
    sverts, sfaces = dent["sverts"], dent["sfaces"]
    inn, it = imp["nodes"], imp["tets"]
    if WITH_CROWN:
        crown_cache = "cache_crown_real.npz" if "crownreal" in job else "cache_crown.npz"
        crw = np.load(f"{OUT}/{crown_cache}")
        cnn, ct = crw["nodes"], crw["tets"]

    # PDL: offset tooth-24 surface outward along vertex normals
    vn = compute_vertex_normals(sverts, sfaces)
    pdl_outer_xyz = sverts + PDL_THICK * vn

    # map surface verts -> dentin node indices (exact coords)
    key2idx = {tuple(np.round(p, 5)): i for i, p in enumerate(dn)}
    svert2dn = np.array([key2idx[tuple(np.round(p, 5))] for p in sverts], dtype=np.int64)
    V = len(sverts)

    # ---- global node numbering ----
    Nb, Nd, Ni = len(bn), len(dn), len(inn)
    off_d = Nb
    off_i = Nb + Nd
    off_p = Nb + Nd + Ni                 # pdl OUTER nodes (inner reuse dentin ids)
    off_c = off_p + V                    # crown nodes (after pdl-outer); 0 of them if no crown
    Nc = len(cnn) if WITH_CROWN else 0
    off_bf = off_c + Nc                  # anatomical-biofilm shell nodes (after crown)
    # CONFORMAL sulcular biofilm: a thin radial-offset shell of each neck surface (PDL-style), so it
    # hugs the irregular natural tooth exactly -- no circular sleeve over/under-shoot.
    if ANAT_BIO:
        def neck_ctr(npts, tets):
            ec = npts[tets].mean(1); s0 = (ec[:, 2] >= SULC_Z0 - 0.3) & (ec[:, 2] <= SULC_Z0 + 0.8)
            return ec[s0, :2].mean(0)
        ic, tc = neck_ctr(inn, it), neck_ctr(dn, dt)
        # biofilm = inner film ON the surface (r 0->0.15); gingiva = conformal cuff STACKED outside it
        # (r 0.18->1.1, up to the higher gum margin) -- both hug the irregular tooth exactly.
        sb_i, st_i, ni_in = offset_shell(inn, it, ic, r0=0.0, r1=0.15)
        sb_t, st_t, nt_in = offset_shell(dn, dt, tc, r0=0.0, r1=0.15)
        Nbi = len(sb_i)
        bf_nodes = np.vstack([sb_i, sb_t]); bf_tets = np.vstack([st_i + off_bf, st_t + off_bf + Nbi])
        bf_inner_imp = set(range(off_bf, off_bf + ni_in))
        bf_inner_too = set(range(off_bf + Nbi, off_bf + Nbi + nt_in))
        off_gum = off_bf + len(bf_nodes)
        gi_i, gtt_i, gii = offset_shell(inn, it, ic, z1=GUM_Z1, r0=0.18, r1=1.1)
        gi_t, gtt_t, gti = offset_shell(dn, dt, tc, z1=GUM_Z1, r0=0.18, r1=1.1)
        Ngi = len(gi_i)
        gum_nodes = np.vstack([gi_i, gi_t]); gum_tets = np.vstack([gtt_i + off_gum, gtt_t + off_gum + Ngi])
        gum_inner_imp = set(range(off_gum, off_gum + gii))
        gum_inner_too = set(range(off_gum + Ngi, off_gum + Ngi + gti))
        # enamel cap on the natural neighbour-tooth clinical crown (CEJ -> occlusal), a conformal shell
        off_en = off_gum + len(gum_nodes)
        enam_nodes, en_t, en_in = offset_shell(dn, dt, tc, z0=ENAM_Z0, z1=38.6, r0=0.0, r1=0.45)
        enam_tets = en_t + off_en
        enam_inner = set(range(off_en, off_en + en_in))
        # cementum: a very thin layer on the dentin root (subgingival, between dentin and PDL)
        off_cem = off_en + len(enam_nodes)
        cem_nodes, cm_t, cm_in = offset_shell(dn, dt, tc, z0=19.0, z1=ENAM_Z0, r0=0.0, r1=0.08)
        cem_tets = cm_t + off_cem; cem_inner = set(range(off_cem, off_cem + cm_in))
        # pulp: solid body of revolution along the dentin centre-line, embedded in the dentin
        off_pulp = off_cem + len(cem_nodes)
        zp = np.linspace(20.5, 36.0, 10); cp = dentin_centerline(dn, dt, zp)
        rp = np.array([0.16, 0.20, 0.24, 0.30, 0.42, 0.58, 0.66, 0.55, 0.30, 0.08])
        pulp_nodes, pl_t = solid_revolve(zp, rp, cp); pulp_tets = pl_t + off_pulp
        # abutment screw: solid down the implant axis, embedded in the Ti
        off_scr = off_pulp + len(pulp_nodes)
        zss = np.linspace(25.0, 32.4, 8); css = [ic] * len(zss); rss = np.full(len(zss), 0.8)
        scr_nodes, sc_t = solid_revolve(zss, rss, css); scr_tets = sc_t + off_scr
        # interdental papilla: a gum dome between the necks (GINGIVA mat); base on the bone crest,
        # tied to it so the connected dome is anchored.
        off_pap = off_scr + len(scr_nodes)
        mid = (np.asarray(ic) + np.asarray(tc)) / 2.0
        zpa = np.linspace(CREST_Z, GUM_Z1 + 0.4, 6); cpa = [mid] * len(zpa)
        rpa = np.array([1.25, 1.35, 1.20, 0.95, 0.55, 0.10])
        pap_nodes, pp_t = solid_revolve(zpa, rpa, cpa); pap_tets = pp_t + off_pap
        print(f"conformal shells: biofilm {len(bf_tets)} / gingiva {len(gum_tets)} / enamel {len(enam_tets)} "
              f"/ cementum {len(cem_tets)} / pulp {len(pulp_tets)} / screw {len(scr_tets)} / papilla {len(pap_tets)}")
    else:
        z0 = np.zeros((0, 3)); z4 = np.zeros((0, 4), np.int64)
        bf_nodes = gum_nodes = enam_nodes = cem_nodes = pulp_nodes = scr_nodes = pap_nodes = z0
        bf_tets = gum_tets = enam_tets = cem_tets = pulp_tets = scr_tets = pap_tets = z4
        bf_inner_imp = bf_inner_too = gum_inner_imp = gum_inner_too = enam_inner = cem_inner = set()

    parts = [bn, dn, inn, pdl_outer_xyz]
    if WITH_CROWN:
        parts.append(cnn)
    if ANAT_BIO:
        parts += [bf_nodes, gum_nodes, enam_nodes, cem_nodes, pulp_nodes, scr_nodes, pap_nodes]
    nodes = np.vstack(parts)

    # ---- element connectivity (global, 0-based) ----
    bt_g = bt
    dt_g = dt + off_d
    it_g = it + off_i
    ct_g = (ct + off_c) if WITH_CROWN else None
    # pdl tets: inner idx<V -> dentin global; outer idx>=V -> pdl-outer global
    pdl_tris = sfaces
    pdl_tets = []
    for (i0, i1, i2) in pdl_tris:
        B0, B1, B2 = off_d + svert2dn[i0], off_d + svert2dn[i1], off_d + svert2dn[i2]
        T0, T1, T2 = off_p + i0, off_p + i1, off_p + i2
        pdl_tets += [[B0, B1, B2, T2], [B0, B1, T2, T1], [B0, T1, T2, T0]]
    pdl_tets = np.array(pdl_tets, dtype=np.int64)

    # ---- biofilm collar: carve from BONE crest tets near each neck ----
    bcent = bn[bt].mean(axis=1)
    r24 = np.hypot(bcent[:, 0] - T24_AXIS[0], bcent[:, 1] - T24_AXIS[1])
    r23 = np.hypot(bcent[:, 0] - T23_AXIS[0], bcent[:, 1] - T23_AXIS[1])
    band = (bcent[:, 2] >= CREST_Z - 2.0) & (bcent[:, 2] <= CREST_Z + 1.5)
    coll = band & ((r24 <= 2.6) | (r23 <= 2.6))
    if ANAT_BIO:
        coll = np.zeros(len(bt), bool)   # bone stays intact; biofilm = supracrestal sulcular sleeves
    print(f"biofilm collar tets (bone-carved): {coll.sum()}  ANAT_BIO={ANAT_BIO}")

    # ---- two-layer bone (level-up): cortical shell / lamina dura vs cancellous core ----
    # only for the generic-implant job; tier2b_real stays single-layer (preserved original).
    bone_block = bt_g[~coll]
    biofilm_block = bf_tets if ANAT_BIO else bt_g[coll]   # sulcular sleeve vs old bone-carved collar
    if job != "tier2b_real":
        from scipy.spatial import cKDTree
        planes = [(0, X0), (0, X1), (1, Y0), (1, Y1), (2, Z0)]   # artificial crop faces to ignore
        surf_c = boundary_face_centroids(bt_g, bn, planes)       # real outer + socket-wall faces
        bc = bn[bone_block].mean(axis=1)
        dist, _ = cKDTree(surf_c).query(bc)
        is_cort = dist < CORTICAL_THICK
        print(f"two-layer bone: cortical={int(is_cort.sum())} cancellous={int((~is_cort).sum())}")
        blocks = [("CORTICAL", bone_block[is_cort]), ("CANCELLOUS", bone_block[~is_cort]),
                  ("BIOFILM", biofilm_block), ("DENTIN", dt_g), ("TI", it_g), ("PDL", pdl_tets)]
    else:
        blocks = [("BONE", bone_block), ("BIOFILM", biofilm_block),
                  ("DENTIN", dt_g), ("TI", it_g), ("PDL", pdl_tets)]
    if WITH_CROWN:
        blocks.append(("CROWN", ct_g))
    if ANAT_BIO and len(gum_tets):
        blocks.append(("GINGIVA", gum_tets))
    if ANAT_BIO and len(enam_tets):
        blocks += [("ENAMEL", enam_tets), ("CEMENTUM", cem_tets), ("PULP", pulp_tets),
                   ("ABUTSCREW", scr_tets), ("PAPILLA", pap_tets)]
    eid = {}; gid = 1; elem_rows = []
    for name, conn in blocks:
        conn, ndrop, nflip = clean_tets(conn, nodes)
        if ndrop or nflip:
            print(f"  clean {name}: flipped {nflip}, dropped {ndrop} degenerate")
        ids = np.arange(gid, gid + len(conn))
        eid[name] = (ids, conn); gid += len(conn)
        elem_rows.append((name, ids, conn))
    print("element counts:", {n: len(eid[n][1]) for n in eid})

    # ---- TIE surfaces ----
    # slave: PDL exterior faces (outer + rim) ; implant exterior faces
    pdl_free_all = free_faces(pdl_tets, eid["PDL"][0])
    # keep ONLY the outer (bone-facing) faces: all 3 nodes are PDL-outer nodes (>= off_p).
    # the inner faces sit on the tooth surface and must NOT be tied to bone (that would kill the
    # 0.25 mm PDL compliance), so they are dropped.
    pdl_free = {k: v for k, v in pdl_free_all.items() if min(k) >= off_p}
    imp_free = free_faces(it_g, eid["TI"][0])
    # master: BONE boundary faces near either socket wall (over all bone elsets: single or two-layer).
    # The old bone-carved BIOFILM collar was part of the socket; the anatomical sulcular sleeve is NOT.
    _bset = ("BONE", "CORTICAL", "CANCELLOUS") if ANAT_BIO else ("BONE", "CORTICAL", "CANCELLOUS", "BIOFILM")
    bone_sets = [s for s in _bset if s in eid]
    bone_all = np.vstack([eid[s][1] for s in bone_sets])
    bone_gids = np.concatenate([eid[s][0] for s in bone_sets])
    bone_free = free_faces(bone_all, bone_gids)
    # restrict master to socket region (near tooth axes, below crest). The tooth-23 alveolar socket
    # is wider (buccolingually) than a standard Ø4.1 implant, so the tooth-23 master radius is enlarged
    # and the implant tie uses a larger position tolerance (the standard screw is bonded to the socket
    # walls -- an explicit idealisation, since a real placement would be in a healed/drilled ridge).
    R23_MASTER = 4.8 if job != "tier2b_real" else 3.2
    master = {}
    for k, v in bone_free.items():
        fc = nodes[list(k)].mean(axis=0)
        rr24 = np.hypot(fc[0] - T24_AXIS[0], fc[1] - T24_AXIS[1])
        rr23 = np.hypot(fc[0] - T23_AXIS[0], fc[1] - T23_AXIS[1])
        if fc[2] <= CREST_Z + 0.5 and (rr24 <= 3.2 or rr23 <= R23_MASTER):
            master[k] = v
    print(f"TIE faces: master(bone socket)={len(master)} slave_pdl={len(pdl_free)} slave_imp={len(imp_free)}")

    # ---- crown <-> abutment TIE ----
    # The ceramic crown is seated on the abutment top: tie the crown intaglio (seat) faces to the
    # abutment-top faces. Crown = slave (finer, more compliant), abutment top = master (stiff Ti).
    # The crown is a hollow cap sheathing the abutment: the load-transfer seat is the bore-roof disk
    # at the abutment-top plane (z_abut, r<RBORE), NOT the crown's lowest faces (the cervical margin).
    crown_seat, abut_top = {}, {}
    if WITH_CROWN:
        z_abut = nodes[off_i:off_i + Ni, 2].max()              # abutment platform top z (= 32.5)
        crown_free = free_faces(ct_g, eid["CROWN"][0])
        for k, v in crown_free.items():
            fc = nodes[list(k)].mean(axis=0)
            r = np.hypot(fc[0] - T23_AXIS[0], fc[1] - T23_AXIS[1])
            if abs(fc[2] - z_abut) < 0.4 and r < 2.0:          # bore-roof seat disk only
                crown_seat[k] = v
        for k, v in imp_free.items():
            if nodes[list(k)][:, 2].min() >= z_abut - 0.4:     # abutment-top disk faces only
                abut_top[k] = v
        print(f"crown TIE: slave_seat={len(crown_seat)} master_abut_top={len(abut_top)} "
              f"(z_abut={z_abut:.2f})")

    # ---- anatomical biofilm sleeve TIE (sulcular biofilm adheres to the implant / tooth neck) ----
    # The sleeve is a free body; its INNER shell faces are tied to the neck surface so it is anchored
    # (and physically adherent). Implant sleeve -> IMP_OUT; tooth sleeve -> the dentin neck faces.
    bio_imp_inner, bio_too_inner, dent_neck = {}, {}, {}
    if ANAT_BIO and "BIOFILM" in eid:
        # the conformal shell's INNER faces are exactly on the body surface -> identify them by node
        # membership (all 3 nodes are inner-shell nodes), then tie them to the implant / tooth surface.
        bf_free = free_faces(eid["BIOFILM"][1], eid["BIOFILM"][0])
        for k, v in bf_free.items():
            ks = set(k)
            if ks <= bf_inner_imp:
                bio_imp_inner[k] = v
            elif ks <= bf_inner_too:
                bio_too_inner[k] = v
        dt_free = free_faces(dt_g, eid["DENTIN"][0])
        for k, v in dt_free.items():
            fc = nodes[list(k)].mean(axis=0)
            if SULC_Z0 - 0.4 <= fc[2] <= GUM_Z1 + 0.4 and np.hypot(fc[0] - tc[0], fc[1] - tc[1]) < 4.5:
                dent_neck[k] = v
        print(f"biofilm conformal TIE: imp_inner={len(bio_imp_inner)} too_inner={len(bio_too_inner)} "
              f"dent_neck(master)={len(dent_neck)}")

    # ---- gingiva conformal cuff TIE (mucosa adheres to the neck, outside the biofilm) ----
    gum_imp_inner, gum_too_inner = {}, {}
    if ANAT_BIO and "GINGIVA" in eid:
        for k, v in free_faces(eid["GINGIVA"][1], eid["GINGIVA"][0]).items():
            ks = set(k)
            if ks <= gum_inner_imp:
                gum_imp_inner[k] = v
            elif ks <= gum_inner_too:
                gum_too_inner[k] = v
        print(f"gingiva conformal TIE: imp_inner={len(gum_imp_inner)} too_inner={len(gum_too_inner)}")

    # ---- enamel cap TIE (enamel adheres to the dentin clinical crown) ----
    enam_inner_f, dent_crown = {}, {}
    if ANAT_BIO and "ENAMEL" in eid:
        for k, v in free_faces(eid["ENAMEL"][1], eid["ENAMEL"][0]).items():
            if set(k) <= enam_inner:
                enam_inner_f[k] = v
        for k, v in free_faces(dt_g, eid["DENTIN"][0]).items():
            fc = nodes[list(k)].mean(axis=0)
            if ENAM_Z0 - 0.5 <= fc[2] <= 38.6 and np.hypot(fc[0] - tc[0], fc[1] - tc[1]) < 4.6:
                dent_crown[k] = v
        print(f"enamel TIE: inner={len(enam_inner_f)} dent_crown(master)={len(dent_crown)}")

    # ---- cementum TIE (thin root layer adheres to the dentin root) ----
    cem_inner_f, dent_root = {}, {}
    if ANAT_BIO and "CEMENTUM" in eid:
        for k, v in free_faces(eid["CEMENTUM"][1], eid["CEMENTUM"][0]).items():
            if set(k) <= cem_inner:
                cem_inner_f[k] = v
        for k, v in free_faces(dt_g, eid["DENTIN"][0]).items():
            fc = nodes[list(k)].mean(axis=0)
            if 18.5 <= fc[2] <= ENAM_Z0 + 0.3 and np.hypot(fc[0] - tc[0], fc[1] - tc[1]) < 4.0:
                dent_root[k] = v
        print(f"cementum TIE: inner={len(cem_inner_f)} dent_root(master)={len(dent_root)}")

    # ---- papilla base TIE (interdental gum dome anchored to the bone crest) ----
    pap_base, bone_crest = {}, {}
    if ANAT_BIO and "PAPILLA" in eid:
        for k, v in free_faces(eid["PAPILLA"][1], eid["PAPILLA"][0]).items():
            if nodes[list(k)].mean(axis=0)[2] < CREST_Z + 0.4:
                pap_base[k] = v
        for k, v in bone_free.items():
            fc = nodes[list(k)].mean(axis=0)
            if abs(fc[2] - CREST_Z) < 0.9 and np.hypot(fc[0] - mid[0], fc[1] - mid[1]) < 1.9:
                bone_crest[k] = v
        print(f"papilla TIE: base={len(pap_base)} bone_crest(master)={len(bone_crest)}")

    # ---- node sets ----
    nx, ny, nz = nodes[:, 0], nodes[:, 1], nodes[:, 2]
    tol = 0.6
    fixed = np.where((np.arange(len(nodes)) < Nb) &
                     ((nz <= Z0 + tol) | (nx <= X0 + tol) | (nx >= X1 - tol) |
                      (ny <= Y0 + tol) | (ny >= Y1 - tol)))[0]
    # crown tops (load): dentin nodes & implant nodes near their max z
    d_ids = np.arange(off_d, off_d + Nd)
    i_ids = np.arange(off_i, off_i + Ni)
    too_top = d_ids[nz[d_ids] >= nz[d_ids].max() - 1.5]
    imp_top = i_ids[nz[i_ids] >= nz[i_ids].max() - 1.5]
    # with a crown, the occlusal load is applied at the crown occlusal table (z~40.8), NOT the
    # abutment top (z=32.5) -> the bite force now acts through the crown-height moment arm.
    if WITH_CROWN:
        c_ids = np.arange(off_c, off_c + len(cnn))
        crown_top = c_ids[nz[c_ids] >= nz[c_ids].max() - 1.0]
        load_top = crown_top
    else:
        load_top = imp_top
    print(f"fixed={len(fixed)} tooth_top={len(too_top)} imp_top={len(imp_top)} "
          f"load_top={len(load_top)} (crown={WITH_CROWN})")

    # ---- write INP ----
    L = []; ap = L.append
    ap("*HEADING"); ap(" Tier-2(b) real-shape implant + tooth24 + PDL + mandible, TIE-coupled")
    ap("*NODE")
    for i, (x, y, z) in enumerate(nodes, start=1):
        ap(" %d, %.5f, %.5f, %.5f" % (i, x, y, z))
    for name, ids, conn in elem_rows:
        ap("*ELEMENT, TYPE=C3D4, ELSET=%s" % name)
        for e, c in zip(ids, conn + 1):
            ap(" %d, %d, %d, %d, %d" % (e, c[0], c[1], c[2], c[3]))
    for m, (E, nu) in MATS.items():
        if m not in eid:                      # only emit materials that have elements
            continue
        ap("*SOLID SECTION, ELSET=%s, MATERIAL=%s" % (m, m))
        ap("*MATERIAL, NAME=%s" % m)
        ap("*ELASTIC"); ap(" %.1f, %.3f" % (E, nu))
        if m == "BIOFILM":
            ap("*EXPANSION"); ap(" %.6f" % EPS_GROWTH)
    if "PAPILLA" in eid:                          # interdental papilla uses the gingiva material
        ap("*SOLID SECTION, ELSET=PAPILLA, MATERIAL=GINGIVA")
    # embedded internal bodies (no overlap meshing needed): pulp in dentin, abutment screw in Ti,
    # papilla base in bone -- each embedded element's response is constrained to its host.
    if ANAT_BIO and "PULP" in eid:
        ap("*EMBEDDED ELEMENT, HOST ELSET=DENTIN"); ap(" PULP")
    if ANAT_BIO and "ABUTSCREW" in eid:
        ap("*EMBEDDED ELEMENT, HOST ELSET=TI"); ap(" ABUTSCREW")

    def surf(name, faces):
        ap("*SURFACE, NAME=%s, TYPE=ELEMENT" % name)
        for (g, fid) in faces.values():
            ap(" %d, S%d" % (g, fid))
    surf("BONE_SOCKET", master)
    surf("PDL_OUT", pdl_free)
    surf("IMP_OUT", imp_free)
    # ADJUST=NO: do NOT move slave nodes onto the master (that collapses the thin PDL / implant-tip
    # tets); the small initial socket gap is absorbed into the rigid tie constraint instead.
    imp_tol = 2.8 if job != "tier2b_real" else 1.0   # generic screw < natural socket -> larger gap
    ap("*TIE, NAME=T_PDL, ADJUST=NO, POSITION TOLERANCE=1.0")
    ap(" PDL_OUT, BONE_SOCKET")
    ap("*TIE, NAME=T_IMP, ADJUST=NO, POSITION TOLERANCE=%.1f" % imp_tol)
    ap(" IMP_OUT, BONE_SOCKET")
    if WITH_CROWN:
        surf("CROWN_SEAT", crown_seat)
        surf("ABUT_TOP", abut_top)
        ap("*TIE, NAME=T_CROWN, ADJUST=NO, POSITION TOLERANCE=0.5")
        ap(" CROWN_SEAT, ABUT_TOP")
    if ANAT_BIO and bio_imp_inner:
        surf("BIO_IMP_IN", bio_imp_inner)
        # ADJUST=NO: do not pull the thin (0.3 mm) sleeve inner nodes onto the neck (that flattens the
        # sleeve tets); the small gap is absorbed into the tie. Biofilm adheres to the Ti neck.
        ap("*TIE, NAME=T_BIO_IMP, ADJUST=NO, POSITION TOLERANCE=0.6")
        ap(" BIO_IMP_IN, IMP_OUT")
    if ANAT_BIO and bio_too_inner and dent_neck:
        surf("BIO_TOO_IN", bio_too_inner)
        surf("DENT_NECK", dent_neck)
        ap("*TIE, NAME=T_BIO_TOO, ADJUST=NO, POSITION TOLERANCE=0.6")    # biofilm adheres to the tooth neck
        ap(" BIO_TOO_IN, DENT_NECK")
    if ANAT_BIO and gum_imp_inner:                                      # gingiva cuff adheres to the necks
        surf("GUM_IMP_IN", gum_imp_inner)
        ap("*TIE, NAME=T_GUM_IMP, ADJUST=NO, POSITION TOLERANCE=1.4")
        ap(" GUM_IMP_IN, IMP_OUT")
    if ANAT_BIO and gum_too_inner and dent_neck:
        surf("GUM_TOO_IN", gum_too_inner)
        ap("*TIE, NAME=T_GUM_TOO, ADJUST=NO, POSITION TOLERANCE=1.4")
        ap(" GUM_TOO_IN, DENT_NECK")
    if ANAT_BIO and enam_inner_f and dent_crown:                       # enamel adheres to the dentin crown
        surf("ENAM_IN", enam_inner_f); surf("DENT_CROWN", dent_crown)
        ap("*TIE, NAME=T_ENAM, ADJUST=NO, POSITION TOLERANCE=0.6")
        ap(" ENAM_IN, DENT_CROWN")
    if ANAT_BIO and cem_inner_f and dent_root:                         # cementum adheres to the dentin root
        surf("CEM_IN", cem_inner_f); surf("DENT_ROOT", dent_root)
        ap("*TIE, NAME=T_CEM, ADJUST=NO, POSITION TOLERANCE=0.6")
        ap(" CEM_IN, DENT_ROOT")
    if ANAT_BIO and pap_base and bone_crest:                           # papilla anchored to the bone crest
        surf("PAP_BASE", pap_base); surf("BONE_CREST", bone_crest)
        ap("*TIE, NAME=T_PAP, ADJUST=NO, POSITION TOLERANCE=0.8")
        ap(" PAP_BASE, BONE_CREST")

    def nset(nm, ids):
        ids = np.unique(np.asarray(ids)) + 1
        ap("*NSET, NSET=%s" % nm)
        for k in range(0, len(ids), 16):
            ap(" " + ",".join(str(int(v)) for v in ids[k:k + 16]))
    nset("FIXED", fixed); nset("TOOTOP", too_top); nset("IMPTOP", imp_top)
    nset("ALLN", np.arange(len(nodes)))

    ap("*INITIAL CONDITIONS, TYPE=TEMPERATURE"); ap(" ALLN, 0.0")
    ap("*BOUNDARY"); ap(" FIXED, 1, 3")
    ap("*STEP, NLGEOM=NO"); ap(" 1) dysbiotic biofilm growth"); ap("*STATIC")
    ap("*TEMPERATURE"); ap(" ALLN, 1.0")
    ap("*OUTPUT, FIELD"); ap("*NODE OUTPUT"); ap(" U")
    ap("*ELEMENT OUTPUT, POSITION=CENTROID"); ap(" S, COORD"); ap("*END STEP")
    # occlusal load: ISO 14801-style 30deg oblique for the leveled-up generic job (lateral/axial =
    # tan30 = 0.577); the preserved tier2b_real keeps its original ~11deg (0.2) load.
    lat = 0.577 if job != "tier2b_real" else 0.2
    ang = "30deg oblique (ISO 14801)" if job != "tier2b_real" else "near-axial"
    ap("*STEP, NLGEOM=NO"); ap(" 2) occlusal load on both crowns (%s)" % ang); ap("*STATIC")
    ap("*CLOAD")
    for ids in (too_top, load_top):
        if len(ids):
            f = 100.0 / len(ids)            # ~100 N occlusal resultant per crown
            for n in np.unique(ids) + 1:
                ap(" %d, 3, %.5f" % (n, -f)); ap(" %d, 1, %.5f" % (n, lat * f))
    ap("*OUTPUT, FIELD"); ap("*NODE OUTPUT"); ap(" U")
    ap("*ELEMENT OUTPUT, POSITION=CENTROID"); ap(" S, COORD"); ap("*END STEP")
    open(f"{OUT}/{job}.inp", "w").write("\n".join(L) + "\n")

    # meta for visualisation
    matarr = np.empty(gid - 1, dtype=object)
    for name, ids, conn in elem_rows:
        matarr[ids - 1] = name
    allconn = np.vstack([c for _, _, c in elem_rows])
    np.savez(f"{OUT}/{job}_meta.npz", nodes=nodes, conn=allconn,
             mat=matarr.astype(str), cent=nodes[allconn].mean(axis=1))
    print("wrote %s.inp  total elems=%d nodes=%d" % (job, gid - 1, len(nodes)))


if __name__ == "__main__":
    main()
