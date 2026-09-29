"""protocol.py — wire schema for the Gauss-point material bridge.

A single material evaluation is one request → one response, newline-delimited
JSON (easy to debug; swap for a binary/MessagePack frame later without touching
the physics). All tensors are row-major length-9 lists; Voigt vectors follow the
**Abaqus order 11,22,33,12,13,23** (the Python reference core's convention).

Request  : {F:[9], Fv:[9], alpha, C10, C01, D1, eta, mtype, dt}
Response : {stress:[6], Fv_new:[9], detFe, dsdePl:[36]}  |  {error:"..."}

A second request "kind" shares the same server/socket: the 0D Hamilton
ecology ODE step (ecology_jax.py), used to evolve a Gauss point's local
composition state live instead of reading a precomputed alpha field
(coupling/README.md next-steps #4). A plain material request carries no
"kind" key, so this is backward compatible with kUsePy=1 clients that
predate it.

Ecology request  : {kind:"ecology", g:[12], theta:[20], dt_h, n_sub?}
Ecology response : {g_new:[12], phi_int}  |  {error:"..."}

`n_sub` sub-divides the increment. The ecology ODE has its own time scale --
measured to turn erratic above dt ~ 1.5e-3 -- which is far finer than the
pseudo-time steps a mechanical deck is built around, and a solver's own
bisection cannot bridge a factor of a thousand. Sub-stepping inside one call
does bridge it, at one socket round trip rather than n_sub of them, and it is
not a new scheme: every reference implementation in `apdl/` already advances
the ODE by chaining `ecology_step`, which is exactly what this does.

`phi_int` is the accumulated integral of the living fraction over those
sub-steps, sum_k dt_sub * phi_tot(g_k). The caller needs it because growth
accumulates as k_alpha * that integral, and with n_sub > 1 the final state's
phi_tot alone no longer stands for the whole increment. At n_sub = 1 it is
dt_h * phi_tot(g_new), which is what a caller computed for itself before this
field existed, so the change is backward compatible in value as well as shape.
"""
from __future__ import annotations

import json

REQ_KEYS = ("F", "Fv", "alpha", "C10", "C01", "D1", "eta", "mtype", "dt")
ECO_REQ_KEYS = ("kind", "g", "theta", "dt_h")   # n_sub is optional


def encode_request(F, Fv, alpha, C10, C01, D1, eta, mtype, dt) -> bytes:
    return (json.dumps({
        "F": list(map(float, F)), "Fv": list(map(float, Fv)),
        "alpha": float(alpha), "C10": float(C10), "C01": float(C01),
        "D1": float(D1), "eta": float(eta), "mtype": float(mtype),
        "dt": float(dt),
    }) + "\n").encode()


def decode_request(line: bytes) -> dict:
    d = json.loads(line)
    missing = [k for k in REQ_KEYS if k not in d]
    if missing:
        raise ValueError(f"request missing keys: {missing}")
    return d


def encode_response(stress, Fv_new, detFe, dsdePl) -> bytes:
    return (json.dumps({
        "stress": list(map(float, stress)),
        "Fv_new": list(map(float, Fv_new)),
        "detFe": float(detFe),
        "dsdePl": list(map(float, dsdePl)),
    }) + "\n").encode()


def encode_error(msg: str) -> bytes:
    return (json.dumps({"error": str(msg)}) + "\n").encode()


def decode_response(line: bytes) -> dict:
    return json.loads(line)


def encode_ecology_request(g, theta, dt_h, n_sub: int = 1) -> bytes:
    """One ecology request. `dt_h` is the whole increment; `n_sub` divides it.

    Omitting `n_sub` sends the pre-existing single-step request, so a client
    built against the older schema keeps working unchanged.
    """
    if n_sub < 1:
        raise ValueError(f"n_sub must be at least 1, got {n_sub}")
    d = {
        "kind": "ecology",
        "g": list(map(float, g)),
        "theta": list(map(float, theta)),
        "dt_h": float(dt_h),
    }
    if n_sub != 1:
        d["n_sub"] = int(n_sub)
    return (json.dumps(d) + "\n").encode()


def decode_ecology_request(line: bytes) -> dict:
    d = json.loads(line)
    missing = [k for k in ECO_REQ_KEYS if k not in d]
    if missing:
        raise ValueError(f"ecology request missing keys: {missing}")
    return d


def encode_ecology_response(g_new, phi_int=None) -> bytes:
    d = {"g_new": list(map(float, g_new))}
    if phi_int is not None:
        d["phi_int"] = float(phi_int)
    return (json.dumps(d) + "\n").encode()
