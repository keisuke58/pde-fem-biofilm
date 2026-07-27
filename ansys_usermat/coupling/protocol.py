"""protocol.py — wire schema for the Gauss-point material bridge.

A single material evaluation is one request → one response, newline-delimited
JSON (easy to debug; swap for a binary/MessagePack frame later without touching
the physics). All tensors are row-major length-9 lists; Voigt vectors follow the
**Abaqus order 11,22,33,12,13,23** (the Python reference core's convention).

Request  : {F:[9], Fv:[9], alpha, C10, C01, D1, eta, mtype, dt}
Response : {stress:[6], Fv_new:[9], detFe, dsdePl:[36]}  |  {error:"..."}
"""
from __future__ import annotations

import json

REQ_KEYS = ("F", "Fv", "alpha", "C10", "C01", "D1", "eta", "mtype", "dt")


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
