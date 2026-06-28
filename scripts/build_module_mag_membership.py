#!/usr/bin/env python3
"""Consolidate the Louvain-module membership of every co-occurrence network into a
single MAG-resolved JSON.

The per-network module-membership files list members as ASV iterativeIDs; the
co-occurrence networks are already restricted to MAG-representative ASV nodes, so
every member maps 1:1 to a MAG (CAN_x_bin.N) via mag_iterativeID_old_to_new.json.
This writes each module's membership as the underlying MAGs (bin name + iterativeID).

Input
  network/network_module_membership_p_value_FDR.json            (aggregate network)
  network/network_module_membership_p_value_FDR_phase{I,III,IV,V}.json
  mag_iterativeID_old_to_new.json                               (MAG bin -> iterativeID)
Output
  network/louvain_module_mag_membership.json

Run:  ~/Documents/py_venv/bin/python scripts/build_module_mag_membership.py
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NET = ROOT / "network"

# iterativeID -> MAG bin name (the network nodes are MAG-representative, 1:1)
iid2mag = {v: k for k, v in json.loads((ROOT / "mag_iterativeID_old_to_new.json").read_text()).items()}

# (output network key, source membership file)
NETWORKS = [
    ("aggregate", "network_module_membership_p_value_FDR.json"),
    ("phaseI",    "network_module_membership_p_value_FDR_phaseI.json"),
    ("phaseIII",  "network_module_membership_p_value_FDR_phaseIII.json"),
    ("phaseIV",   "network_module_membership_p_value_FDR_phaseIV.json"),
    ("phaseV",    "network_module_membership_p_value_FDR_phaseV.json"),
]

out = {
    "_description": (
        "MAG-resolved Louvain-module membership of each co-occurrence network. "
        "Members are the MAG-representative ASV nodes (iterativeID) with their "
        "underlying MAG bin name. Networks: aggregate + per-phase (phase II has no "
        "FDR-passing co-occurrence network, so no modules)."
    ),
    "networks": {},
}

for key, fn in NETWORKS:
    src = json.loads((NET / fn).read_text())
    mod_keys = sorted((k for k in src if k.startswith("module")),
                      key=lambda k: int(k.split("_")[1]))
    net = {}
    n_drop = 0
    for mk in mod_keys:
        members = []
        for iid in src[mk]["members"]:
            mag = iid2mag.get(iid)            # MAG-filter: keep only members that map to a MAG
            if mag is None:
                n_drop += 1
                continue
            members.append({"iterativeID": iid, "MAG": mag})
        net[mk] = {"size": len(members), "type": src[mk].get("type", "louvain"),
                   "members": members}
    out["networks"][key] = net
    total = sum(m["size"] for m in net.values())
    print(f"{key:9}: {len(net)} modules, {total} MAG members"
          + (f"  ({n_drop} non-MAG members dropped)" if n_drop else ""))

dest = NET / "louvain_module_mag_membership.json"
dest.write_text(json.dumps(out, indent=2) + "\n")
print(f"\nwrote {dest.relative_to(ROOT)}")
