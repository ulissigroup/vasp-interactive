#!/usr/bin/env python3
"""Convert the raw output of compatibility test stored in gist
https://gist.github.com/alchem0x2A/afede700c2b7703c77e10e51333bfa75
to a readme table with coloring from shield.io
"""
import requests
import os

default_gist_id = "afede700c2b7703c77e10e51333bfa75"
if os.environ.get("GIST_ID", None) is not None:
    default_gist_id = os.environ.get("GIST_ID")

tests = {
    "ulissi_docker": "Docker images (**)",
    "cori_hsw": "Cori Haswell (†)"
    "cori_knl": "Cori KNL (†)"
    "perlmutter_cpu": "Perlmutter CPU (†)"
    "perlmutter_gpu": "Perlmutter GPU (†)"
}

all_states = {
    "incompatible": "red",
    "minimal_support": "orange",
    "partial_pass": "olive",
    "all_pass": "green",
    "not_available": "lightgrey",
}

vasp_versions = ("5.4", "6.1", "6.2", "6.3")

def gen_badge(msg):
    assert msg in all_states.keys()
    color = all_states[msg]
    link = f"![](https://img.shields.io/badge/-{msg}-{color})"
    return link

def parse_txt(txt):
    lines = txt.split("\n")
    all_matches = []
    for line in lines:
        if "," not in line:
            continue
        version, state = list(map(lambda s: s.strip(), line.split(",")))
        state = state.strip()
        state = state.replace(" ", "_")
        tpc = ("-tpc" in version)
        ver_maj = None
        for ver_dig in vasp_versions:
            if ver_dig in version:
                ver_maj = ver_dig
                break
        assert ver_maj is not None, "Cannot read version!"
        all_matches.append((ver_maj, tpc, state))
    return all_matches

def read_gist(gist_id=None, file_name=None):
    # Read from gist
    if gist_id is None:
        gist_id = default_gist_id
    url = f"https://gist.githubusercontent.com/alchem0x2A/{gist_id}/raw/{file_name}.txt"
    req = requests.get(url)
    if req.status_code != requests.codes.ok:
        raise RuntimeError(f"Getting url {url} failed")
    text = req.text
    return text

def parse_all_states():
    res_dict = {}
    # Init dict
    for ver in vasp_versions:
        d = {}
        for tpc in [False, True]:
            dd = {}
            for sys in tests.keys():
                dd[sys] = "not_available"
            d[tpc] = dd
        res_dict[ver] = d
    
    for fn in tests.keys():
        text = read_gist(file_name=fn)
        all_matches = parse_txt(text)
        for res in all_matches:
            ver, tpc, state = res
            res_dict[ver][tpc][fn] = state
        
    header = "|" + " " * 10 + [f"| {disp_name} " for fn, disp_name in tests.items()] + " |"
    body = []
    for ver in vasp_versions:
        for tpc in (False, True):
            row = res_dict[ver][tpc]
            if all([row[sys] is None for sys in tests]):
                print(f"{ver}-{tpc} does not exist. pass")
                continue
            vasp_name = f"VASP {ver}.x"
            if tpc:
                vasp_name += " - TPC (*)"
            line = f"| {tpc} "
            for sys in tests:
                state = row[sys]
                if state is None:
                    state = "not_available"
                line += "| " + gen_badge(badge) + " "
            line += " |"
    text = "\n".join([header] + body)
    return text
     

