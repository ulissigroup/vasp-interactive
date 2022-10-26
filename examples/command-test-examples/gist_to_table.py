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
    "cori_hsw": "Cori Haswell (†)",
    "cori_knl": "Cori KNL (†)",
    "perlmutter_cpu": "Perlmutter CPU (†)",
    "perlmutter_gpu": "Perlmutter GPU (†)",
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
        if ("," not in line) or (line.startswith("#")):
            continue
        version, state = list(map(lambda s: s.strip(), line.split(",")))
        state = state.strip()
        state = state.replace(" ", "_")
        tpc = "-tpc" in version
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
    return res_dict


def render(res_dict):
    header = (
        "|"
        + " " * 10
        + "".join([f"| {disp_name} " for fn, disp_name in tests.items()])
        + " |"
    )
    hline = "".join(["| ---------- " for i in range(len(tests) + 1)]) + " |"
    body = []
    for ver in vasp_versions:
        for tpc in (False, True):
            row = res_dict[ver][tpc]
            if all([row[sys] == "not_available" for sys in tests]):
                # print(f"{ver}-{tpc} does not exist. pass")
                continue
            vasp_name = f"VASP {ver}.x"
            if tpc:
                vasp_name += " - TPC (*)"
            line = f"| {vasp_name} "
            for sys in tests:
                state = row[sys]
                if state is None:
                    state = "not_available"
                line += "| " + gen_badge(state) + " "
            line += " |"
            body.append(line)
    text = "\n".join([header, hline] + body) + "\n"
    return text


def main():
    import argparse
    from pathlib import Path
    curdir = Path(__file__).parent
    readme = curdir.parents[1] / "README.md"
    readme_bk = curdir.parents[1] / "README.md.bk"
    print(readme.is_file())
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--update", action="store_true")
    args = parser.parse_args()
    res_dict = parse_all_states()
    update_text = render(res_dict)
    if not args.update:
        print("Output replacement only")
        print(update_text)
        return
    readme_content = open(readme, "r").readlines()
    new_content = []
    begin, end = 0, 0
    for i, line in enumerate(readme_content):
        if "PLACEHOLDER BEGIN" in line:
            begin = i
        elif "PLACEHOLDER END" in line:
            end = i
        else:
            pass
    head = readme_content[:begin + 1]
    tail = readme_content[end:]
    new_content = head + ["\n" + update_text + "\n"] + tail
    # print(new_content)

    with open(readme_bk, "w") as fd:
        fd.writelines(readme_content)
    print("Backed up readme")

    with open(readme, "w") as fd:
        fd.writelines(new_content)
    print("Updated readme")

if __name__ == "__main__":
    main()
        
        
