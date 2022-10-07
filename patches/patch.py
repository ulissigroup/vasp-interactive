#!/usr/bin/env python
import re
import os
from argparse import ArgumentParser
from pathlib import Path
from copy import copy
import shutil

fpath = Path(__file__)
fdir = fpath.parent


def patch_txt(old_file, pattern, patch_file, replace_func):
    """Backup the old file and write the patch"""
    old_file, patch_file = Path(old_file), Path(patch_file)
    txt = open(old_file, "r").read()
    repl = replace_func(open(patch_file, "r").read())
    matches = re.findall(pattern, txt, flags=(re.MULTILINE | re.IGNORECASE))
    assert len(matches) == 1
    # print(matches[0])
    txt_new = re.sub(pattern, repl, txt, count=1, flags=(re.MULTILINE | re.IGNORECASE))
    backup_file = old_file.with_suffix(old_file.suffix + ".bak")
    shutil.copy(old_file, backup_file)
    with open(old_file, "w") as fd:
        fd.write(txt_new)
    return


patches = []
patches.append(
    {
        "name": "main.F",
        "desc": """Replace pattern
                    !-----------------------------------------------------------------------
                    ! interactive mode
                    !-----------------------------------------------------------------------
                        ELSE IF (DYN%IBRION==11) THEN
                        ...
                        ENDIF
                    """,
        "pattern": (
            r"\!\-*\n.*?interactive\smode.*?\n"
            r"[\s\S]*?ELSE\sIF\s\(DYN\%IBRION==11\)"
            r"[\s\S]*?ENDIF[\s\S]*?\n"
        ),
        "file": "main.F.patch",
        "replace_func": lambda s: s,
    }
)
patches.append(
    {
        "name": "poscar.F",
        "desc": """Replace pattern
                    END SUBROUTINE INPOS
                    
                    ... 
                    
                    !*************************SUBROUTINE OUTPOS_TRAIL  *********************
                    """,
        "pattern": (
            r"(^\s*?END\sSUBROUTINE\sINPOS.*?$\n{2})"
            r"([\s\S\n]*?)"
            r"(^\!\*+?SUBROUTINE\sOUTPOS_TRAIL)"
        ),
        "file": "poscar.F.patch",
        "replace_func": lambda s: r"\1" + s + r"\3",
    }
)


def main():
    parser = ArgumentParser(
        description="Patch VASP Fortran files for better interactive mode integration"
    )
    parser.add_argument("src", type=str, help="Path to VASP source code")
    args = parser.parse_args()
    src = Path(os.path.expanduser(args.src)).resolve()
    for patch in patches:
        old_file = src / patch["name"]
        patch_file = fdir / patch["file"]
        pattern = patch["pattern"]
        func = patch["replace_func"]
        print(f"Making patch for {old_file.as_posix()}")
        patch_txt(old_file, pattern, patch_file, func)
    print("Success")
    return


if __name__ == "__main__":
    main()
