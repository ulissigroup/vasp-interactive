#!/usr/bin/env python
import re
import os
from argparse import ArgumentParser
from pathlib import Path
from copy import copy
import shutil

fpath = Path(__file__)
fdir = fpath.parent

patch_main_F = """!-----------------------------------------------------------------------
! interactive mode
!-----------------------------------------------------------------------
        ELSE IF (DYN%IBRION==11) THEN
io_begin
           IF (IO%LOPEN) CALL WFORCE(IO%IU6)            ! write of OUTCAR
           IF (IO%LOPEN) CALL WFORCE(17)                ! write of OSZICAR
           IF (IO%LOPEN) CALL XML_FLUSH                 ! force flush xml file
io_end
           CALL INPOS(LATT_CUR, T_INFO, DYN, IO%IU6, IO%IU0, INFO%LSTOP, WDES%COMM)
           CALL INLATT(LATT_CUR, T_INFO, DYN, IO%IU6, IO%IU0, INFO%LSTOP, WDES%COMM)
        ENDIF

"""

patch_poscar_F = """

!=======================================================================
!
! read direct lattice from stdin
! Patch by T.Tian (alchem0x2a) for upgrading VASP's interactive mode
! The implementation should not interfere with normal VASP usage but 
! do tests at your own risk!
!
!=======================================================================
      SUBROUTINE INLATT(LATT_CUR, T_INFO, DYN, IU6, IU0, LSTOP, MYCOMM)
      USE prec
      USE lattice
      USE main_mpi
      IMPLICIT NONE

      TYPE (latt)::       LATT_CUR
      TYPE (type_info) :: T_INFO
      TYPE (dynamics)  :: DYN
      INTEGER :: IU6, IU0
      LOGICAL :: LSTOP
      TYPE (communic) :: MYCOMM
    ! local
      INTEGER NI, IERR, I, J
      REAL(q) :: LATT_A_OLD(3, 3), LATT_B_OLD(3, 3)

      ! Local copy of the old lattice parameters
      LATT_A_OLD = LATT_CUR%A
      LATT_B_OLD = LATT_CUR%B
      ! Update the lattice on I/O rank else 0
      ! Use M_sum_b to gather results
      IF (IU6>=0) THEN
         IF (IU0>=0) WRITE(IU0,'(A)') 'LATTICE: reading from stdin'
         DO NI=1,3
            READ(*,*,  IOSTAT=IERR) LATT_CUR%A(1,NI), LATT_CUR%A(2,NI), LATT_CUR%A(3,NI)
            IF (IERR/=0) EXIT
         ENDDO

         IF (IERR/=0) THEN
            LSTOP=.TRUE.
         ELSE
            LSTOP=.FALSE.
         ENDIF
         CALLMPI( M_sum_i(MYCOMM, IERR, 1))
         IF (IERR==0) THEN
            CALLMPI( M_sum_d(MYCOMM, LATT_CUR%A(1,1), 9))
         ENDIF
      ELSE
         IERR=0
         CALLMPI( M_sum_i(MYCOMM, IERR, 1))
         IF (IERR==0) THEN
            LATT_CUR%A(:,1:3)=0
            CALLMPI( M_sum_d(MYCOMM, LATT_CUR%A(1,1), 9))
         ENDIF
      ENDIF
      ! Update other lattice params (reciprocal, volume, norm vect) on all ranks
      CALL LATTIC(LATT_CUR)
      ! If current rank allows output, write the positions
      IF (IU0 >= 0) THEN
            WRITE(IU0, '(A,3X,A,19X,A)') 'New', 'direct lattice vectors', 'reciprocal lattice vectors'
            WRITE(IU0, 8900) ((LATT_CUR%A(I,J),I=1,3),(LATT_CUR%B(I,J),I=1,3),J=1,3)
            WRITE(IU0, '(A,3X,A,19X,A)') 'Old', 'direct lattice vectors', 'reciprocal lattice vectors'
            WRITE(IU0, 8900) ((LATT_A_OLD(I,J),I=1,3),(LATT_B_OLD(I,J),I=1,3),J=1,3)
            WRITE(IU0, '(A)') 'LATTICE: read from stdin'
      ENDIF
8900  FORMAT(3(2(3X,3F13.7)/)/)
      END SUBROUTINE INLATT



"""


def patch_txt(old_file, pattern, patch_content, replace_func):
    """Backup the old file and write the patch"""
    old_file = Path(old_file)
    txt = open(old_file, "r").read()
    repl = replace_func(patch_content)
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
            r"^\!\-*\n.*?interactive\smode.*?\n"
            r"[\s\S]*?ELSE\sIF\s\(DYN\%IBRION==11\)"
            r"[\s\S]*?ENDIF[\s\S]*?\n{2}"
        ),
        "patch_content": patch_main_F,
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
        "patch_content": patch_poscar_F,
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
        patch_content = patch["patch_content"]
        pattern = patch["pattern"]
        func = patch["replace_func"]
        print(f"Making patch for {old_file.as_posix()}")
        patch_txt(old_file, pattern, patch_content, func)
    print("Success")
    return


if __name__ == "__main__":
    main()
