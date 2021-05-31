# Interactive Vasp calculator
This repo provides a patched version of the `VaspInteractive` calculator, 
originally from Atomic Simulation Environment (`ase`).

Improvement from original `VaspInteractive` include:
- Accepting same constructor as normal `Vasp` calculator.
- Better process management to avoid orphan process when exiting.
- A context mode to run the `VaspInteractive` for improved garbage collection.
- Parsing OUTCAR for ionic step results instead of vasprun.xml which can be corrupted

## How does it work?
`VaspInteractive` invokes the interactive mode in VASP calculations by setting the keyword `INTERACTIVE = .TRUE.` in INCAR file.
This option is equivalent to `IBRION = 11` (although no force information printed to `stdout`).

In the interactive mode, after reading the initial input files (INCAR, POTCAR, POSCAR, KPOINTS), 
VASP pauses at the end of the first ionic step
and asks for the scaled coordinates of the next ionic step. Once user types the scaled coordinates to `stdin`, 
VASP performs the next ionic step using previous density and wavefunction. 

The ionic cycles of interactive mode VASP can be terminated by any of the following:
i) setting `NSW` values
ii) writing STOPCAR file to the calculation directory
iii) invalid inputs to stdin (such as `Ctrl+C`)

`VaspInteractive` uses method ii) to stop the ionic cycles. In general, `VaspInteractive` can save up to 50% of the wall 
time compared with classic `Vasp` calculator (combined with ASE optimizers such as BFGS), since less electronic steps are
required and program spin-up time is drastically reduced.

## How to use
Only dependency is `ase`. Hopefully in the near future the `VaspInteractive` calculator will be merged with ASE upstream.

- Install via `pip`

    ```sh
        pip install git+https://github.com/ulissigroup/vasp-interactive-test.git
    ```
    
- Minimal example
    
    The following code shows how to run an BFGS-optimization of H2 molecule using `VaspInteractive`:
    
    ```python
    from ase.optimize import BFGS
    from ase.build import molecule
    from vasp_interactive import VaspInteractive
    atoms = molecule("H2", vacuum=4, pbc=True)
    with VaspInteractive(xc="pbe") as calc:
        atoms.calc = calc
        dyn = BFGS(atoms)
        dyn.run(fmax=0.05)
    ```
    
    The only difference using `VaspInteractive` compared with pure VASP+BFGS, 
    is that `VaspInteractive` is wrapped inside a context manager. The `with` statemen is 
    the recommended way to run `VaspInteractive` as it can seamlessly stop the VASP process
    associated with the calculator outside the context. Although you can also use the classic way of ASE-calculators, 
    see next section
    
- Context mode vs classic mode
    