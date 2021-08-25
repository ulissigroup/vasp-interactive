"""Plot the scf per parent call in online mlp 
"""
import matplotlib.pyplot as plt
import numpy as np

file_names = ["vasp_bfgs_nocache", "vasp_bfgs_cache", "vpi_bfgs", "vasp_al", "vpi_al"]
label_names = [
    "VASP+BFGS nocache",
    "VASP+BFGS cache",
    "VaspInteractive+BFGS",
    "VASP+MLP",
    "VaspInteractive+MLP",
]


def get_scf(i):
    import re

    txt_file = "bench_" + file_names[i] + ".txt"
    elec_steps = []
    with open(txt_file, "r") as fd:
        lines = fd.readlines()
    pat = r"DAV\:\s+([\d]+)"
    prev = -1
    elec_steps = []
    for l in lines:
        m = re.match(pat, l)
        print(m)
        if m:
            nex = int(m[1])
            if nex > prev:
                prev = nex
            else:
                elec_steps.append(prev)
                prev = -1
    elec_steps.append(nex)
    return elec_steps


def main():
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    times = np.load("times.npy")
    ax1 = ax[0]
    ax1.bar(np.arange(3), times[:3], label="DFT + BFGS")
    ax1.bar(np.arange(3, 5), times[3:], label="ΔCalc@MLP + BFGS")
    ax1.set_xticks(np.arange(5))
    ax1.set_xticklabels(label_names, rotation=45, ha="right")
    ax1.set_ylabel("Wall Time (s)")
    ax1.legend()

    ax2 = ax[1]

    for i in range(5):
        elec_steps = get_scf(i)
        print(i, elec_steps)
        name = label_names[i]
        tot_elec = np.hstack([0, np.cumsum(elec_steps)])
        ax2.plot(np.arange(len(elec_steps) + 1), tot_elec, label=name)
    ax2.set_xlabel("DFT Calls")
    ax2.set_ylabel("Cumulative electronic scf steps")
    ax2.legend()

    fig.tight_layout()
    fig.savefig("mlp_online_parent_scf.png")


if __name__ == "__main__":
    main()
