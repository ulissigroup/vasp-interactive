"""Plot the scf per parent call in online mlp 
"""
import matplotlib.pyplot as plt
import numpy as np
fig, ax = plt.subplots(1, 1, figsize=(8, 4))
d1 = np.load("elec_steps_vasp.npy")
d2 = np.load("elec_steps_VaspInteractive.npy")

# Use the cutoff to trim data. The mlp demo seems to give different parent calls
ax.plot(d1[:12], label="VASP single point")
ax.plot(d2[:12], label="VASP Interactive")
ax.set_xlabel("No. Parent Calls")
ax.set_ylabel("Electronic SCF Steps")
ax.legend()
ax.set_title("Machine Learning Potential Optimizer for Cu$_7$ clusters")
fig.tight_layout()
fig.savefig("mlp_online_parent_scf.png")