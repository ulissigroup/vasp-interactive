"""Simple example showing isolating vasp process in separate kubernetes pod
The example requires to be run within an actual kubernetes pod with `kubectl` access.
Dask-kubernetes serves as the backend for pod scheduling

The example uses shared volume mounting on both local and remote pods, 
so it may not work when setting directory to tmpdir.

Assume the volume is mounted to /home/jovyan/dat
"""
import time
from dask_kubernetes import KubeCluster
from ase.build import molecule
from vasp_interactive.kubernetes import KubeVaspInteractive, generate_kubecluster_spec, get_kubecluster_pods
from pathlib import Path


mol1 = molecule("CH4", vacuum=5, pbc=True)
mol2 = mol1.copy()
mol2.rattle()
params = dict(xc="pbe", kpts=1, encut=350)

# Generate pods
pod_spec = generate_kubecluster_spec(cpu=8, memory="2Gi")
cluster = KubeCluster(pod_spec)
print("Created a KubeCluster for scheduling...")
image = pod_spec["spec"]["containers"][0]["image"]
print(f"Pods running on image {image}")

cluster.scale(2)
print("Scale to 2 worker pods...")
# Refresh to get worker status
while True:
    if cluster.workers != {}:
        break
    else:
        time.sleep(0.5)
        

# Get individual pod name & namespace
worker_pods = get_kubecluster_pods(cluster)
for number, pod in worker_pods.items():
    print(f"Pod {number}, name: {pod['name']}, namespace: {pod['namespace']}")



cluster.close()
                      



