"""The helper module uses dask_kubernetes to deploy a scalable KubeCluster and pods
The pods are only used for passing names now but can be used for much more complex situations
You can also create your own k8s pods to be used for this purpose

Idea for getting the local pod resources are taken from `genkubejob` script
"""

import os
import json
from subprocess import run, PIPE, Popen
from dask_kubernetes import KubeCluster


def get_local_kube_config(podname=None, namespace=None):
    """Use kubectl to detect the kube config spec as a dict
    if podname and namespace are not given, use the local $HOSTNAME as the podname
    """
    args_ns = [f"--namespace={namespace}"] if namespace is not None else []
    hostname = os.environ.get("HOSTNAME", None) if podname is None else podname
    if hostname is None:
        raise ValueError(
            "You should either specify the podname, or set the HOSTNAME environmental variable"
        )

    cmds = ["kubectl", "get", "pods", hostname] + args_ns + ["-o", "json"]
    proc = run(cmds, stdout=PIPE)
    spec_json = json.loads(proc.stdout.decode("utf8"))
    return spec_json
