"""The helper module uses dask_kubernetes to deploy a scalable KubeCluster and pods
The pods are only used for passing names now but can be used for much more complex situations
You can also create your own k8s pods to be used for this purpose

Idea for getting the local pod resources are taken from `genkubejob` script
"""

import os
import json
import time
from subprocess import run, PIPE, Popen
from dask_kubernetes.core import KubeCluster


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


def generate_kubecluster_spec(
    name="vasp-pod", cpu=4, memory="2Gi", extra_envs={}, local_pod={}
):
    """Generate a minimal spec for a dask cluster scheduler from local spec
    The following values in the "spec" field is inherited:
    - containers
    - volumes
    - volumeMounts

    To make a cluster, simple do something like:
    ```
    cluster = KubeCluster(generate_kubecluster_spec(**kwargs))
    cluster.scale(2)
    ```
    """
    local_podname = local_pod.get("podname", None)
    namespace = local_pod.get("namespace", None)
    local_spec_dict = get_local_kube_config(local_podname, namespace)

    # Now construct a minimal spec for KubeCluster
    minimal_spec = {"metadata": {}, "spec": {}}

    # Container settings
    containers = local_spec_dict["spec"]["containers"]
    if len(containers) > 1:
        raise NotImplementedError("Can only handle 1 container per pod now")
    container = containers[0]
    resources = {
        "limits": {"cpu": cpu, "memory": memory},
        "requests": {"cpu": cpu, "memory": memory},
    }
    container["resources"] = resources
    container["name"] = name
    envs = container["env"].copy()
    for k, v in extra_envs.items():
        to_append = True
        # Do we need to update the old env?
        for env in envs:
            if env["name"] == k:
                env["value"] = v
                to_append = False
                break
        if to_append:
            envs.append({"name": k, "value": v})
    container["env"] = envs

    # Add dask-specific args to enable scheduling
    dask_args = [
        "dask-worker",
        "--nthreads",
        "1",
        "--no-dashboard",
        "--death-timeout",
        "60",
    ]
    container["args"] = dask_args
    container["command"] = None

    # Add containers back to spec
    minimal_spec["spec"]["containers"] = containers
    minimal_spec["spec"]["volumes"] = local_spec_dict["spec"]["volumes"]

    # Additional flags
    minimal_spec["spec"]["restartPolicy"] = "Never"
    return minimal_spec


def get_kubecluster_pods(cluster):
    """Get series of worker pods name / namespace dict from current cluster"""
    worker_pods = cluster.workers
    pods_dict = {}
    # following code from dask source
    for number, pod in worker_pods.items():
        podname = pod._pod.metadata.name
        namespace = pod.namespace
        pods_dict[number] = {"name": podname, "namespace": namespace}
    return pods_dict


def is_pod_ready(name, namespace):
    pod_dict = get_local_kube_config(name, namespace)
    #     print(pod_dict["status"])
    conditions = pod_dict["status"]["conditions"]
    if len(conditions) == 0:
        all_ready = False
    else:
        try:
            all_ready = all([cond["status"] == "True" for cond in conditions])
        except KeyError:
            all_ready = False
    return all_ready


def create_kube_pods(
    scale=1,
    name="vasp-pod",
    cpu=4,
    memory="2Gi",
    extra_envs={},
    local_pod={},
    wait_for_pods=True,
    timeout=60,
):
    pod_spec = generate_kubecluster_spec(
        name=name, cpu=cpu, memory=memory, extra_envs=extra_envs, local_pod=local_pod
    )
    cluster = KubeCluster(pod_spec, n_workers=scale)
    # wait for pods to be ready in a loop
    #     print(cluster, cluster.workers)
    pods_dict = get_kubecluster_pods(cluster)
    #     print(pods_dict)
    if wait_for_pods:
        wait_time = 0
        all_ready = False
        while wait_time < timeout:
            pod_statuses = []
            for number, pod in pods_dict.items():
                pod_statuses.append(is_pod_ready(pod["name"], pod["namespace"]))
            #             print(pod_statuses)
            if all(pod_statuses) and (len(pod_statuses) > 0):
                all_ready = True
                break
            else:
                t_sleep = 0.5
                time.sleep(t_sleep)
                wait_time += t_sleep
        if not all_ready:
            raise RuntimeError(
                f"Not all worker pods on scheduler {cluster} are ready within the timeout!"
            )
    return cluster, pods_dict
