"""Running vasp interactive in separate remote resources, other than the local processes
   Currently the following remote(s) are implemented
   - kubernetes
   
   Be careful when using this module for production purpose.
"""
import sys
from pathlib import Path
from vasp_interactive import VaspInteractive
from subprocess import run, Popen, PIPE


def _gen_kubectl_cmd(cmd, pod_name, pod_namespace=None):
    """easy wrapper to generate cmds for kubectl"""
    if isinstance(cmd, (list, tuple)):
        cmd_string = " ".join(cmd)
    elif isinstance(cmd, str):
        cmd_string = cmd
    else:
        raise TypeError("Can only take command as string / list tuple")

    args_namespace = (
        [f"--namespace={pod_namespace}"] if pod_namespace is not None else []
    )
    commands = (
        ["kubectl", "exec", "-i", pod_name]
        + args_namespace
        + ["--", "bash", "-c", cmd_string]
    )
    return commands


class KubeVaspInteractive(VaspInteractive):
    """Running VaspInteractive in separate kubernetes pod using server-client mode
    instead of using local process of N-cores,
    this calculator use `kubectl exec` to open pipe with a pod and read stdin from there.
    `KubeVaspInteractive` can be useful for situations where individual VASP relaxations need to be gathered,
    e.g. NEB or Monte Carlo Samping.

    If synchronization between individual relaxation processes are not important, original VaspInteractive running solely inside
    kubernete pods (e.g. wrapped by dask worker) may be more straightforward.

    Most parameters should be similar with local VaspInteractive, but note the following:
    `directory`: the directory on local machine, by default where the input files can be generated
    `command`: the command to be run inside the pod. Additional kubectl directives will be added automatically

    `remote_directory`: default the same as on local if not set. The directory inside the pod
    `pod`: dict with keys `name` and `namespace` to work with kubectl exec

    Note the calculator does not take care of how the kubernetes pods are created / managed.
    Currently the calculator only supports that the pod and local file systems are mounted on same volume with same hierachy.
    supporting on syncing / remote transfer will be added later.
    """

    def __init__(
        self,
        atoms=None,
        directory=".",
        remote_directory=None,
        label="vasp-interactive",
        command="$VASP_COMMAND",
        txt="vasp.out",
        allow_restart_process=True,
        pod={},
        **kwargs,
    ):
        super(KubeVaspInteractive, self).__init__(
            atoms=atoms,
            directory=directory,
            command=command,
            txt=txt,
            allow_restart_process=allow_restart_process,
            **kwargs,
        )
        pod_name = pod.get("name", None)
        if pod_name is None:
            raise ValueError("Name for the pod must be provided.")
        self.pod_name = pod_name
        self.pod_namespace = pod.get("namespace", None)

        if remote_directory is None:
            # makesure absolute path
            self.remote_directory = Path(self.directory).resolve()
        else:
            raise NotImplementedError(
                "KubeVaspInteractive not (yet) supporting remote_directory that differs from local."
            )

        # For running kubectl, the actual command is wrapped around kubectl
        # for lazy testing use the $VASP_COMMAND environ in the pod

        # TODO: how about another make_cmds?
        # warning, no verbatim brack needed for this part
        pod_command = f"cd {self.remote_directory.as_posix()} && {command}"
        self._args = _gen_kubectl_cmd(pod_command, self.pod_name, self.pod_namespace)

        return

    def _kubectl_exec(self, command):
        """Use subprocess to run a command via kubectl, blocking
        currently, output not redirected
        """
        cmds = _gen_kubectl_cmd(command, self.pod_name, self.pod_namespace)
        proc = run(cmds)
        return proc.returncode

    def _force_kill_process(self, proc_pattern=".*vasp.*"):
        """kill the process via kubectl can be tricky. As the last resort, try to use killall"""
        try:
            self.close()
        except Exception as e:
            print(
                (
                    f"Trying to close the VASP stream but encountered error: \n"
                    f"{e}\n"
                    "Will now force closing the VASP process. "
                    "The OUTCAR and vasprun.xml outputs may be incomplete"
                ),
                file=sys.stderr,
            )
            # match all processes in current pod
            returncode = self._kubectl_exec(f"killall -v -r '{proc_pattern}'")
            return returncode
