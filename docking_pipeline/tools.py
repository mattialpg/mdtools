import subprocess

def run_in_terminal(command):
    return subprocess.run(command,
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,)
