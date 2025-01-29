import subprocess


def bgzip_and_tabix(filename):
    run_command(f"bgzip -f {filename} && tabix -f -S 1 -p bed {filename}.gz")


def run_command(cmd):
    subprocess.run(
        cmd,
        shell=True,
        check=True,
    )
