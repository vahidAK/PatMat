import subprocess


def bgzip_and_tabix(filename: str) -> None:
    """Compress file with bgzip and create tabix index.

    Args:
        filename: Path to the file to compress and index.
            File will be compressed to filename.gz and indexed to filename.gz.tbi

    Raises:
        subprocess.CalledProcessError: If bgzip or tabix commands fail
    """
    run_command(f"bgzip -f {filename} && tabix -f -S 1 -p bed {filename}.gz")


def run_command(cmd: str) -> None:
    """Run a shell command and check for successful execution.

    Args:
        cmd: Shell command to execute

    Raises:
        subprocess.CalledProcessError: If the command returns non-zero exit status
    """
    subprocess.run(
        cmd,
        shell=True,
        check=True,
    )
