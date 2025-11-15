"""
CEDAR: exploration of the tree space in a maximum likelihood framework
Running Raxml for tree space exploration
Adapted from https://github.com/Neclow/phylo2vec/blob/main/py-phylo2vec/phylo2vec/opt/_hc_losses.py
Accessed on 2025.09.30
"""

__author__ = "Cedric Chauve"
__credits__ = ["https://github.com/Neclow/phylo2vec/tree/main"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Cedric Chauve"
__email__ = "cedric.chauve@sfu.ca"
__status__ = "Release"

import os
import re
import subprocess
import sys
from pathlib import PurePosixPath

# Test if the current platform is Windows or not
IS_WINDOWS = sys.platform.startswith("win")
# Regex for a negative float
NEG_FLOAT_PATTERN = re.compile(r"-\d+.\d+")

def raxml_loss(
    fasta_path,
    tree,
    model,
    tree_folder_path,
    outfile="tmp.tree",
    **kwargs,
):
    """Compute loss for a given tree via RaXML-NG.

    Parameters
    ----------
    fasta_path : str
        Path to fasta file
    tree: TreeVec object
    model : str
        DNA/AA substitution model
    tree_folder_path : str
        Path to a folder which will contain all intermediary and best trees
    outfile : str, optional
        Path to a temporary tree written in Newick format, by default 'tmp.tree'

    Returns
    -------
    float
        Negative log-likelihood computed using RaXML-NG
    """
    newick = tree.treevec2newick()

    with open(
        os.path.join(tree_folder_path, outfile), "w", encoding="utf-8"
    ) as nw_file:
        nw_file.write(newick)

    return exec_raxml_ng(
        tree_path=str(
            PurePosixPath(tree_folder_path.replace("C:", "/mnt/c/"), outfile)
        ),
        fasta_path=str(PurePosixPath(fasta_path.replace("C:", "/mnt/c"))),
        model=model,
        **kwargs,
    )


def exec_raxml_ng(fasta_path, tree_path, model, cmd="raxml-ng", no_files=True):
    """Optimize branch lengths and free model parameters on a fixed topology
    using RaxML-NG (https://github.com/amkozlov/raxml-ng)

    Parameters
    ----------
    fasta_path : str
        Path to FASTA file (MSA)
    tree_path : str
        Path to tree file (Newick representation of the tree)
    model : str
        DNA evolution model
    cmd : str, optional
        Location of the RAxML-nG executable, by default "raxml-ng"
    no_files : bool, optional
        If True, add the "nofiles" option to raxml

    Returns
    -------
    float
        Negative log-likelihood after optimization
    """
    commands = [
        cmd,
        "--evaluate",
        "--msa",
        fasta_path,
        "--tree",
        tree_path,
        "--model",
        model,
        "--brlen",
        "scaled",
        "--log",
        "RESULT",
        "--threads",
        "1",
    ]

    if no_files:
        commands.append("--nofiles")

    if IS_WINDOWS:
        commands.insert(0, "wsl")  # Use Windows Subsystem for Linux
    else:
        commands = " ".join(commands)  # For Linux

    try:
        output = subprocess.run(
            commands, capture_output=True, check=True, shell=not IS_WINDOWS
        )
    except subprocess.CalledProcessError as _:
        # pylint: disable=subprocess-run-check
        output = subprocess.run(commands, capture_output=True, shell=not IS_WINDOWS)
        # pylint: enable=subprocess-run-check

        raise RuntimeError(output) from _

    stdout = output.stdout.decode("ascii")

    lik_line = [
        line for line in stdout.split("\n") if line.startswith("Final LogLikelihood")
    ][0]

    nll = -1 * float(re.findall(NEG_FLOAT_PATTERN, lik_line)[0])

    return nll
