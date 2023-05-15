#!/usr/bin/env python3
"""Adds all distribution files to a zip file for upload to Bitbucket."""

import argparse
import os
import sys
import zipfile
import subprocess

# Command-line arguments.
parser = argparse.ArgumentParser(description=__doc__, add_help=False)
parser.add_argument("--help", help="print this help", action="help")
parser.add_argument("--dist", help="which files to distribute",
                    choices=["octave", "matlab", "both"], default="both")
parser.add_argument("--root", help="name for root folder in zip file")
parser.add_argument("--name", help="specify name for zip file",
                    default="mpctools.zip")
parser.add_argument("files", default=[], nargs="+",
                    help="Files to include")
kwargs = vars(parser.parse_args(sys.argv[1:]))


# Get files and name of zip file, and decide if there should be a root folder.
files = set(kwargs["files"])
zipname = kwargs["name"]
root = kwargs["root"]
if root is None:
    root = ""

# Filter out example files depending on distribution.
filters = {
    "both" : lambda s : True,
    "octave" : lambda s : not s.startswith("examples-matlab"),
    "matlab" : lambda s : not s.startswith("examples-octave"),
}
files = list(filter(filters[kwargs["dist"]], files))
  
# Decide whether to include README.
README = "README.md"
includereadme = (README in files)
if includereadme:
    files.remove(README)

# Helper functions to translate the version file.
def getgitid():
    """Returns the current git changeset ID as a string."""
    gitid = subprocess.check_output(["git", "rev-parse", "HEAD"])
    return gitid.decode().strip()

def cleanversion(infile):
    """
    Cleans +mpctools/version.m for distribution. Returns cleaned string.
    """
    skip = False
    gitid = getgitid()
    lines = []
    with open(infile, "r") as read:
        for line in read:
            if skip:
                if line.startswith("% -->"):
                    skip = False
            elif line.startswith("% <--"):
                skip = True
            else:
                lines.append(line.replace("%#CHANGESET_ID%", gitid))
    return "".join(lines)

# Now add files.
with zipfile.ZipFile(zipname, "w", zipfile.ZIP_DEFLATED) as z:
    # Recurse through VI directories.
    for fl in files:
        readfile = fl
        writefile = os.path.join(root, fl)
        if readfile == "+mpctools/version.m":
            lines = cleanversion(readfile)
            z.writestr(writefile, lines)
        else:
            z.write(readfile, writefile)    
    # Also add readme with txt extension to play nice with Windows.
    if includereadme:
        z.write(README, os.path.join(root, os.path.splitext(README)[0]
                                               + ".txt"))
    print("Wrote zip file '%s'." % z.filename)