#!/usr/bin/env python3
"""Compares two source files that have been split into chunks."""
import sys
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("left", help="file to show on the left")
parser.add_argument("right", help="file to show on the right")
parser.add_argument("tex", help="name of tex file to write")

ENDCHUNK = set(["#<<ENDCHUNK>>","%<<ENDCHUNK>>"])

def comparison(infiles, outfile, styles=None):
    """
    Opens finfiles and spits out comparison tex file.
    """
    
    files = [open(f, "r") for f in infiles]
    out = open(outfile, "w")
    if styles is None:
        styles = [""]*len(infiles)
    
    keepgoing = True
    while keepgoing:
        keepgoing = False
        out.write(r"\sepline" + "\n\n")
        out.write(r"\begin{parcolumns}{%d}" % len(files) + "\n")
        for (f, style) in zip(files, styles):
            lines = []
            stopchunk = False
            while not stopchunk:
                line = f.readline()
                if line.strip() in ENDCHUNK or len(line) == 0:
                    stopchunk = True
                    if len(line) > 0:
                        keepgoing = True
                else:
                    lines.append(line)
            outputlines(out,lines, style)
        out.write(r"\end{parcolumns}" + "\n\n")
    
    # Clean up.    
    out.close()
    for f in files:
        f.close()
        
        
def outputlines(f, lines, style=""):
    """
    Clean up lines and then write a colchunk environment to f.
    """
    if len(lines) > 0:
        # First remove any leading or trailing blank lines.        
        nhead = 0
        while nhead < len(lines) and len(lines[nhead].strip()) == 0:
            nhead += 1
        ntail = len(lines) - 1
        while ntail > 0 and len(lines[ntail].strip()) == 0:
            ntail -= 1
        lines = lines[nhead:ntail + 1]
        
        # Now add lstlisting environment.
        lines.insert(0,r"\begin{lstlisting}" + style + "\n")
        lines.append(r"\end{lstlisting}" + "\n")

    # Write colchunk begin and end.
    lines.insert(0,r"\colchunk{" + "\n")
    lines.append("}\n")    
    
    # Finally write.
    f.writelines(lines)

if __name__ == "__main__":
    # Read files from command line.
    args = vars(parser.parse_args(sys.argv[1:]))
    comparison([args["left"], args["right"]], args["tex"])
   