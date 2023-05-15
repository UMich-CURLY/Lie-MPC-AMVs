#!/usr/bin/env python3
"""
Gets documentation from the first comment string of an Octave/Matlab m-file.
"""

import sys
import re
import argparse
import os
import time

parser = argparse.ArgumentParser(add_help=False, description=__doc__)
parser.add_argument("--class", help="Treat file as a classdef file",
                    action="store_true")
parser.add_argument("--out", help="Name of output file.")
parser.add_argument("mfile", help="Source m-files", nargs="+")

FUNCTION_RE = r"^function (?:.* = )?([A-Za-z]\w*)\(.*\)\s*$"
CLASS_RE = r"^classdef ([A-Za-z]\w*)"

class OneGroupRE:
    """Holds a regex with one matching group."""
    def __init__(self, regex):
        """Initializes the object with the current regex."""
        self.__regex = re.compile(regex)
    
    def __call__(self, line):
        """Searches the current line for the given regex."""
        match = self.__regex.search(line)
        if match is not None:
            match = match.group(1)
        return match

class ClassFileDoc:
    """Class for reading documentation from a classdef m-file."""
    def __init__(self, outfile=None):
        """Opens outfile for writing and prepares to process lines."""
        if outfile is None:
            self.__outfile = sys.stdout
        else:
            self.__outfile = outfile
        self.write = self.__outfile.write
        self.inclass = False
        self.inmethods = False
        self.infunc = False
        self.indoc = False
        self.funcstart = OneGroupRE(FUNCTION_RE)
        self.classstart = OneGroupRE(CLASS_RE)
        self.funcprefix = ""
        self.classname = None
        self.funcname = None
        self.comment = "%"
        self.nfunctions = 0
    
    def process(self, line):
        """Processes the current line, writing to the file if necessary."""
        # Strip leading spaces and write line (logic inside writedoc).        
        line = line.lstrip()
        self.writedoc(line)        

        # Look for beginnings or endings of various blocks.        
        if self.inclass:
            if self.inmethods:
                if self.infunc:
                    self.checkendfunc(line)
                else:
                    self.checkstartfunc(line)
                self.checkendmethods(line)
            else:
                self.checkstartmethods(line)
            self.checkendclass(line)
        else:
            self.checkstartclass(line)
    
    def checkstartclass(self, line):
        """Looks for the beginning of a classdef."""
        found = False
        classname = self.classstart(line)
        if classname is not None:
            self.inclass = True
            self.classname = classname
            self.funcprefix = classname + "."
            self.funcname = classname
            self.indoc = True
            found = True
        return found
    
    def checkendclass(self, line):
        """Looks for the end of a classdef."""
        found = False
        if line.startswith("end%classdef"):
            self.inclass = False            
            self.classname = None
            self.funcprefix = ""
            found = True
        return found
    
    def checkstartmethods(self, line):
        """Checks if the current line is the start of the methods block."""
        found = False    
        if line.startswith("methods (Access=public)"):
            self.inmethods = True
            found = True
        return found
    
    def checkendmethods(self, line):
        """Checks if the current line ends the methods block."""
        found = False
        if line.startswith("end%methods"):
            self.inmethods = False
            found = True
        return found
    
    def checkstartfunc(self, line):
        """Checks the current line for the start of a function."""
        found = False
        func = self.funcstart(line)
        if func is not None and func != self.classname:
            self.nfunctions += 1
            self.infunc = True
            self.indoc = True
            self.funcname = self.funcprefix + func
            found = True
        return found

    def checkendfunc(self, line):
        """Checks if the current line ends the function block."""
        found = False
        if line.startswith("end%function"):
            self.infunc = False
            found = True
        return found
        
    def writedoc(self, line):
        """Checks the current line for a comment and prints if found."""
        wrote = False
        if self.indoc:
            if line.startswith(self.comment):
                if self.funcname is not None:
                    self.write("## %s ##\n\n" % self.funcname)
                    self.funcname = None
                line = line[len(self.comment):]
                if line.startswith(" "):
                    line = line[1:]
                self.write(line)
                wrote = True
            else:
                self.write("\n")
                self.indoc = False
        return wrote

class MFileDoc(ClassFileDoc):
    """Class for printing the documentation for a single m-file."""
    def __init__(self, outfile=None):
        """Initializes the object."""
        super().__init__(outfile)
        self.inclass = True
        self.inmethods = True
    
    def checkendfunc(self, line):
        """Stop looking after first function is over."""
        found = super().checkendfunc(line)
        if found:
            self.inmethods = False
            self.inclass = False
        return found
    
    def noop(*args):
        """Does nothing."""
        pass
    
    checkstartmethods = noop
    checkendmethods = noop
    checkstartclass = noop
    checkendclass = noop

if __name__ == "__main__":
    # Read arguments.
    args = vars(parser.parse_args(sys.argv[1:]))
    
    # Choose documentation object.
    Doc = ClassFileDoc if args["class"] else MFileDoc
    
    # Open output file.
    if args["out"] is None:
        outfile = sys.stdout
    else:
        outfile = open(args["out"], "w")
    outfile.write("<!--- Written by getmfiledoc.py on %s -->\n\n"
                  % time.strftime("%Y-%m-%d %H:%M:%S"))
    
    # Loop through files.
    error = None
    try:
        for mfile in args["mfile"]:
            doc = Doc(outfile)
            with open(mfile, "r") as read:
                for line in read:
                    doc.process(line)
            if doc.nfunctions == 0:
                raise IOError("No functions found in '%s'!" % mfile)
    except Exception as err:
        # Something went wrong. Save exception for later
        error = err
    finally:
        # Close open file
        if outfile is not sys.stdout:
            outfile.close()
            if error is not None:
                # Give an old timestamp so Make thinks its out of date.
                os.utime(args["out"], (1, 1))
        
        # Re-raise exception.
        if error is not None:
            raise RuntimeError("Error processing file {}!"
                               .format(outfile)) from error
                