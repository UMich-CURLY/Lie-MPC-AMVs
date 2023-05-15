#!/usr/bin/env python3
"""
Cleans an Octave script, possibly translating to Octave.

To specify output language, one of --octave and --matlab must be provided. If
--octave, the implicit defaults are --functions-in-place and --no-wrap-script.
If --matlab, the defaults are --functions-at-end and --wrap-script. Note that
these defaults can be overridden by specifying the option.

Any conflicting arguments will use the value that was passed last in the list
of arguments.
"""
import os
import re
import sys
import argparse

class YesNoAction(argparse.Action):
    """Action to support --toggle/--no-toggle pairs of arguments."""
    def __init__(self, option_strings, nargs=0, **kwargs):
        """Make sure argument pairs are valid and initialize."""
        if len(option_strings) != 2:
            raise ValueError("Must have exactly two option names!")
        elif nargs != 0:
            raise ValueError("nargs must be zero!")
        [self.YES, self.NO] = option_strings
        super().__init__(option_strings, nargs=nargs, **kwargs)
        
    def __call__(self, parser, namespace, values, option):
        """Checks whether option is yes or no string."""
        yesnoval = (True if option == self.YES else
                    (False if option == self.NO else None))
        setattr(namespace, self.dest, yesnoval)

[descr, doc] = __doc__.split("\n", 1)
usage = "%(prog)s --octave/--matlab [options] script out"
parser = argparse.ArgumentParser(description=descr, epilog=doc, usage=usage)
parser.add_argument("--octave", "--matlab", required=True,
                    help="Set target language", action=YesNoAction)
parser.add_argument("--functions-at-end", "--functions-in-place",
                    help="Whether to move function definitions to end",
                    action=YesNoAction)
parser.add_argument("--wrap-script", "--no-wrap-script",
                    help="Whether to wrap the example as a script",
                    action=YesNoAction)
parser.add_argument("script", help="Script file to transform")
parser.add_argument("out", help="Name for output file")

# List delimiters for blocks to omit from examples.
omitblocks = {"<--": "-->"}
commentchar = "%"

def cleanscript(name, read, write, wrapscript=True, funcsatend=False,
                matlab=True):
    """
    Reads the script read, cleans it, and writes to write.
    
    Options decide whether to wrap the script into a function (for older
    Matlab), whether to move functions to the end of the file, and whether to
    change from Octave to Matlab.    
    """
    subfuncs = []
    subfunc = None
    subfuncline = None
    omit = None
    omitlineno = None
    for (i, line) in enumerate(read):
        if i == 0:
            # Check whether script needs to be wrapped.
            wrapscript &= not line.startswith("function")
            if wrapscript:
                write.write("function %s()\n" % name)
        
        # Check for omit prefix.
        sline = line.lstrip()
        if sline.startswith(commentchar):
            sline = sline[len(commentchar):].lstrip()
            if omit is None:
                for check in omitblocks:
                    if sline.startswith(check):
                        omit = check
                        omitlineno = i
            else:
                if sline.startswith(omitblocks[omit]):
                    omit = None
                    omitlineno = None
                continue
        
        # Skip line if currently inside an omit block.
        if omit is not None:
            continue
        
        # Transform Octave-only syntax to Matlab.
        line = transform_source(line, matlab=matlab)
        
        # Break out subfunctions if necessary.
        if subfunc is None:
            if funcsatend and i != 0 and line.strip().startswith("function"):
                subfuncline = i
                subfunc = [line]
            else:
                write.write(line)
        else:
            subfunc.append(line)
        
        # Check if end of subfunction.
        if subfunc is not None:
            if line.strip().startswith("end%function"):
                subfuncs.append(subfunc)
                subfunc = None
                subfuncline = None
    
    # Make sure we closed any subfunction.
    if subfunc is not None:
        raise ValueError("No end%%function found for '%s' on line %d!"
                         % (subfunc[0].strip(), subfuncline + 1))
    
    # Make sure any omit block was ended.
    if omit is not None:
        raise ValueError("No '%s' found to close '%s' from line %d!"
                         % (omitblocks[omit], omit, omitlineno + 1))
    
    # Close function at bottom.
    if wrapscript:
        write.write("end%function\n")
    
    # Dump subfunctions at end.
    if len(subfuncs) > 0:
        write.write("\n\n% {0}\n% Subfunctions\n% {0}\n".format("*"*64))
        for subfunc in subfuncs:
            write.write("\n\n")
            for line in subfunc:
                #write.write("    ")
                write.write(line)

def transform_source(src, matlab=True):
    src = src.replace("endfunction", "end%function")
    src = re.sub("#+", "%", src) # For consistency.    
    if matlab:
        src = re.sub(r"randn?\('(?:seed|state)',\s*(.*)\)", r"rng(\1)", src)
        src = re.sub("end(for|if|switch|while)", "end", src)
        src = re.sub("%%+", "%", src) # Prevents weird things with %% in Matlab editor.
        src = re.sub("!=", "~=", src) # Matlab doesn't know !=.
        src = re.sub(r"^.*pkg\(.*\).*$", "", src) # Don't call pkg.
    return src

# Do main logic.
if __name__ == "__main__":
    # Decide arguments.
    args = vars(parser.parse_args(sys.argv[1:]))
    out = args["out"]
    outdir = os.path.dirname(out)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    scriptname = os.path.splitext(os.path.basename(out))[0]
    
    # Choose language and set defaults.
    matlab = not args["octave"]
    if matlab:
        defaults = {"wrap_script" : True, "functions_at_end" : True}
    else:
        defaults = {"wrap_script" : False, "functions_at_end" : False}
    for (name, default) in defaults.items():
        if args[name] is None:
            args[name] = default
    
    # Actually set options.    
    kwargs = dict()
    kwargs["matlab"] = matlab
    kwargs["wrapscript"] = args["wrap_script"]
    kwargs["funcsatend"] = args["functions_at_end"]
    
    # Call main function.
    with open(args["script"], "r") as read, open(out, "w") as write:
        cleanscript(scriptname, read, write, **kwargs)
    