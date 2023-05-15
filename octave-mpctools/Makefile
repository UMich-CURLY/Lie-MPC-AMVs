# Constants.
TOPDIR := $(shell pwd)
DOC := doc
PKG := +mpctools
EXDIR := examples-src
UTIL := util

CLEANSCRIPT := $(UTIL)/cleanscript.py
GETMFILEDOC := $(UTIL)/getmfiledoc.py
MAKEZIP := $(UTIL)/makezip.py
LATEX2PDF := $(UTIL)/latex2pdf.py
DOSOURCECOMPARISON := $(UTIL)/dosourcecomparison.py
MATPLOTLIB := python3

# Source files. Note that the order of the files in DOCUMENTED_MFILES is the
# order they are listed in the documentation.
DOCUMENTATION_TEX := documentation.tex cheatsheet.tex introslides.tex comparison.tex install.tex
DOCUMENTED_MFILES := getCasadiFunc.m getCasadiIntegrator.m getLinearizedModel.m \
                     nmpc.m nmhe.m sstarg.m mpcplot.m spdinv.m rk4.m ekf.m c2d.m \
                     nlfilter.m collocweights.m getCasadiDAE.m version.m \

DOCUMENTED_CLASSFILES := ControlSolver.m MHESolver.m KalmanFilter.m

MFILES := $(PKG)/* import_mpctools.m

# List example files without directory. They must all be inside $(EXDIR).
EXAMPLES := \
    ballmaze.m \
    blockstructure.m \
    cstr.m \
    cstr_nmpc_nmhe.m \
    cstr_startup.m \
    daeexample.m \
    demcharge.m \
    econmpc.m \
    fishing.m \
    fullinformationexample.m \
    gainscheduling.m \
    mheexample.m \
    mhe_discretization.m \
    nmheexample.m \
    nmpcexample.m \
    onenorm.m \
    predatorprey.m \
    priorupdates.m \
    rocpenalty.m \
    rounding.m \
    runall.m \
    simplemhe.m \
    simplempc.m \
    softconstraints.m \
    springpendulummhe.m \
    sstargexample.m \
    timevaryingmhe.m \
    timevaryingmpc.m \
    vdposcillator.m \

EXTRA_FILES := README.md COPYING.txt

DIST_NAME := mpctools

# List example files (still without directory) that produce data files used in
# the documentation.
EXAMPLES_MAT_M := cstr_startup.m cstr.m vdposcillator.m ballmaze.m econmpc.m

# List Python/Matplotlib scripts found inside doc.
MATPLOTLIB_SRC := cstr_startup.py cstr.py vdposcillator.py ballmaze.py econmpc.py

# Path manipulation and file aggregation.
DOCUMENTATION_TEX := $(addprefix $(DOC)/, $(DOCUMENTATION_TEX))
DOCUMENTATION_PDF := $(DOCUMENTATION_TEX:.tex=.pdf)
TEX_SUFFIXES := .log .aux .synctex.gz .toc .nav .out .snm .vrb
TEX_MISC := $(foreach doc, $(DOCUMENTATION_TEX), $(addprefix $(basename $(doc)), $(TEX_SUFFIXES)))
DOCUMENTED_MFILES := $(addprefix $(PKG)/, $(DOCUMENTED_MFILES))

CLASS_DOC_MD := $(addprefix $(DOC)/, $(DOCUMENTED_CLASSFILES:.m=.md))
CLASS_DOC_TEX := $(CLASS_DOC_MD:.md=.tex)
CLASS_DOC_M := $(addprefix $(PKG)/, $(DOCUMENTED_CLASSFILES))

OCTAVE_EXAMPLES := $(addprefix examples-octave/, $(EXAMPLES))
MATLAB_EXAMPLES := $(addprefix examples-matlab/, $(EXAMPLES))

EXAMPLES_MAT_MAT := $(addprefix $(DOC)/, $(EXAMPLES_MAT_M:.m=.mat))
EXAMPLES_MAT_M := $(addprefix $(EXDIR)/, $(EXAMPLES_MAT_M))

MATPLOTLIB_SRC := $(addprefix $(DOC)/, $(MATPLOTLIB_SRC))
MATPLOTLIB_PDF := $(MATPLOTLIB_SRC:.py=.pdf)

DIST_FILES := $(MFILES) $(DOCUMENTATION_PDF) $(EXTRA_FILES) \
              $(OCTAVE_EXAMPLES) $(MATLAB_EXAMPLES)

DIST_ZIP := $(addsuffix .zip, $(DIST_NAME))

# Set default target.
.DEFAULT_GOAL := $(DIST_ZIPS)

# Recipes.
$(DOCUMENTATION_PDF) : %.pdf : %.tex
	@echo "Making $@."
	@$(LATEX2PDF) --dir $(DOC) $<

$(DOC)/mdoc.md : $(DOCUMENTED_MFILES)
	@echo "Making $@."
	@$(GETMFILEDOC) --out $@ $^

$(CLASS_DOC_MD) : $(DOC)/%.md : $(PKG)/%.m
	@echo "Making $@."
	@$(GETMFILEDOC) --class --out $@ $<

$(DOC)/mdoc.tex $(CLASS_DOC_TEX) : %.tex : %.md Makefile
	@echo "Making $@."
	@pandoc $< --from=markdown --to=latex -o $@

$(DIST_ZIP) : %.zip : $(DIST_FILES) Makefile
	@echo "Making $@."
	@$(MAKEZIP) --name $@ --root $(basename $@) $(DIST_FILES)
.PHONY : $(DIST_ZIP)
# DIST_ZIP depends on the current hg changeset ID, so it always needs to be
# out of date.

$(OCTAVE_EXAMPLES) : examples-octave/%.m : $(EXDIR)/%.m
	@$(CLEANSCRIPT) --octave $< $@

$(MATLAB_EXAMPLES) : examples-matlab/%.m : $(EXDIR)/%.m
	@$(CLEANSCRIPT) --matlab $< $@

$(EXAMPLES_MAT_MAT) : $(DOC)/%.mat : $(EXDIR)/%.m
	@echo "Making $@."
	@cd $(@D) && octave $(TOPDIR)/$<

$(MATPLOTLIB_PDF) : %.pdf : %.py
	@echo "Making $@."
	@cd $(<D) && $(MATPLOTLIB) $(<F)

dist : $(DIST_ZIP)
.PHONY : dist

examples : $(MATLAB_EXAMPLES) $(OCTAVE_EXAMPLES)
.PHONY : examples

# Cleanup rules.
clean :
	@echo "Cleaning up."
	@rm -f $(DOC)/mdoc.md $(DOC)/mdoc.tex $(CLASS_DOC_MD) $(CLASS_DOC_TEX)
	@rm -f $(MATPLOTLIB_PDF) $(EXAMPLES_MAT_MAT)
	@rm -f $(TEX_MISC)
	@rm -f $(DOC)/casadivsmtc.tex
.PHONY : clean

cleanexamples :
	@rm -rf examples-octave examples-matlab
.PHONY : cleanexamples

realclean : clean cleanexamples
	@rm -f $(DOCUMENTATION_PDF)
	@rm -f $(DIST_ZIP)
.PHONY : realclean

# Bitbucket upload rule.
UPLOAD_COMMAND := POST https://api.bitbucket.org/2.0/repositories/rawlings-group/octave-mpctools/downloads
upload : mpctools.zip
	echo -n "Enter bitbucket username: " && read bitbucketuser && curl -v -u $$bitbucketuser -X $(UPLOAD_COMMAND) -F files=@"$<"
.PHONY : upload

# Extra dependencies.
$(DOC)/documentation.pdf : $(DOC)/mdoc.tex $(CLASS_DOC_TEX)

$(DOC)/introslides.pdf : $(MATPLOTLIB_PDF)

$(DOC)/casadivsmtc.tex : $(EXDIR)/comparison_casadi.m $(EXDIR)/comparison_mtc.m
	@echo "Making $@."
	@$(DOSOURCECOMPARISON) $^ $@

$(DOC)/comparison.pdf : $(DOC)/casadivsmtc.tex

simplefigures = $(addprefix $(DOC)/, cstr_startup.pdf cstr.pdf vdposcillator.pdf ballmaze.pdf econmpc.pdf)
$(simplefigures) : %.pdf : %.mat