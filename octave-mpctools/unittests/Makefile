# Makefile to help with multiple versions of Casadi.
CASADI_VERSIONS := 3.1.1   3.2.3   3.3.0   3.4.0   3.4.3

CASADI_ZIPS := $(patsubst %, casadi-%.tar.gz, $(CASADI_VERSIONS))
CASADI_DIRS := $(patsubst %, casadi-%, $(CASADI_VERSIONS))
CASADI_MEX := $(addsuffix /casadiMEX.mex, $(CASADI_DIRS))

# Casadi versions.
$(CASADI_MEX) : casadi-%/casadiMEX.mex : casadi-%.tar.gz
	@echo "Unzipping $<"
	@rm -rf $(@D)
	@mkdir $(@D)
	@cd $(@D) && tar -x -f ../$<
	@touch $@

# Downloads for zipped archives.
SOURCEFORGE := https://sourceforge.net/projects/casadi/files/CasADi
GITHUB := https://github.com/casadi/casadi/releases/download/

casadi-3.1.1.tar.gz casadi-3.2.3.tar.gz : casadi-%.tar.gz :
	@echo "Downloading $@"
	@wget -O $@ "$(SOURCEFORGE)/$*/linux/casadi-octave-v$*.tar.gz"

casadi-3.3.0.tar.gz : casadi-%.tar.gz :
	@echo "Downloading $@"
	@wget -O $@ "$(GITHUB)/$*/linux__casadi-octave-v$*.tar.gz"

casadi-3.4.0.tar.gz casadi-3.4.3.tar.gz : casadi-%.tar.gz :
	@echo "Downloading $@"
	@wget -O $@ "$(GITHUB)/$*/casadi-linux-octave-v$*.tar.gz"

# Set version.
$(CASADI_VERSIONS) : % : casadi-%/casadiMEX.mex
	@echo "Using CasADi version $*"
	@ln -sfT $(<D) casadi
.PHONY : $(CASADI_VERSIONS)

# Run all tests.
tests : $(CASADI_MEX)
	@./rununittests $(CASADI_VERSIONS)
.PHONY : alltests

# Cleanup rules.
clean :
	@echo "Cleaning up."
	@rm -rf $(CASADI_DIRS) casadi
realclean : clean
	@rm -rf $(CASADI_ZIPS)
.PHONY : clean realclean