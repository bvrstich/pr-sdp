###############################################################################
#
#  Makefile template for the sources
#
###############################################################################

# -----------------------------------------------------------------------------
#   Sources for all modules
# -----------------------------------------------------------------------------
BINNAME = sdp
CPPSRC	= sdp.cpp\
            TPM.cpp\
            Matrix.cpp\
            SPM.cpp\
            PHM.cpp\
            DPM.cpp\
            SUP/SUP_PQ.cpp\
            SUP/SUP_PQG.cpp\
            SUP/SUP_PQGT1.cpp\
            EIG/EIG_PQ.cpp\
            EIG/EIG_PQG.cpp\
            EIG/EIG_PQGT1.cpp

OBJ	= $(CPPSRC:.cpp=.o)

# -----------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -----------------------------------------------------------------------------

BRIGHT_ROOT= .

MKL_ROOT = /opt/intel/Compiler/11.1/069/mkl

INCLUDE = -I${MKL_ROOT}/include

LIBS= -L${MKL_ROOT}/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

CC	= icc
CXX	= icpc

# -----------------------------------------------------------------------------
#   Compiler & Linker flags
# -----------------------------------------------------------------------------
CFLAGS	= $(INCLUDE) -g -Wall -O2 -ipo -xSSE4.2
LDFLAGS	= $(LIBS) -g -Wall -O2 -xSSE4.2

# =============================================================================
#   Targets & Rules
# =============================================================================
all:
	@echo
	@echo '  +++ Building $(BINNAME)...'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQG"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built successfully!'; \
	   echo; \
	 fi

#------------------------------------------------------------------------------
#  Compile with only P and Q conditions activated
#------------------------------------------------------------------------------

PQ:
	@echo
	@echo '  +++ Building $(BINNAME) with P and Q conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQ"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P and Q conditions successfully!'; \
	   echo; \
	 fi

PQG:
	@echo
	@echo '  +++ Building $(BINNAME) with P, Q and G conditions active'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQG"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P, Q and G conditions successfully!'; \
	   echo; \
	 fi

PQGT1:
	@echo
	@echo '  +++ Building $(BINNAME) with P, Q, G and T_1 conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQGT1"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P, Q, G and T_1 conditions successfully!'; \
	   echo; \
	 fi

# -----------------------------------------------------------------------------
#   The default way to compile all source modules
# -----------------------------------------------------------------------------
%.o:	%.for Makefile
	@echo; echo "Compiling $(@:.o=.for) ..."
	$(FF) -c $(FFLAGS) $(SFLAGS) $(@:.o=.for) -o $@

%.o:	%.c Makefile
	@echo; echo "Compiling $(@:.o=.c) ..."
	$(CC) -c $(CFLAGS) $(SFLAGS) $(@:.o=.c) -o $@

%.o:	%.cpp Makefile
	@echo; echo "Compiling $(@:.o=.cpp) ..."
	$(CXX) -c $(CFLAGS) $(SFLAGS) $(DEFS) $(@:.o=.cpp) -o $@


# -----------------------------------------------------------------------------
#   Link everything together
# -----------------------------------------------------------------------------
$(BRIGHT_ROOT)/$(BINNAME):	Makefile $(OBJ) 
	@echo; echo "Linker: creating $(BRIGHT_ROOT)/$(BINNAME) ..."
	$(CXX) $(LDFLAGS) $(SFLAGS) -o $(BRIGHT_ROOT)/$(BINNAME) $(OBJ) $(LIBS)

# -----------------------------------------------------------------------------
#   Create everything newly from scratch
# -----------------------------------------------------------------------------
new:	clean all

# -----------------------------------------------------------------------------
#   Clean up all object files
# -----------------------------------------------------------------------------
clean:
	@echo -n '  +++ Cleaning all object files ... '
	@echo -n $(OBJ)
	@rm -f $(OBJ)
	@echo 'Done.'

doc:
	doxygen doc-config


# ====================== End of file 'Makefile.in' ========================== #
