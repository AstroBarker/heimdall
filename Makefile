MACHINE = $(THORNADO_MACHINE)

OPT_LEVEL = DEBUG
FLAGS     = $(FLAGS_$(OPT_LEVEL))

MICROPHYSICS   = WEAKLIB
#GRAVITY_SOLVER = POSEIDON_NEWTON

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles

all: EoS_jl

EoS_jl: \
	$(weaklib) \
	$(thornado) \
        EoS_jl.o
	$(FLINKER) $(FLAGS) -dynamiclib -shared -fPIC -o EoS_jl.so \
        $(weaklib) \
        $(thornado) \
        EoS_jl.o \
        $(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

realclean:
	rm -f *.o *.mod *.ld *.so

# makes executable.
# EoS_jl: \
# 	$(weaklib) \
# 	$(thornado) \
#         EoS_jl.o
# 	$(FLINKER) $(FLAGS) -o EoS_jl_$(MACHINE) \
#         $(weaklib) \
#         $(thornado) \
#         EoS_jl.o \
#         $(LIBRARIES)