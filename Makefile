FC=gnu95
F2PY=f2py

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

all: EoS_Py

clean:
	rm -f *.o *.mod *.ld *.so

EoS_Py: \
	$(weaklib) \
	$(thornado) \
        EoS_Py.o
	$(FLINKER) $(FLAGS) -dynamiclib -shared -fPIC -o EoS_Py.so \
        $(weaklib) \
        $(thornado) \
        EoS_Py.o \
        $(LIBRARIES)

# makes executable.
# EoS_Py: \
# 	$(weaklib) \
# 	$(thornado) \
#         EoS_Py.o
# 	$(FLINKER) $(FLAGS) -o EoS_Py_$(MACHINE) \
#         $(weaklib) \
#         $(thornado) \
#         EoS_Py.o \
#         $(LIBRARIES)