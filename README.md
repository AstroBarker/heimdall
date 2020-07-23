# heimdall
Helpful EOS Investigative Multi-Dimensional thornAdo anaLysis tooL (HEIMDALL)

Complementary equation of state (EoS) analysis code for [thornado](https://github.com/endeve/thornado) [1](https://iopscience.iop.org/article/10.1088/1742-6596/1225/1/012014)[2](https://trace.tennessee.edu/cgi/viewcontent.cgi?article=3333&context=utk_chanhonoproj). 

_Heimdall_ contains FORTRAN routines to call _thornado_'s EoS routines. The main purpose is to obtain thermodynamic derivatives as interpolated from the 
tabulated nuclear EoS. These routines are called by a Julia wrapper to construct the derivatives and other quantities such as the eigenvector matrices for 
characteristic limiting. 

# Setup

_Heimdall_ is built on Julia, so that needs to be [installed](https://julialang.org/downloads/). After that, _heimdall_ just needs 
_thornado_'s environment variables to be set:

* `THORNADO_DIR` - path to this _THORNADO_ directory, e.g. export THORNADO_DIR=${HOME}/path/to/thornado
* `WEAKLIB_MODELS` - path to weaklib directory 

# Useage

Assuming that _thornado_ and _weaklib_ are set up correctly, you should just be able to `make`.
...

# TODO:

* Implement optional arguement to functions - "remove units"
 * May need a new "composite type" for holding data that is unit - aware?
 * Have to modify FORTRAN routines to accept unitless inputs

* Implement Multi-D stuff and things