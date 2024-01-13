# Conformal Cubic Atmospheric Model

The Conformal Cubic Atmospheric Model (CCAM) is a global stretched grid climate
model used for weather and climate research.  CCAM is typically used for downscaling General Circulation Models for regional projections and tailored weather
forecasting.  CCAM's variable resolution grid can be multiply nested for
downscaling without lateral boundaries.  The model also supports various
Earth system components including an implicitly coupled ocean model, prognostic
aerosols and terrestial carbon cycle.


## Website

For documentation, see our website at
[https://confluence.csiro.au/display/CCAM/CCAM]


## Dependencies

CCAM requires the NetCDF and Message Passing Interface (MPI) libraries.
Although CCAM is written in fortran, it links with the NetCDF C library.


## Building CCAM

To build CCAM with the intel, gnu and cray fortran compiler use

```
make
make GFORTRAN=yes
make CRAY=yes
```

CCAM also has limited support for GPUs and OpenMP, which can be enabled with

```
make GPU=yes OMP=yes
```

Debugging is also enabled with

```
make TEST=yes
```





