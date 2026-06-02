# nlmixr2 6.0.1

- Fix LTO violation as requested by CRAN by adding
  -DARMA_DONT_USE_OPENMP to PKG_CXXFLAGS in src/Makevars.in

- Require rxode2 5.1.2 which has the fixed M1-san issues observed
  here.

- Guard against null pointer arithmetic in inner.cpp
