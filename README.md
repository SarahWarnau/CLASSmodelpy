CLASS - Python
-------

This adapted version of the CLASS model has several additional functionalities. It was designed to simulate Technology-Enhanced Atmospheric Moistening scenarios (1), however, the atmospheric profile, sea surface, and distance instead of time modules can be used for other purposes as well.

# Atmospheric profile
With the atmospheric profile switch `sw_ap=True` CLASS will run in a different mode, with a discretized vertical grid. This has been implemented to facilitate running the model with non linear Free Atmospheric slopes.

As an input, a `pandas.DataFrame` needs to be provided which should at least contain a column labeled `'z'` containing the atmospheric levels in meters, and optionally vertical profiles of `'theta'`, `'q'`, `'u'`, and `'v'`. While running the model, the vertical profile will be updated with the new mixed layer values at the end of every integration step. Then, the profile is stored in a NetCDF4 file, which is called `out_NetCDF`.

During the run, free atmosphere slopes of the provided variables will be derived from the values above the ABL top. To make sure the atmospheric profile and the model variables are corresponding, the $\Delta \phi$ values are derived from the updated profile as well.

# Spray evaporation
With the spray switch `sw_sp=True` evaporation of sprayed water droplets into the ABL can be simulated. At the provided moment `tspray` the wet bulb potential temperature and specific humidity will be calculated at the spray height `zspray`. The sprayed layer will become the new mixed layer with the wet bulb potential temperature and specific humidity as mixed layer values. For this module to work, the atmospheric profile module needs to be active.

# Sea surface
With the sea surface switch `sw_ss=True` the surface has a constant temperature given by the `SST` input with corresponding saturation specific humidity. Roughness lengths for momemtum and scalars specific to the sea need to be specified when the surface layer switch is on.

# Solar evaporation surface
With the solar evaporation surface switch `sw_so=True` the surface fluxes are calculated for a black (no shortwave reflection), well insulated (no ground heat flux), freely evaporating (no surface resistance) surface. At every timestep, the surface temperature and corresponding saturation specific humidity are calculated from the net incoming radiation `Qin`, which can be constant, or variable when the radiation switch is on. Roughness lengths for momemtum and scalars specific to the solar technology need to be specified when the surface layer switch is on.

# Advecting an ABL over different surfaces
With the distance instead of time integration switch `sw_x=True` CLASS will run for a provided distance instead of a provided time. This has been implemented to facilitate running the model over different types of surfaces in one run. The velocity of ABL (column velocity `col_vel`) needs to be provided and in combination with the `numpy.array` `x` the `runtime` and `dt` are determined. From the `numpy.array` `X`, the codes are used to change the surface modules where 0=sea surface, 1=land surface, 2=solar evaporator surface.

# Example setup
`runmodel_atmprof.py` contains an example of how to use all these functionalities, and retrieve and plot the outputs.

NOTE: the `Makefile` for compiling `ribtol.cpp` (fast Bulk Richardson - Obukhov length iterative solver) is not required to run CLASS. By default a (somewhat slower) Python solver is used, and the fast c++ version is only interesting for computationally expensive experiments.


# References
(1) Warnau, S. N.; Theeuwen, J. J. E.; Sadeghi, G.; Benedict, I.; Hamelers, B. H. V. M.; van Heerwaarden, C. C. Technology-Enhanced Atmospheric Moistening (TEAM) for More Precipitation: A Perspective. Environ. Sci. Technol. 2026, 60 (2), 1612–1620. https://doi.org/10.1021/acs.est.5c06428.
