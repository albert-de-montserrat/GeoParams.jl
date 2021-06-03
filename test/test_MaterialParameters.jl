using Test
using GeoParams


# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1000km);

# Define a linear viscous creep law
x1      =   GeoParams.MaterialParameters.LinearViscous(eta=1e18Pa*s)
@test x1.eta.val == 1e18Pa*s

x1_ND   =   GeoParams.MaterialParameters.LinearViscous(eta=1e18Pa*s)
Nondimensionalize!(x1_ND,CharUnits_GEO)        # check that we can nondimensionalize all entries within the struct
@test x1_ND.eta*1.0==0.1                            

# perform a computation with the viscous creep laws 
@test CreepLaw_EpsII(1e6*Pa, x1)==5.0e-13/s             # dimensional input       
@test CreepLaw_EpsII(1e0, x1_ND)==5.0                   # non-dimensional
@test CreepLaw_EpsII([1e0; 2.0], x1_ND)==[5.0; 10.0]    # vector input

#x2  =   GeoParams.MaterialParameters.LinearViscous(eta=1e3)

# Define a struct for a first phase
#Phase = MaterialParams(Name="test1",
#                CreepLaws=(GeoParams.MaterialParameters.LinearViscous, GeoParams.MaterialParameters.LinearViscous(eta=1e20)),
#                Density = 1000kg/m^3);

# Evaluate viscous strainrate, given a stress

