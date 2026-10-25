############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  The code is written by J.J.C. Remmers, C.V. Verhoosel and R. de Borst.  #
#                                                                          #
#  The latest stable version can be downloaded from the web-site:          #
#     http://www.wiley.com/go/deborst                                      #
#                                                                          #
#  A github repository, with the most up to date version of the code,      #
#  can be found here:                                                      #
#     https://github.com/jjcremmers/PyFEM                                  #
#                                                                          #
#  The code is open source and intended for educational and scientific     #
#  purposes only. If you use PyFEM in your research, the developers would  #
#  be grateful if you could cite the book.                                 #  
#                                                                          #
#  Disclaimer:                                                             #
#  The authors reserve all rights but do not guarantee that the code is    #
#  free from errors. Furthermore, the authors shall not be liable in any   #
#  event caused by the use of the program.                                 #
############################################################################

#
# Example: Pressure-assisted sintering (hot pressing)
#
# This example demonstrates viscous densification during hot pressing
# of a ceramic powder compact. External pressure accelerates densification
# compared to free sintering.
#
# The model captures:
# - Accelerated densification under applied pressure
# - Combined effects of sintering stress and external load
# - Time-dependent material response
# - Evolution towards full density
#
# The Skorohod-Olevsky model describes pressure-assisted sintering.
#

input = "pressure_sintering.dat";

Continuum =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "SkorohodOlevsky";
    
    # Viscosity parameters
    eta0 = 1.0e12;         # Pa·s - Reference viscosity
    Q = 500000.0;          # J/mol - Activation energy
    R = 8.314;             # J/(mol·K) - Gas constant
    T = 1600.0;            # K - Hot pressing temperature
    
    # Density parameters
    rho0 = 0.6;            # Initial relative density
    
    # Sintering stress
    sigma_sint = 1.0e6;    # Pa - Sintering stress
    
    # Viscosity exponents
    n_vol = 2.0;           # Volumetric exponent
    n_shear = 1.0;         # Shear exponent
  };
};

solver =
{
  type = "NonlinearSolver";
  
  maxCycle = 10;
  
  # Time stepping for hot pressing
  tmax = 100.0;
  dtime = 2.0;
};

outputModules = ["vtk", "GraphWriter"];

vtk =
{
  type = "MeshWriter";
  interval = 5;
};

GraphWriter =
{
  onScreen = true;

  columns = ["time", "density", "height", "pressure"];

  time =
  {
    type = "time";
  };

  density =
  {
    type = "Density";
    node = 13;
  };
 
  height =
  {
    type = "state";
    node = 23;
    dof = 'v';
  };

  pressure =
  {
    type = "S22";
    node = 13;
  };
};
