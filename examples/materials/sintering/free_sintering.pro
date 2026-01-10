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
# Example: Free sintering of ceramic sample
#
# This example demonstrates viscous densification during free sintering
# of a ceramic powder compact. The model captures:
#
# - Progressive densification driven by sintering stress
# - Temperature-dependent viscous flow
# - Evolution of material stiffness with density
# - Time-dependent shrinkage
#
# The Skorohod-Olevsky model describes the viscous sintering process.
#

input = "free_sintering.dat";

Continuum =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "SkorohodOlevsky";
    
    # Viscosity parameters
    eta0 = 1.0e12;         # Pa·s - Reference viscosity at reference temp
    Q = 500000.0;          # J/mol - Activation energy for sintering
    R = 8.314;             # J/(mol·K) - Gas constant
    T = 1600.0;            # K - Sintering temperature
    
    # Density parameters
    rho0 = 0.6;            # Initial relative density (60% of theoretical)
    
    # Sintering stress
    sigma_sint = 1.0e6;    # Pa - Capillary-driven sintering stress
    
    # Viscosity exponents (Skorohod functions)
    n_vol = 2.0;           # Volumetric viscosity exponent
    n_shear = 1.0;         # Shear viscosity exponent
  };
};

solver =
{
  type = "NonlinearSolver";
  
  maxCycle = 10;
  
  # Time stepping for sintering process
  # Total time: 100 seconds of simulated sintering
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

  columns = ["time", "density", "shrinkage", "stress"];

  time =
  {
    type = "time";
  };

  density =
  {
    type = "Density";
    node = 9;
  };
 
  shrinkage =
  {
    type = "state";
    node = 17;
    dof = 'u';
    factor = -1.0;  # Negative because shrinkage is negative displacement
  };

  stress =
  {
    type = "S11";
    node = 9;
  };
};
