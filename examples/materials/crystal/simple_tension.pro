############################################################################
#  Crystal Plasticity Example - Simple Tension Test                       #
#                                                                          #
#  Description: Single crystal aluminum under uniaxial tension             #
#               Demonstrates crystal plasticity with {110}<111> slip      #
#               systems and rate-dependent behavior.                       #
#                                                                          #
#  Usage:       pyfem simple_tension.pro                                   #
############################################################################

input = "simple_tension.dat";

ContElem = 
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "Crystal";
    
    # Elastic properties (cubic crystal - aluminum)
    c11 = 108.0e3;   # MPa
    c12 = 62.0e3;    # MPa
    c44 = 28.0e3;    # MPa
    
    # Alternative: use isotropic elastic properties
    # E  = 70.0e3;   # MPa
    # nu = 0.33;
    
    # Slip system sets
    nsets = 1;
    
    # Set 1: {110}<111> slip systems (FCC crystal)
    plane1 = [1.0, 1.0, 0.0];   # {110} plane normal
    dir1   = [1.0, 1.0, 1.0];   # <111> slip direction
    
    # Crystal orientation (identity - no rotation)
    orient1_local  = [1.0, 0.0, 0.0];
    orient1_global = [1.0, 0.0, 0.0];
    orient2_local  = [0.0, 1.0, 0.0];
    orient2_global = [0.0, 1.0, 0.0];
    
    # Viscoplasticity parameters
    gamma0 = 0.001;   # Reference shear strain rate (1/s)
    m      = 20.0;    # Rate sensitivity exponent
    
    # Hardening parameters
    tau0 = 16.0;      # Initial CRSS (MPa)
    h0   = 180.0;     # Initial hardening modulus (MPa)
    taus = 148.0;     # Saturation stress (MPa)
    a    = 2.25;      # Hardening exponent
    qab  = 1.4;       # Latent hardening ratio
    
    # Integration parameters
    theta   = 0.5;    # Implicit integration parameter
    nlgeom  = false;  # Use small strain formulation
    maxiter = 20;     # Max Newton-Raphson iterations
    tol     = 1.0e-6; # Convergence tolerance
  };
};

solver =
{
  type = "NonlinearSolver";
  
  maxCycle = 50;
  dtime    = 0.1;   # Time increment (important for rate-dependent plasticity)
};

outputModules = ["vtk", "GraphWriter", "OutputWriter"];

vtk =
{
  type     = "MeshWriter";
  prefix   = "simple_tension";
  interval = 1;
};

GraphWriter =
{
  type     = "GraphWriter";
  onScreen = true;
  
  columns = ["time", "disp", "force"];
  
  disp =
  {
    type = "state";
    node = 5;
    dof  = "u";
  };
  
  force =
  {
    type = "fint";
    node = 5;
    dof  = "u";
  };
};

OutputWriter =
{
  type     = "OutputWriter";
  onScreen = true;
};
