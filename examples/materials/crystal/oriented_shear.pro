############################################################################
#  Crystal Plasticity Example - Shear Test with Oriented Crystal          #
#                                                                          #
#  Description: Single crystal with 45-degree orientation under shear     #
#               Multiple slip system sets to capture realistic behavior    #
#                                                                          #
#  Usage:       pyfem oriented_shear.pro                                   #
############################################################################

input = "oriented_shear.dat";

ContElem = 
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "Crystal";
    
    # Elastic properties (cubic crystal - copper)
    c11 = 168.0e3;   # MPa
    c12 = 121.0e3;   # MPa
    c44 = 75.0e3;    # MPa
    
    # Multiple slip system sets
    nsets = 3;
    
    # Set 1: {110}<111> slip
    plane1 = [1.0, 1.0, 0.0];
    dir1   = [1.0, 1.0, 1.0];
    
    # Set 2: {110}<111> slip (different variant)
    plane2 = [1.0, 0.0, 1.0];
    dir2   = [1.0, 1.0, 1.0];
    
    # Set 3: {110}<111> slip (different variant)
    plane3 = [0.0, 1.0, 1.0];
    dir3   = [1.0, 1.0, 1.0];
    
    # Crystal orientation (45-degree rotation about z-axis)
    # This will align the slip systems favorably for shear
    orient1_local  = [1.0, 0.0, 0.0];
    orient1_global = [0.707, 0.707, 0.0];  # cos(45°), sin(45°), 0
    orient2_local  = [0.0, 1.0, 0.0];
    orient2_global = [-0.707, 0.707, 0.0]; # -sin(45°), cos(45°), 0
    
    # Viscoplasticity parameters
    gamma0 = 0.001;   # Reference shear strain rate (1/s)
    m      = 30.0;    # Higher rate sensitivity
    
    # Hardening parameters (copper)
    tau0 = 60.0;      # Initial CRSS (MPa)
    h0   = 200.0;     # Initial hardening modulus (MPa)
    taus = 180.0;     # Saturation stress (MPa)
    a    = 2.0;       # Hardening exponent
    qab  = 1.4;       # Latent hardening ratio
    
    # Integration parameters
    theta   = 0.5;
    nlgeom  = false;
    maxiter = 25;
    tol     = 1.0e-6;
  };
};

solver =
{
  type = "NonlinearSolver";
  
  maxCycle = 100;
  dtime    = 0.05;  # Smaller time step for better convergence
  tol      = 1.0e-4;
};

outputModules = ["vtk", "GraphWriter", "OutputWriter"];

vtk =
{
  type     = "MeshWriter";
  prefix   = "oriented_shear";
  interval = 2;
};

GraphWriter =
{
  type     = "GraphWriter";
  onScreen = true;
  
  columns = ["time", "disp", "force", "gamma"];
  
  disp =
  {
    type = "state";
    node = 10;
    dof  = "u";
  };
  
  force =
  {
    type = "fint";
    node = 10;
    dof  = "u";
  };
  
  gamma =
  {
    type = "material";
    label = "GammaTotal";
    elem = 1;
    point = 1;
  };
};

OutputWriter =
{
  type     = "OutputWriter";
  onScreen = true;
};
