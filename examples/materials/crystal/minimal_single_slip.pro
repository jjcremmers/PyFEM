############################################################################
#  Crystal Plasticity Example - Minimal Single-Slip Tension               #
#                                                                          #
#  Description: Smallest crystal plasticity example in the repository.     #
#               Uses a single quad element and one slip-system set.        #
#                                                                          #
#  Usage:       pyfem minimal_single_slip.pro                              #
############################################################################

input = "minimal_single_slip.dat";

ContElem =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "Crystal";

    # Isotropic elastic fallback to keep the example compact.
    E  = 70000.0;
    nu = 0.33;

    # One slip-system set.
    nsets = 1;
    plane1 = [1.0, 1.0, 0.0];
    dir1   = [1.0, 1.0, 1.0];

    # Identity crystal orientation.
    orient1_local  = [1.0, 0.0, 0.0];
    orient1_global = [1.0, 0.0, 0.0];
    orient2_local  = [0.0, 1.0, 0.0];
    orient2_global = [0.0, 1.0, 0.0];

    # Rate-dependent slip law.
    gamma0 = 0.001;
    m      = 20.0;

    # Hardening parameters.
    tau0 = 25.0;
    h0   = 150.0;
    taus = 120.0;
    a    = 2.0;
    qab  = 1.0;

    theta   = 0.5;
    maxiter = 20;
    tol     = 1.0e-6;
  };
};

solver =
{
  type = "NonlinearSolver";
  maxCycle = 20;
  dtime    = 0.1;
};

outputModules = ["OutputWriter"];

OutputWriter =
{
  type     = "OutputWriter";
  onScreen = true;
};