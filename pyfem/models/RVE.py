# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

import sys
from venv import logger
from pyfem.models.BaseModel import BaseModel
import numpy as np


class RVE( BaseModel ):
    """Representative Volume Element (RVE) model for periodic boundary conditions.

    This model implements periodic boundary conditions on a rectangular RVE by
    constraining nodes on opposite boundaries (Left-Right, Top-Bottom) to maintain
    periodic displacement fields under prescribed macroscopic strain.

    Attributes
    ----------
    dx : float
        Width of the RVE (x-direction offset).
    dy : float
        Height of the RVE (y-direction offset).
    crdsLeft, crdsRight : numpy.ndarray
        Coordinates of left and right boundary nodes, sorted by y-coordinate.
    crdsTop, crdsBottom : numpy.ndarray
        Coordinates of top and bottom boundary nodes, sorted by x-coordinate.
    nodesLeft, nodesRight : numpy.ndarray
        Node IDs of left and right boundary nodes, sorted by y-coordinate.
    nodesTop, nodesBottom : numpy.ndarray
        Node IDs of top and bottom boundary nodes, sorted by x-coordinate.
    """

    def __init__ ( self, props , globdat ):
        """Initialize the RVE model and set up periodic boundary conditions.

        Parameters
        ----------
        props : Properties
            Model properties and configuration.
        globdat : GlobalData
            Global data structure containing nodes, elements, and DOFs.
        """
   
        BaseModel.__init__( self , props , globdat )

        if not hasattr( self , "boundaryType" ):
            self.boundaryType = "Periodic"

        if self.boundaryType not in [ "Periodic" , "Prescribed"]:
            sys.exit("Error: RVE boundaryType must be 'Periodic' or 'Prescribed'")  
                       
        self.getBoundaries( props , globdat )

        globdat.dofs.createConstrainer()   
    
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    def getTangentStiffness(self, props, globdat, mbuilder):
        """
        Apply periodic boundary constraints based on prescribed macroscopic strain.

        This method applies displacement constraints to corner nodes and periodic
        constraints to boundary nodes to enforce a prescribed macroscopic strain state.
        The strain vector follows Voigt notation: [ε_xx, ε_yy, γ_xy].

        Parameters
        ----------
        props : Properties
            Model properties (currently unused).
        globdat : GlobalData
            Global data structure containing nodes, elements, and DOFs.
        mbuilder : MatrixBuilder
            MatrixBuilder instance (not used in this model, included for interface compatibility).

        Notes
        -----
        The method:
        1. Defines macroscopic strain components (currently hardcoded).
        2. Calculates corner node displacements from strain.
        3. Constrains the bottom-left corner node (node 1) to zero displacement.
        4. Applies calculated displacements to the three other corner nodes.
        5. Applies periodic constraints to all interior boundary nodes.
        """
        strain = np.zeros(3)
        strain[0] = 0.01
        strain[1] = 0.02
        strain[2] = 0.02

        if self.boundaryType == "Prescribed":
            self.applyPrescribedBC(strain, props, globdat)
        elif self.boundaryType == "Periodic":
            self.applyPeriodicBC(strain, props, globdat)
        else:
            sys.exit("Error: boundaryType must be 'Periodic' or 'Prescribed'")

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getBoundaries( self , props , globdat ):
        """Extract and validate boundary nodes for periodic RVE constraints.

        This method retrieves boundary node coordinates and IDs from node groups
        (Left, Right, Top, Bottom), sorts them for proper pairing, and validates
        the RVE geometry to ensure:
        - Matching numbers of nodes on opposite boundaries.
        - Proper alignment of node coordinates.
        - Consistent spacing between paired nodes.

        Parameters
        ----------
        props : Properties
            Model properties (currently unused).
        globdat : GlobalData
            Global data structure containing nodes, elements, and DOFs.

        Raises
        ------
        Warning
            If boundary nodes don't match in number, alignment, or spacing.

        Notes
        -----
        Sets the following instance attributes:
        - crdsLeft, crdsRight: Sorted coordinates of left/right boundaries.
        - crdsTop, crdsBottom: Sorted coordinates of top/bottom boundaries.
        - nodesLeft, nodesRight: Sorted node IDs of left/right boundaries.
        - nodesTop, nodesBottom: Sorted node IDs of top/bottom boundaries.
        - dx: Width of the RVE.
        - dy: Height of the RVE.
        """

        crdL = globdat.nodes.getNodeCoords( 'Left' )
        crdR = globdat.nodes.getNodeCoords( 'Right' )

        # Sort both left and right coordinates by y-coordinate (column 1)
        sortIdxLeft = np.argsort(crdL[:, 1])
        sortIdxRight = np.argsort(crdR[:, 1])
        
        self.crdsLeft = crdL[sortIdxLeft]
        self.crdsRight = crdR[sortIdxRight]

        nodL = np.array(globdat.nodes.getNodeIDs( 'Left' ))
        nodR = np.array(globdat.nodes.getNodeIDs( 'Right' ))
        
        self.nodesLeft = nodL[sortIdxLeft]
        self.nodesRight = nodR[sortIdxRight]
        # Verify that y-coordinates match (within tolerance)

        x_offsets = self.crdsRight[:, 0] - self.crdsLeft[:, 0]
        
        tol = 1e-6 * np.max( np.abs(x_offsets) )

        if len(self.crdsLeft) != len(self.crdsRight):
            logger.warning(f"Warning: Different number of nodes - Left: {len(self.crdsLeft)}, Right: {len(self.crdsRight)}")
        elif not np.allclose(self.crdsLeft[:, 1], self.crdsRight[:, 1], atol=tol):
            logger.warning("Warning: Y-coordinates do not match between Left and Right")
            logger.warning("Left y-coords:", self.crdsLeft[:, 1])
            logger.warning("Right y-coords:", self.crdsRight[:, 1])

        # Calculate x-offset for all pairs
        self.dx = x_offsets[0]
                
        if not np.allclose(x_offsets, self.dx, atol=tol):
            logger.warning("Warning: X-offsets vary between pairs")
            for i, offset in enumerate(x_offsets):
                print(f"  Pair {i}: offset = {offset} (diff = {offset - self.dx})")

        # Process Top and Bottom node groups
        crdT = globdat.nodes.getNodeCoords( 'Top' )
        crdB = globdat.nodes.getNodeCoords( 'Bottom' )

        # Sort both top and bottom coordinates by x-coordinate (column 0)
        sortIdxTop = np.argsort(crdT[:, 0])
        sortIdxBottom = np.argsort(crdB[:, 0])
        
        self.crdsTop = crdT[sortIdxTop]
        self.crdsBottom = crdB[sortIdxBottom]

        nodT = np.array(globdat.nodes.getNodeIDs( 'Top' ))
        nodB = np.array(globdat.nodes.getNodeIDs( 'Bottom' ))
        
        self.nodesTop = nodT[sortIdxTop]
        self.nodesBottom = nodB[sortIdxBottom]    

        # Verify that x-coordinates match (within tolerance)
        if len(crdT) != len(crdB):
            logger.warning(f"Warning: Different number of nodes - Top: {len(crdT)}, Bottom: {len(crdB)}")
        elif not np.allclose(self.crdsTop[:, 0], self.crdsBottom[:, 0], atol=tol):
            logger.warning("Warning: X-coordinates do not match between Top and Bottom")
            logger.warning("Top x-coords:", self.crdsTop[:, 0])
            logger.warning("Bottom x-coords:", self.crdsBottom[:, 0])
  
        # Calculate y-offset for all pairs
        y_offsets = self.crdsTop[:, 1] - self.crdsBottom[:, 1]
        self.dy   = y_offsets[0]
                
        if not np.allclose(y_offsets, self.dy, atol=tol):
            print("Warning: Y-offsets vary between pairs")
            for i, offset in enumerate(y_offsets):
                print(f"  Pair {i}: offset = {offset} (diff = {offset - self.dy})")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

    def applyPeriodicBC(self, strain: np.ndarray, props, globdat) -> None:
        
        """Apply periodic boundary constraints for prescribed macroscopic strain.

        Implements periodic boundary conditions by:
        1. Fixing the bottom-left corner node (node 1) to zero displacement.
        2. Prescribing displacements at the three other corner nodes based on strain.
        3. Coupling interior boundary nodes on opposite sides with periodic constraints.

        Parameters
        ----------
        strain : numpy.ndarray
            Macroscopic strain tensor in Voigt notation [ε_xx, ε_yy, γ_xy].
            Components define the prescribed deformation of the RVE.
        props : Properties
            Model properties (currently unused).
        globdat : GlobalData
            Global data structure containing nodes, elements, and DOFs.

        Notes
        -----
        Corner node displacements are computed from:
        - u12 (bottom-right): u_x = dx * ε_xx, u_y = 0.5 * dx * γ_xy
        - u14 (top-left): u_x = 0.5 * dy * γ_xy, u_y = dy * ε_yy
        - u34 (top-right): u = u12 + u14

        Interior nodes are coupled via linear constraints:
        - Right nodes: u_right = u12 + u_left
        - Top nodes: u_top = u14 + u_bottom
        """

        u12 = np.zeros(2)
        u14 = np.zeros(2)

        u12[0] =      self.dx * strain[0]
        u12[1] = 0.5* self.dx * strain[2] 

        u14[0] = 0.5* self.dy * strain[2]  
        u14[1] =      self.dy * strain[1]

         # Constrain node 1 (left bottom corner)
        
        dofIDs = globdat.dofs.get( [self.nodesLeft[0]] )
        for i in dofIDs:
            globdat.dofs.cons.addConstraint(i,0.0,"main")
            
        # Constrain node 2 bottom right corner    

        dofIDs = globdat.dofs.get( [self.nodesRight[0]] )
        for i in range(2):
            globdat.dofs.cons.addConstraint(dofIDs[i],u12[i],"main")

        # Constrain node 4 top left corner
        dofIDs = globdat.dofs.get( [self.nodesTop[0]] )
        for i in range(2):
            globdat.dofs.cons.addConstraint(dofIDs[i],u14[i],"main")  

        # Constrain node 3 top right corner
        dofIDs = globdat.dofs.get( [self.nodesTop[-1]] )
        for i in range(2):
            globdat.dofs.cons.addConstraint(dofIDs[i],u12[i]+u14[i],"main")              

        # Connect top and bottom nodes

        for nodeBottom, nodeTop in zip( self.nodesBottom[1:-1] , self.nodesTop[1:-1] ):
            dofIDsBottom  = globdat.dofs.get( nodeBottom )
            dofIDsTop = globdat.dofs.get( nodeTop )

            for i in range(2):
                globdat.dofs.cons.addConstraint( dofIDsTop[i] , [u14[i] , dofIDsBottom[i] , 1.0 ] , "main" )

        # Connect left and right nodes

        for nodeLeft, nodeRight in zip( self.nodesLeft[1:-1] , self.nodesRight[1:-1] ):
            dofIDsLeft  = globdat.dofs.get( nodeLeft )
            dofIDsRight = globdat.dofs.get( nodeRight )

            for i in range(2):
                globdat.dofs.cons.addConstraint( dofIDsRight[i] , [u12[i] , dofIDsLeft[i] , 1.0 ] , "main" )

        globdat.dofs.cons.flush()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

    def applyPrescribedBC(self, strain: np.ndarray, props, globdat) -> None:

        """Apply prescribed displacement boundary conditions for uniform strain.

        Prescribes displacements at all boundary nodes based on a linear displacement
        field corresponding to the given macroscopic strain. This approach is simpler
        than periodic BCs but does not allow periodic fluctuations in the displacement
        field.

        Parameters
        ----------
        strain : numpy.ndarray
            Macroscopic strain tensor in Voigt notation [ε_xx, ε_yy, γ_xy].
        props : Properties
            Model properties (currently unused).
        globdat : GlobalData
            Global data structure containing nodes, elements, and DOFs.

        Notes
        -----
        For each boundary node at position (x, y), the prescribed displacements are:
        - u_x = ε_xx * x + 0.5 * γ_xy * y
        - u_y = ε_yy * y + 0.5 * γ_xy * x

        All boundary nodes (Left, Right, Top, Bottom) are collected and constrained,
        with corner nodes appearing only once due to the use of np.unique().
        """

        allBoundaryNodes = np.unique(np.concatenate([
                                        self.nodesLeft,
                                        self.nodesRight,
                                        self.nodesTop,
                                        self.nodesBottom]))
        
        crds = globdat.nodes.getNodeCoords( allBoundaryNodes.tolist() )
        
        for nodeID,crd in zip(allBoundaryNodes,crds):
            dofIDs = globdat.dofs.get( [nodeID] )

            uX = strain[0]*crd[0] + 0.5*strain[2]*crd[1]
            uY = strain[1]*crd[1] + 0.5*strain[2]*crd[0]

            globdat.dofs.cons.addConstraint( dofIDs[0] , uX , "main" )
            globdat.dofs.cons.addConstraint( dofIDs[1] , uY , "main" )

        globdat.dofs.cons.flush()                        
  
     
