################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2023. The code is written in 2011-2012 by                #
#  Joris J.C. Remmers, Clemens V. Verhoosel and Rene de Borst and since        #
#  then augmented and maintained by Joris J.C. Remmers.                        #
#  All rights reserved.                                                        #
#                                                                              #
#  A github repository, with the most up to date version of the code,          #
#  can be found here:                                                          #
#     https://github.com/jjcremmers/PyFEM/                                     #
#     https://pyfem.readthedocs.io/                                            #	
#                                                                              #
#  The original code can be downloaded from the web-site:                      #
#     http://www.wiley.com/go/deborst                                          #
#                                                                              #
#  The code is open source and intended for educational and scientific         #
#  purposes only. If you use PyFEM in your research, the developers would      #
#  be grateful if you could cite the book.                                     #    
#                                                                              #
#  Disclaimer:                                                                 #
#  The authors reserve all rights but do not guarantee that the code is        #
#  free from errors. Furthermore, the authors shall not be liable in any       #
#  event caused by the use of the program.                                     #
################################################################################

from numpy import zeros, eye, cross
from scipy.linalg import norm

from pyfem.util.utilFunctions import tensor2array, array2tensor
from pyfem.util.utilFunctions import sintot, costot2, sintot3

#===============================================================================
#
#===============================================================================


class BeamRotation:

    def __init__( self , nItm ):
    
        self.lamNew = zeros(shape=(3, 3, nItm))
        self.lamOld = zeros(shape=(3, 3, nItm))
        self.omeNew = zeros(shape=(3, nItm))
        self.omeOld = zeros(shape=(3, nItm))

        for iItm in range(nItm):
            self.lamOld[0, 0, iItm] = 1.0
            self.lamOld[1, 1, iItm] = 1.0
            self.lamOld[2, 2, iItm] = 1.0

    def __del__(self):
        pass

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def update( self , iInt , kin ):
    
        tnorm = norm(kin.theta)

        theta_tensor = array2tensor( kin.theta )
        theta_tensor2 = theta_tensor @ theta_tensor

        texp = sintot( tnorm ) * theta_tensor + costot2( tnorm ) * theta_tensor2 + eye( 3 )

        kin.lambda_ = texp @ self.lamOld[:, :, iInt]

        g = sum( kin.dtheta * kin.theta )
        gv = cross( kin.dtheta, kin.theta )

        ka = (
            sintot( tnorm ) * kin.dtheta
            + costot2( tnorm ) * gv
            + sintot3( tnorm ) * g * kin.theta
        )

        kappa_tensor = array2tensor( ka )
        omega_tensor = array2tensor( self.omeOld[:, iInt] )

        kappa_tensor = ( kappa_tensor + omega_tensor ) @ texp.T
        kappa_tensor = texp @ kappa_tensor

        omearr = tensor2array( kappa_tensor )
        kin.kappa = kin.lambda_.T @ omearr

        self.lamNew[:, :, iInt] = kin.lambda_
        self.omeNew[:, iInt]    = omearr

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def commit(self):
    
        self.lamOld = self.lamNew.copy()
        self.omeOld = self.omeNew.copy()

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------


def glob2loc(loc, trans, glob):

    for i in range(4):
        loc[3*i:3*(i+1)] = trans @ glob[3*i:3*(i+1)]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def loc2glob(glob, trans, loc):

    if loc.ndim == 1:
        for i in range(4):
            glob[3*i:3*(i+1)] = trans.T @ loc[3*i:3*(i+1)]
    elif loc.ndim == 2:
        tt = zeros(shape=(12, 12))
    
        for i in range(4):
            tt[3*i:3*(i+1), 3*i:3*(i+1)] = trans
            
        glob[:] = tt.T @ (loc @ tt)
