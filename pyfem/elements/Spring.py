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

from .Element import Element
from pyfem.util.transformations import toElementCoordinates, toGlobalCoordinates

from numpy import zeros, eye, array

class Spring ( Element ):

  #Number of dofs per element
  dofTypes = ['u','v']

  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )


  def __type__ ( self ):
    return name


  def getTangentStiffness ( self, elemdat ):
    
    #Compute the current state vector

    a  = toElementCoordinates( elemdat.state  , elemdat.coords )
    Da = toElementCoordinates( elemdat.Dstate , elemdat.coords )

    #Compute the elongation of the spring
    elong = a[2]-a[0] 

    #Compute the force in the spring
    Fs = elong * elemdat.props.k

    #Compute the element internal force vector in the element coordinate system
    elFint = array([-Fs,0.,Fs,0])

    #Determine the element tangent stiffness in the element coordinate system
    elKbar = zeros( (4,4) )

    elKbar[:2,:2] =  elemdat.props.k*eye(2)
    elKbar[:2,2:] = -elemdat.props.k*eye(2)

    elKbar[2:,:2] = elKbar[:2,2:]
    elKbar[2:,2:] = elKbar[:2,:2]

    #Rotate element tangent stiffness to the global coordinate system
    elemdat.stiff = toGlobalCoordinates( elKbar, elemdat.coords )
    elemdat.fint = toGlobalCoordinates( elFint, elemdat.coords )
  
#------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):

    #Compute the current state vector

    a  = toElementCoordinates( elemdat.state  , elemdat.coords )
    Da = toElementCoordinates( elemdat.Dstate , elemdat.coords )

    #Compute the elongation of the spring
    elong = a[2]-a[0] 

    #Compute the force in the spring
    Fs = elong * elemdat.props.k

    #Compute the element internal force vector in the element coordinate system
    elFint = array([-Fs,0.,Fs,0])

    #Rotate element fint to the global coordinate system
    elemdat.fint = toGlobalCoordinates( elFint, elemdat.coords )
