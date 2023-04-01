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

from numpy import zeros, ones, ix_ , append, repeat, array
from scipy.sparse import coo_matrix
from pyfem.util.dataStructures import Properties
from pyfem.util.dataStructures import elementData


#######################################
# General array assembly routine for: # 
# * assembleInternalForce             #
# * assembleTangentStiffness          #
#######################################

def assembleArray ( props, globdat, rank, action ):

  #Initialize the global array A with rank 2

  B = zeros( len(globdat.dofs) * ones(1,dtype=int) )
  cc = 0.0

  val   = array([],dtype=float)
  row   = array([],dtype=int)
  col   = array([],dtype=int)

  nDof  = len(globdat.dofs)

  if action != 'commit':
    globdat.resetNodalOutput()

  #Loop over the element groups
  for elementGroup in globdat.elements.iterGroupNames():

    #Get the properties corresponding to the elementGroup
    el_props = getattr( props, elementGroup )

    #Loop over the elements in the elementGroup
    for iElm,element in enumerate(globdat.elements.iterElementGroup( elementGroup )):

      #Get the element nodes
      el_nodes = element.getNodes()

      #Get the element coordinates
      el_coords = globdat.nodes.getNodeCoords( el_nodes )

      #Get the element degrees of freedom
      el_dofs = globdat.dofs.getForTypes( el_nodes , element.dofTypes )
      
      #Get the element state
      el_a  = globdat.state [el_dofs]
      el_Da = globdat.Dstate[el_dofs]

      #Create the an element state to pass through to the element
      #el_state = Properties( { 'state' : el_a, 'Dstate' : el_Da } )
      elemdat = elementData( el_a , el_Da )

      elemdat.coords   = el_coords
      elemdat.nodes    = el_nodes
      elemdat.props    = el_props
      elemdat.iElm     = iElm 

      element.globdat  = globdat
      
      if hasattr( element , "matProps" ):
        elemdat.matprops = element.matProps

      if hasattr( element , "mat" ):
        element.mat.reset()

      #Get the element contribution by calling the specified action
      if hasattr( element , action ):
        getattr( element, action )( elemdat )

      #for label in elemdat.outlabel:	
      #  element.appendNodalOutput( label , globdat , elemdat.outdata )

      #Assemble in the global array
      if rank == 1:
        B[el_dofs] += elemdat.fint
        cc         += elemdat.diss
      elif rank == 2 and action == "getTangentStiffness":  

        row = append(row,repeat(el_dofs,len(el_dofs)))

        for i in range(len(el_dofs)):
          col=append(col,el_dofs)        

        val = append(val,elemdat.stiff.reshape(len(el_dofs)*len(el_dofs)))

        B[el_dofs] += elemdat.fint
      elif rank == 2 and action == "getMassMatrix": 

        row = append(row,repeat(el_dofs,len(el_dofs)))

        for i in range(len(el_dofs)):
          col=append(col,el_dofs)        

        val = append(val,elemdat.mass.reshape(len(el_dofs)*len(el_dofs))) 

        B[el_dofs] += elemdat.lumped
  #    else:
  #      raise NotImplementedError('assemleArray is only implemented for vectors and matrices.')

  if rank == 1:
    return B,cc
  elif rank == 2:

    '''if globdat.contact.flag:
      row , val , col = globdat.contact.checkContact( row , val , col , B , globdat )      '''

    return coo_matrix((val,(row,col)), shape=(nDof,nDof)),B


##########################################
# Internal force vector assembly routine # 
##########################################

def assembleInternalForce ( props, globdat ):
  fint = assembleArray( props, globdat, rank = 1, action = 'getInternalForce' )
  return fint[0]

##########################################
# External force vector assembly routine # 
##########################################

def assembleExternalForce ( props, globdat ):
  fext = assembleArray( props, globdat, rank = 1, action = 'getExternalForce' )   

  return fext[0] + globdat.fhat * globdat.solverStatus.lam
  
def assembleDissipation ( props, globdat ):
  return assembleArray( props, globdat, rank = 1, action = 'getDissipation' )   
 
#############################################
# Tangent stiffness matrix assembly routine # 
#############################################

def assembleTangentStiffness ( props, globdat ):
  return assembleArray( props, globdat, rank = 2, action = 'getTangentStiffness' )

#############################################
# Mass matrix assembly routine              # 
#############################################

def assembleMassMatrix ( props, globdat ):
  return assembleArray( props, globdat, rank = 2, action = 'getMassMatrix' )

def commit ( props, globdat ):
  return assembleArray( props, globdat, rank = 0, action = 'commit' )

#
#
#

def getAllConstraints ( props , globdat ):

  #Loop over the element groups
  for elementGroup in globdat.elements.iterGroupNames():

    #Get the properties corresponding to the elementGroup
    el_props = getattr( props, elementGroup )

    #Loop over the elements in the elementGroup
    for element in globdat.elements.iterElementGroup( elementGroup ):

      #Get the element nodes
      el_nodes = element.getNodes()

      elemdat.nodes    = el_nodes
      elemdat.props    = el_props
      
      #Get the element contribution by calling the specified action
      getattr( element, 'getConstraints', None )( elemdat )

