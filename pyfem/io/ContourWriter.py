################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2025. The code is written in 2011-2012 by                #
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

from pyfem.util.BaseModule import BaseModule

class ContourWriter( BaseModule ):
 
  def __init__( self , props , globdat ):

    self.prefix       = globdat.prefix
    self.interval     = 1

    BaseModule.__init__( self , props )
    
    self.k            = 0
    self.columndata   = []

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
	
  def run( self , props , globdat ):

    if not globdat.solverStatus.cycle%self.interval == 0:
      return
     
    self.writeHeader()
 
    crd = globdat.nodes.getNodeCoords(self.nodes[0])
    outfile = open( self.prefix + '-contour-' + str(self.k) + '.out' ,'w' )
        
    outfile.write(f'#Node  {"x-coor":<10} {"y-coor":<10}')
  
    if len(crd) == 3:
      outfile.write(f'{"z-coor":<10} ')

    for dofType in globdat.dofs.dofTypes:
      outfile.write(f'{dofType:<10} ')

    for name in globdat.outputNames:
      outfile.write(f'{name:<10} ')

    outfile.write('\n')

    for iNod in self.nodes:
      crd = globdat.nodes.getNodeCoords(iNod)
      outfile.write(f'{iNod:4d} {crd[0]:10.3e} {crd[1]:10.3e}')

      if len(crd) == 3:
        outfile.write(f' {crd[2]:10.3e}')

      for dofType in globdat.dofs.dofTypes:
        outfile.write(f' {globdat.state[globdat.dofs.getForType(iNod, dofType)]:10.3e}')

      for name in globdat.outputNames:
        stress = globdat.getData(name, list(range(len(globdat.nodes))))
        outfile.write(f' {stress[iNod]:10.3e}')

      outfile.write('\n')
                
    outfile.close()
  
    self.k = self.k+1
