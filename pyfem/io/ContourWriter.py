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

from pyfem.util.BaseModule import BaseModule
from pyfem.util.logger   import getLogger

logger = getLogger()

class ContourWriter( BaseModule ):
 
  def __init__( self , props , globdat ):

    self.prefix       = globdat.prefix
    self.interval     = 1
    self.k            = 0

#    self.stresslabels = [ "sxx" , "syy" , "sxy" ]

    BaseModule.__init__( self , props )

    self.columndata = []

#    for i,col in enumerate ( self.columns ):
#      self.columndata.append( colProps )
	
  def run( self , props , globdat ):

    if not globdat.solverStatus.cycle%self.interval == 0:
      return
     
    logger.info("Writing contour file ......\n")
 
    crd = globdat.nodes.getNodeCoords(self.nodes[0])
    outfile = open( self.prefix + '-contour-' + str(self.k) + '.out' ,'w' )
        
    #tractions = globdat.getData( "tractions" , range(len(globdat.nodes)) )

    outfile.write( '#Node  %-10s %-10s' % ('x-coor','y-coor') )
  
    if len(crd) == 3:
      outfile.write( '%-10s ' % 'z-coor' )

    for dofType in globdat.dofs.dofTypes:
      outfile.write( '%-10s ' % dofType )

    for name in globdat.outputNames:
      outfile.write('%-10s ' % name )

    outfile.write('\n')

    for iNod in self.nodes:
      crd = globdat.nodes.getNodeCoords(iNod)
      outfile.write('%4i %10.3e %10.3e' % (iNod,crd[0],crd[1]))
      
      if len(crd) == 3:
        outfile.write(' %10.3e' % crd[2] )
       
      for dofType in globdat.dofs.dofTypes:
        outfile.write(' %10.3e' % (globdat.state[globdat.dofs.getForType(iNod,dofType)]))
      
      for name in globdat.outputNames:
        stress = globdat.getData( name , list(range(len(globdat.nodes))) )    
        outfile.write(' %10.3e' % stress[iNod] )

      outfile.write('\n')
                
    outfile.close()
  
    self.k = self.k+1
