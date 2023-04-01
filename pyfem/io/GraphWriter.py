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
from pyfem.util.dataStructures import Properties
from pylab import plot, show, xlabel, ylabel, draw, ion, figure, gcf
from numpy import ndarray,zeros

class GraphWriter( BaseModule ):

  def __init__ ( self, props , globdat ):

    self.prefix    = globdat.prefix
    self.extension = ".out"
    self.onScreen  = False

    BaseModule.__init__( self , props )

    if not hasattr( self , "filename" ):
      self.filename  = self.prefix + self.extension
    
    self.columndata = []

    for i,col in enumerate ( self.columns ):

      if hasattr( self , col ):
        colProps = getattr( self , col )    
      else:
        colProps = Properties()
        
      if not hasattr( colProps , "type" ):
        colProps.type = col   
      
      if not hasattr( colProps , "factor" ):
        colProps.factor = 1.0
        
      if hasattr( colProps , "node" ):
        if type(colProps.node) == str:
          colProps.node = globdat.nodes.groups[colProps.node]

      self.columndata.append( colProps )

    if self.onScreen and hasattr( globdat , "onScreen" ):
      self.onScreen = False
    else:
      globdat.onScreen = True

      self.fig = gcf()
      self.fig.show()
      self.fig.canvas.draw()

    self.outfile = open( self.filename ,'w' )

    if self.onScreen:
      self.output = []

      xlabel(self.columns[0])
      ylabel(self.columns[1])
  
      ion()

    self.run( props , globdat ) 

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def run( self , props , globdat ):

    a = []

    for i,col in enumerate(self.columndata):
            
      if col.type in globdat.outputNames:
        data = globdat.getData( col.type , col.node )
        
      elif hasattr(globdat,col.type):
        b = getattr( globdat , col.type )
        if type(b) is ndarray:
          if type(col.node) is list:
            data = 0.0
            for nod in col.node:
              data += b[globdat.dofs.getForType(int(nod),col.dof)]
          else:
            data = b[globdat.dofs.getForType(col.node,col.dof)]
        else:
          data = b
      elif col.type in globdat.outputNames:
        data = globdat.getData( col.type , col.node )
      elif hasattr(globdat.solverStatus,col.type):
        data = getattr(globdat.solverStatus,col.type)
   
      data = data * col.factor

      a.append(data)
   
      self.outfile.write(str(data)+' ',)
      self.outfile.flush()

    self.outfile.write('\n')

    if self.onScreen: 
      self.output.append( a )

      plot( [x[0] for x in self.output], [x[1] for x in self.output], 'ro-' )
      self.fig.canvas.draw()
    
    if not globdat.active:
      self.outfile.close
