############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  The code is written by J.J.C. Remmers, C.V. Verhoosel and R. de Borst.  #
#                                                                          #
#  The latest stabke version can be downloaded from the web-site:          #
#     http://www.wiley.com/go/deborst                                      #
#                                                                          #
#  A github repository, with the most up to date version of the code,      #
#  can be found here:                                                      #
#     https://github.com/jjcremmers/PyFEM                                  #
#                                                                          #
#  The code is open source and intended for educational and scientific     #
#  purposes only. If you use PyFEM in your research, the developers would  #
#  be grateful if you could cite the book.                                 #  
#                                                                          #
#  Disclaimer:                                                             #
#  The authors reserve all rights but do not guarantee that the code is    #
#  free from errors. Furthermore, the authors shall not be liable in any   #
#  event caused by the use of the program.                                 #
############################################################################
from numpy import array
from pyfem.util.itemList import itemList
import re,sys

class NodeSet( itemList ):

  def __init__( self ):
    self.rank = -1

  def getNodeCoords( self, nodeIDs ):
    return array( self.get( nodeIDs ) )
    
  def readFromFile( self, fname ):
    
    print("  Reading nodes ................")

    fin = open( fname )
  
    while True:
    
      line = fin.readline()  
  
      if line.startswith('<Nodes>') == True:
      
        while True:
          line = fin.readline()  

          if line.startswith('</Nodes>') == True:
            return
  
          line = re.sub('\s{2,}',' ',line)
          a = line.split(';')
     
          for a in a[:-1]:
            b = a.strip().split(' ')
            
            if b[0].startswith("//") or b[0].startswith("#"):
              break
            if len(b) > 1 and type(eval(b[0])) == int:
              if self.rank == -1:
                self.rank = len(b)-1

              self.add( eval(b[0]), [eval(crd) for crd in b[1:]] )
            

      elif line.startswith('gmsh') == True:
        ln = line.replace('\n','').replace('\t','').replace(' ','').replace('\r','').replace(';','')
        ln = ln.split('=',1)
        self.readGmshFile( ln[1][1:-1] )
        return

  def readGmshFile( self, fname ):
    
    fin = open( fname )

    while True:
    
      line = fin.readline()  
  
      if line.startswith('$MeshFormat') == True:
      
        while True:
          line = fin.readline()  
          line = re.sub('\s{2,}',' ',line)

          a = line.split(';')
          b = a[0].strip().split(' ')
      
          if eval(b[0]) < 2.0:
            print("error")
            sys.exit()
          
          break

      if line.startswith('$Nodes'):
  
        nNodes = eval(fin.readline())

        for i in range(nNodes):
          line = fin.readline()
          line = re.sub('\s{2,}',' ',line)
          b    = line.strip().split(' ')

          if len(b) > 1 and type(eval(b[0])) == int:
            self.add( eval(b[0]), [eval(crd) for crd in b[1:3]] )
      
      if line.startswith('$EndNodes'):
        return 

#------
#
#-------

  def __repr__( self ):
    return "Nodeset contains %i nodes.\n" % len(self)
