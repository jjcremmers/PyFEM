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

from pyfem.util.dataStructures import Properties
import re

def containsValue( db , val ):

  '''
  Addition of two numbers

  :param db: a
  :type db:  integer
  :param val: btted
  :type val:  integer
  :returns: new value
  :rtype:   integer
  '''

  keys = list(db.keys())
  
  for key in keys:

    if type(db[key]) == dict:
      if containsValue(db[key],val):
        return True
    else:
      if db[key] == val:
        return True
  return False
  
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def cleanVariable( a ):

  if a == 'true':
    return True
  elif a == 'false':
    return False
  else:
    try:
      return eval(a)
    except:
      return a
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def isNodeDof( nodeDof ):

  return type(nodeDof) == str and '[' in nodeDof

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def decodeNodeDof( nodeDof , nodes ):
            
  a       = nodeDof.split('[')
  dofType = a[0]  
  nodeID  = cleanVariable(a[1].split(']')[0])
  
  if type(nodeID) == str:
    return dofType,nodes.groups[nodeID]
  else:
    return dofType,[nodeID]
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def getType( a ):

  return type(cleanVariable(a))

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def storeValue( db , key , a ):

  if type(a) == list:
    tmp=[]
    for v in a:
      tmp.append(cleanVariable(v))
    db.store( key , tmp )
  else:
    db.store( key , cleanVariable(a) )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def readItem( l1 , db ):

  if '.' in l1[0]:
    l2 = l1[0].split('.',1)

    if l2[0] in db:
      if type(db[l2[0]]) == dict:
        child=db[l2[0]]
      else:
        child=Properties()
    else:
      child=Properties()

    l1[0]=l2[1]

    ln = readItem( l1 , child )

    db[l2[0]] = child

    return ln

  else:
    l2 = l1[1].split(';',1)

    if l2[0][0] == '[':
      l3 = l2[0][1:-1].split(',')
      storeValue( db , l1[0] , l3 )
    else:
      storeValue( db , l1[0] , l2[0] )
   
    return l2[1]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def readBlock( ln , db ):

  while True:

    if ln[0:7] == 'include':
      l1 = ln.split(';',1)
      deepFileParser( l1[0][8:-1] , db )
      ln = l1[1]
      continue

    l1 = ln.split('=',1)

    if len(l1) == 1:
      return ln

    if l1[0][0:2] == '};':
      return ln[2:]

    if l1[0][0:2] == '//':
      ln = l1[1].split(';',1)[1]
      continue

    #if l1[0][0:1] == '#':
    #  ln = l1[1].split(';',1)[1]
    #  continue

    if l1[1][0] == '{':
      child = Properties()
      ln = l1[1][1:]
     
      ln = readBlock( ln , child )

      db.store( l1[0] , child )

    else:     
      ln = readItem( l1 , db )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def fileParser( fileName ):

  db = Properties()

  f = open(fileName)
  
  f2 = ''
 
  for line in f:
    if not line.startswith('#'):
      f2 = f2+line
    
  ln = open(fileName).read().replace('\n','').replace('\t','').replace(' ','').replace('\r','')
  ln = f2.replace('\n','').replace('\t','').replace(' ','').replace('\r','')

  readBlock( ln , db )

  return db

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def deepFileParser( fileName , db ):

  ln = open(fileName).read().replace('\n','').replace('\t','').replace(' ','').replace('\r','')

  readBlock( ln , db )

  return db

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class nodeTable :

  def __init__( self , label , subLabel = "None" ):
  
    self.label    = label   
    self.subLabel = subLabel     
    self.data     = []

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def readNodeTable( fileName , label , nodes = None ):

  fin        = open( fileName , 'r' )

  startLabel = str('<'+label)
  endLabel   = str('</'+label)

  output = []

  for line in fin:

    if line.strip().startswith(startLabel) == True:

      nt = nodeTable( label )

      if 'name' in line:
        subLabel = line.split('=')[1].replace('\n','').replace('>','').replace(' ','').replace('\"','').replace('\'','')
        nt.subLabel = subLabel

      for line in fin:

        if line.strip().startswith(endLabel) == True:
          output.append(nt)
          break
        
        fullRel = line.strip().split(';')           

        if len(fullRel) == 2:
          splitRel = fullRel[0].split('=')

          if len(splitRel) == 2:
            lhs = splitRel[0]
            rhs = splitRel[1]

            if not isNodeDof(lhs):
              raise RuntimeError(str(lhs) + ' is not a NodeDof')
            
            dofType,nodeIDs = decodeNodeDof(lhs,nodes)
            
            if getType(rhs) is float or getType(rhs) is int:  
              for nodeID in nodeIDs:             
                nt.data.append([ dofType,int(nodeID),float(eval(rhs)) ])
            else:              
              rhs = rhs.replace(" ","").replace("+"," +").replace("-"," -")
              splitrhs = rhs.split(" ")
              rhs = 0.0
              for irhs in splitrhs:
                if irhs == "":
                  continue
                if '[' not in irhs:
                  for nodeID in nodeIDs:
                    nt.data.append([ dofType,int(nodeID),float(eval(irhs)) ])
                else:
                  eq_rhs = irhs.split("*")
                  factor = 1.0
                  
                  for ieq_rhs in eq_rhs:
                    if (getType(ieq_rhs) is float) or (getType(ieq_rhs) is int):
                      factor = cleanVariable(ieq_rhs)
                    else:
                      if isNodeDof(ieq_rhs):
                        if '-' in ieq_rhs:
                          factor = -1.0;
                        ieq_rhs = ieq_rhs.replace("-","").replace("+","")
                        slaveDofType,slaveNodeID = decodeNodeDof( ieq_rhs , nodes)
                        
                  for nodeID in nodeIDs:
                    dt = [ dofType,int(nodeID),rhs,slaveDofType,slaveNodeID,factor ]
                    nt.data.append(dt)
                                                    
  return output

