############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  The code is written by J.J.C. Remmers, C.V. Verhoosel and R. de Borst.  #
#                                                                          #
#  The latest stable version can be downloaded from the web-site:          #
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
from pyfem.util.dataStructures import Properties

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

def getType( a ):

  if a == 'true':
    return True
  elif a == 'false':
    return False
  else:
    try:
      return eval(a)
    except:
      return a

def storeValue( db , key , a ):

  if type(a) == list:
    tmp=[]
    for v in a:
      tmp.append(getType(v))
    db.store( key , tmp )
  else:
    db.store( key , getType(a) )

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

def deepFileParser( fileName , db ):

  ln = open(fileName).read().replace('\n','').replace('\t','').replace(' ','').replace('\r','')

  readBlock( ln , db )

  return db

