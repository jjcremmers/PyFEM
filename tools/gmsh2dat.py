############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  Copyright (C) 2011-2024. The code is written in 2011-2012 by            #
#  Joris J.C. Remmers, Clemens V. Verhoosel and Rene de Borst and since    #
#  then augmented and  maintained by Joris J.C. Remmers.                   #
#  All rights reserved.                                                    #
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

import sys,os,meshio,getopt

def printGroup( outFile , group , length = 10 ):

  print('{',file=outFile,end=(" "))
  counter = 0
  for itemID in group: 
    if counter < length:  
      print(itemID,file=outFile,end=(" "))
      counter += 1
    else:
      print(itemID,file=outFile,end=("\n"))
      counter = 0 
      
  print('}',file=outFile)
  
def printHelp():

  print(" gsmh2dawn\n")
  print(" Usage\n")
  print("   gmsh2dawn filename.msh\n")
  print("   -o output file\n")
  
#-------------------------------------------------------------------------------
#  checkRank
#-------------------------------------------------------------------------------

def checkRank( mesh ):

  rank = 2
  
  for key in mesh.cell_sets_dict:
    for typ in mesh.cell_sets_dict[key]:
      if typ[:4] in ["pris","pyra","hexa","wedg","tetr"]:
        rank = 3  
        
  return rank
  
#-------------------------------------------------------------------------------
#  printNodes
#-------------------------------------------------------------------------------

def printNodes( of , mesh ):

  rank = checkRank( mesh )
  
  print("<Nodes>",file=of) 
  
  for nodeID,p in enumerate(mesh.points):
    if rank == 2:
      print(nodeID,p[0],p[1],' ;',file=of)
    else:
      print(nodeID,p[0],p[1],p[2],' ;',file=of)
    
  print("\n</Nodes>\n",file=of)
  
  print("  Number of nodes ............ %6d" %len(mesh.points) )
  
#-------------------------------------------------------------------------------
#  printElements
#-------------------------------------------------------------------------------

def printElements( of , mesh ):

  elemID     = 0

  print("<Elements>",file=of)
  
  for key in mesh.cell_sets_dict:
    if key == "gmsh:bounding_entities":
      pass
    else: 
      for typ in mesh.cell_sets_dict[key]:
        for idx in mesh.cell_sets_dict[key][typ]:         
          iNodes = mesh.cells_dict[typ][idx]
          print("  %i '%s'" % (elemID,key),end=' ',file=of)
          for nodeID in iNodes:
            print(nodeID,end=' ',file=of)
          print(";",file=of)

          elemID += 1

  print("</Elements>\n",file=of) 
  
  print("  Number of elements ......... %6d" %elemID )

#-------------------------------------------------------------------------------
#  printGroups
#-------------------------------------------------------------------------------

def printGroups( of , mesh ):

  elemCount  = []
  nodeCount  = []
  groupNames = []
  
  for key in mesh.cell_sets_dict:
    iElm = 0
    if key == "gmsh:bounding_entities":
      pass
    else:
      grp = []
      for typ in mesh.cell_sets_dict[key]:
        for idx in mesh.cell_sets_dict[key][typ]:         
          iNodes = mesh.cells_dict[typ][idx]
          for nodeID in iNodes:
            grp.append(nodeID)
          iElm += 1
          
      print("<NodeGroup name='%s'>" %(key),file=of)
                
      printGroup( of , set(list(grp)) )
      
      print("</NodeGroup>\n",file=of)  
      
      elemCount .append(iElm)
      nodeCount .append(len(set(list(grp))))
      groupNames.append(key)
      
  print("  Number of  groups .......... %6d" %len(groupNames) )
  print("  -----------------------------------")
  print("    name               #nodes  #elems")
  print("    ---------------------------------")
  for name,nCount,eCount in zip(groupNames,nodeCount,elemCount):
    print("    %-16s    %6d %6d " %(name[:16],nCount,eCount))      
   
#-------------------------------------------------------------------------------
#  printRest
#-------------------------------------------------------------------------------

def printRest( of , mesh ):

  print("<NodeConstraints>\n",file=of)  
  print("</NodeConstraints>\n",file=of)  

  print("<ExternalForces>\n",file=of)  
  print("</ExternalForces>\n",file=of)        

#===============================================================================
#  Main file
#===============================================================================


slist = 'i:o:h'
llist = ['input=','output=','help','version']
  
argv = sys.argv

options, remainder = getopt.getopt( argv[1:] , slist , llist )

inputFileName  = None
outputFileName = None
  
if len(argv) == 1:
  print('Please specify input file')
  
if len(options) == 0:
  inputFileName  = argv[1]
  options, remainder = getopt.getopt( argv[2:] , slist, llist )
    
for opt, arg in options:      
  if opt in ('-i', '--input'):
    inputFileName  = arg
  elif opt in ('-o', '--output'):
    outputFileName  = arg
  elif opt in ('-h', '--help'):
    printHelp()
    sys.exit()
      
if outputFileName == None:
  outputFileName = inputFileName[:-4] + '.dat'
  
mesh = meshio.read(inputFileName,file_format="gmsh")

of = open(outputFileName, "w")

printNodes( of , mesh )

printElements( of , mesh )

printGroups( of , mesh )

printRest( of , mesh )

print("The mesh described in '%s' has been saved as '%s'." %(inputFileName,outputFileName))
