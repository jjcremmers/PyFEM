################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2024. The code is written in 2011-2012 by                #
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

class itemList ( dict ):

    """
    Class to construct a list of items that have a coninuous local number, and
    a global ID.
    """
  
    def add ( self, ID: int, item ):
    
        """
        Adds an item with an ID to the list. This item will be stored in the list.
    
        Args:
            ID (int): the ID of the item to be stored.
            item: the value(s) of the item to be stored.
        """
        
        if ID in self:
            raise RuntimeError( 'ID ' + str(ID) + ' already exists in ' + type(self).__name__ )

        self[ID] = item

    def get ( self, IDs ):

        """
        Returns the index / indices of an ID or list of IDs of items in the list.
    
        Args:
            IDs (list[int]|int,optional): the ID/IDs. If ommited, a list with all indces
            will be returned.
        Returns:
            list[int]: a list with the indices. In the case of a single ID, this list has
            length 1.
        """
        
        if isinstance(IDs,int):
            return self[IDs]
        elif isinstance(IDs,list):
            return [ self[ID] for ID in IDs ]
      
        raise RuntimeError('illegal argument for itemList.get')
    
    def getIndices ( self, IDs : list[int] | int = -1 ) -> list[int]:
        
        """
        Returns the index / indices of an ID or list of IDs of items in the list.
    
        Args:
            IDs (list[int]|int,optional): the ID/IDs. If ommited, a list with all indces
            will be returned.
        Returns:
            list[int]: a list with the indices. In the case of a single ID, this list has
            length 1.
        """
        
        if IDs == -1:
            return list(self.keys())
        elif isinstance(IDs,int):
            return list(self.keys()).index( IDs )
        elif isinstance(IDs,list):
            return [ list(self.keys()).index( ID ) for ID in IDs ]
      
        raise RuntimeError('illegal argument for itemList.getIndices')  
    
    def findID( self , index : int ) -> int:
  
        """
        Returns the ID of an index in the list.
    
        Args:
            index (int): the index of the item
            
        Returns:
            int: the ID of the item
        """
        
        return list(self.keys())[index]
