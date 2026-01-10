# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from typing import Any, Union


class itemList(dict):
    """
    A dictionary-based container that maintains items with both global IDs and
    continuous local indices.

    This class extends the built-in dict to provide a mapping from global identifiers
    (IDs) to items, and also provides methods to query items by their position in the
    dictionary. It is commonly used in finite element analysis to manage collections
    of elements, nodes, or other entities that require both global and local numbering.

    Attributes:
        Items are stored with integer IDs as keys and arbitrary values as items.

    Examples:
        Create and populate an itemList:
        >>> items = itemList()
        >>> items.add(1, "element_1")
        >>> items.add(5, "element_5")
        >>> items.get(1)
        'element_1'
        >>> items.findID(0)
        1

    Raises:
        RuntimeError: If attempting to add an item with a duplicate ID, or if
            invalid arguments are provided to methods.
    """

    def add(self, ID: int, item: Any) -> None:
        """
        Add an item with an associated ID to the list.

        If the ID already exists in the itemList, a RuntimeError is raised to
        prevent accidental overwriting of existing items.

        Args:
            ID (int): The unique global identifier for the item. Must be unique
                within this itemList.
            item (Any): The value or object to store. Can be any type (int, str,
                object instance, etc.).

        Raises:
            RuntimeError: If the ID already exists in the itemList.

        Examples:
            >>> items = itemList()
            >>> items.add(10, "data")
            >>> items.add(10, "other")  # doctest: +SKIP
            Traceback (most recent call last):
                ...
            RuntimeError: ID 10 already exists in itemList
        """
        if ID in self:
            raise RuntimeError('ID ' + str(ID) + ' already exists in ' + type(self).__name__)

        self[ID] = item

    def get(self, IDs: Union[int, list]) -> Any:
        """
        Retrieve item(s) from the list by ID or list of IDs.

        Args:
            IDs (int | list[int]): A single ID (returns one item) or a list of
                IDs (returns a list of items).

        Returns:
            Any: If IDs is an int, returns the single item. If IDs is a list,
                returns a list of items corresponding to the requested IDs.

        Raises:
            RuntimeError: If IDs is neither an int nor a list.
            KeyError: If a requested ID does not exist in the itemList.

        Examples:
            >>> items = itemList()
            >>> items.add(1, "first")
            >>> items.add(2, "second")
            >>> items.get(1)
            'first'
            >>> items.get([1, 2])
            ['first', 'second']
        """
        if isinstance(IDs, int):
            return self[IDs]
        elif isinstance(IDs, list):
            return [self[ID] for ID in IDs]

        raise RuntimeError('illegal argument for itemList.get')

    def getIndices(self, IDs: Union[list, int] = -1) -> Union[list, int]:
        """
        Get the local index/indices corresponding to one or more global IDs.

        Returns the position(s) of items in the dictionary's key sequence. Useful
        for converting between global IDs and continuous local numbering (0, 1, 2, ...).

        Args:
            IDs (int | list[int], optional): A single ID, list of IDs, or -1 to
                retrieve all indices. Defaults to -1 (all indices).

        Returns:
            list[int]: If IDs is -1, returns list of all keys. If IDs is an int,
                returns the index (int) of that ID. If IDs is a list, returns a
                list of indices for those IDs.

        Raises:
            RuntimeError: If IDs is not of type int, list, or -1.
            ValueError: If a requested ID does not exist.

        Examples:
            >>> items = itemList()
            >>> items.add(100, "first")
            >>> items.add(200, "second")
            >>> items.getIndices()
            [100, 200]
            >>> items.getIndices(100)
            0
            >>> items.getIndices([100, 200])
            [0, 1]
        """
        if IDs == -1:
            return list(self.keys())
        elif isinstance(IDs, int):
            return list(self.keys()).index(IDs)
        elif isinstance(IDs, list):
            return [list(self.keys()).index(ID) for ID in IDs]

        raise RuntimeError('illegal argument for itemList.getIndices')

    def findID(self, index: int) -> int:
        """
        Retrieve the global ID of an item at a given local index.

        This is the inverse operation of getIndices() for a single index, converting
        from continuous local numbering (0, 1, 2, ...) back to global IDs.

        Args:
            index (int): The local index (position) of the item in the ordered
                sequence of items.

        Returns:
            int: The global ID of the item at the specified index.

        Raises:
            IndexError: If the index is out of range.

        Examples:
            >>> items = itemList()
            >>> items.add(10, "first")
            >>> items.add(20, "second")
            >>> items.findID(0)
            10
            >>> items.findID(1)
            20
        """
        return list(self.keys())[index]
