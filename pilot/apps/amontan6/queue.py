import collections
from typing import Any, Deque

class Queue:

    def __init__(self) -> None:
        """Initializes a new, empty queue."""
        # The internal data structure is a deque, optimized for this purpose.
        self._items: Deque[Any] = collections.deque()

    def enqueue(self, item: Any) -> None:
        """Adds an item to the back of the queue."""
        self._items.append(item)

    def dequeue(self) -> Any:
        """Removes and returns the item from the front of the queue."""
        if not self._items:
            raise IndexError("Cannot dequeue from an empty queue")
        return self._items.popleft()

    def size(self) -> int:
        """Returns the total number of items in the queue."""
        return len(self._items)

    def is_empty(self) -> bool:
        """Checks if the queue contains no items."""
        return len(self._items) == 0

    def __len__(self) -> int:
        """Allows `len(queue)` to return the size of the queue."""
        return self.size()

    def __repr__(self) -> str:
        """Provides a string representation for debugging."""
        return f"Queue({list(self._items)})"

