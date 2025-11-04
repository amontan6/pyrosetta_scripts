from queue import Queue

def test_first():
    print("Running my first unit tests!")
    assert True


def test_enqueue(): #checking about adding an item to the end of the queue
    imaginary = Queue()
    items = ["a","b"]
    initial_length = len(imaginary._items)
    
    for i in range(len(imaginary)):
        imaginary.enqueue(items)
        final_length = len(imaginary._items)
        assert final_length > initial_length

#def test_enqueue2(): #checking about adding an item to the front of the queue
#   imaginary = Queue()
#   items = ["a","b"]
#   imaginary.enqueue(items)
#   first_index = items[0]
#   assert imaginary._items[first_index] == items

def test_dequeue():
    imaginary = Queue()
    items = ["a","b"]
    initial_length = len(imaginary._items)

    for i in range(len(imaginary)):
        imaginary.dequeue(items)
        final_length = len(imaginary._items)
        assert final_length < initial_length

#def test_dequeue2():
#    imaginary = Queue()
#    for i in range(len(imaginary)):
#        first_index = imaginary.pop(-1)
#        if first_index == imaginary.dequeue():
#            assert False

"""From here, I need to come back and redo these."""

def test_size():
    imaginary = Queue()
    for i in range(len(imaginary)):
        if len(imaginary) == imaginary.size():
            assert True

def test_size2():
    imaginary = Queue()
    for i in range(len(imaginary)):
        if len(imaginary) != imaginary.size():
            assert False

def test_isempty():
    imaginary = Queue()
    for i in range(len(imaginary)):
        if len(imaginary) == 0:
            assert True

def test_isempty2():
    imaginary = Queue()
    for i in range(len(imaginary)):
        if len(imaginary) != 0:
            assert False

