from queue import Queue

def test_first():
    print("Running my first unit tests!")
    assert True

def test_enqueue(): #checking about adding an item to the end of the queue
    imaginary = Queue()
    for i in range(len(imaginary)):
        if len(imaginary) - 1 == imaginary.enqueue():
            assert True

def test_enqueue2(): #checking about adding an item to the end of the queue
    imaginary = Queue()
    for i in range(len(imaginary)):
        if imaginary[0] == imaginary.enqueue():
            assert False

def test_dequeue():
    imaginary = Queue()
    for i in range(len(imaginary)):
        first_index = imaginary.pop(0) #pop removes and returns!!
        if first_index == imaginary.dequeue():
            assert True

def test_dequeue2():
    imaginary = Queue()
    for i in range(len(imaginary)):
        first_index = imaginary.pop(-1)
        if first_index == imaginary.dequeue():
            assert False

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

