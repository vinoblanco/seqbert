import time

class Node:
  def __init__(self, data):
    self.data = data
    self.next = None

class LinkedList:
    def __init__(self):
        self.head = None
        self.tail = None

    def append(self, value):
        new_node = Node(value)

        if self.head is None:
            # Liste ist leer, Head und Tail sind derselbe Node
            self.head = new_node
            self.tail = new_node
        else:
            # Den aktuellen Tail auf den neuen Knoten verlinken
            self.tail.next = new_node
            # Tail nachziehen
            self.tail = new_node

    def iterate(self):
        node = self.head
        while node is not None:
            print(node.data)
            node = node.next

    def lenght(self):
        count = 0
        node = self.head
        while node is not None:
            node = node.next
            count += 1
        return count

    def pairwise(self):
        node = self.head
        while node and node.next:
            yield node, node.next
            node = node.next