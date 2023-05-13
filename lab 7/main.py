import sys
import time
import matplotlib.pyplot as plt
import random

class Graph:
    def __init__(self, vertices):
        self.V = vertices
        self.graph = [[0 for column in range(vertices)]
                    for row in range(vertices)]

    def add_edge(self, u, v, w):
        self.graph[u][v] = w
        self.graph[v][u] = w


class Prim(Graph):
    def __init__(self, vertices):
        super().__init__(vertices)

    def primMST(self):
        key = [sys.maxsize] * self.V
        parent = [None] * self.V
        key[0] = 0
        mstSet = [False] * self.V

        parent[0] = -1

        for cout in range(self.V):
            u = self.minKey(key, mstSet)
            mstSet[u] = True

            for v in range(self.V):
                if (self.graph[u][v] > 0 and mstSet[v] == False and key[v] > self.graph[u][v]):
                    key[v] = self.graph[u][v]
                    parent[v] = u

        self.printMST(parent)

    def minKey(self, key, mstSet):
        min = sys.maxsize
        for v in range(self.V):
            if key[v] < min and mstSet[v] == False:
                min = key[v]
                min_index = v

        return min_index

    def printMST(self, parent):
        print("Edge \tWeight")
        for i in range(1, self.V):
            print(parent[i], "-", i, "\t", self.graph[i][parent[i]])


class Kruskal:
    def __init__(self, vertices):
        self.V = vertices
        self.graph = []

    def add_edge(self, u, v, w):
        self.graph.append([u, v, w])

    def find(self, parent, i):
        if parent[i] == i:
            return i
        return self.find(parent, parent[i])

    def union(self, parent, rank, x, y):
        xroot = self.find(parent, x)
        yroot = self.find(parent, y)

        if rank[xroot] < rank[yroot]:
            parent[xroot] = yroot
        elif rank[xroot] > rank[yroot]:
            parent[yroot] = xroot
        else:
            parent[yroot] = xroot
            rank[xroot] += 1

    def kruskalMST(self):
        result = []
        i, e = 0, 0
        self.graph = sorted(self.graph, key=lambda item: item[2])
        parent = []
        rank = []

        for node in range(self.V):
            parent.append(node)
            rank.append(0)

        while e < self.V - 1:
            u, v, w = self.graph[i]
            i = i + 1
            x = self.find(parent, u)
            y = self.find(parent, v)

            if x != y:
                e = e + 1
                result.append([u, v, w])
                self.union(parent, rank, x, y)

        self.printMST(result)

    def printMST(self, result):
        print("Edges in the constructed MST")
        for u, v, weight in result:
            print("%d -- %d == %d" % (u, v, weight))


def analyze_algorithms():
    vertices_input = input("Enter the number of vertices: ")
    try:
        v = int(vertices_input)
    except ValueError:
        print("Invalid input! Please enter an integer.")
        return

    vertices = range(10, v, 20)
    times_prim = []
    times_kruskal = []

    for v in vertices:
        graph_prim = Prim(v)
        graph_kruskal = Kruskal(v)
        for i in range(v):
            for j in range(i+1, v):
                weight = random.randint(1, 100)
                graph_prim.add_edge(i, j, weight)
                graph_kruskal.add_edge(i, j, weight)

        start_time = time.time()
        graph_prim.primMST()
        times_prim.append(time.time() - start_time)

        start_time = time.time()
        graph_kruskal.kruskalMST()
        times_kruskal.append(time.time() - start_time)

    plt.plot(vertices, times_prim, label='Prim')
    plt.plot(vertices, times_kruskal, label='Kruskal')
    plt.xlabel('Vertices')
    plt.ylabel('Time (s)')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    analyze_algorithms()