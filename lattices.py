import numpy as np
import networkx as nx
import itertools as it
from pairlist import pairs_iter
import random
# 拡張ice ruleを満たすような氷を生成するプログラムが必要。

def simplecubic(N):
    rpos = np.array([[x,y,z] for x,y,z in it.product(range(N),range(N),range(N))]) / N
    cellmat = np.diag([N,N,N]) * 1.0
    return rpos, cellmat

def bodycenter(N):
    rpos = np.array([[x,y,z] for x,y,z in it.product(range(N),range(N),range(N))])
    rpos = np.vstack([rpos, rpos+0.5]) / N
    cellmat = np.diag([N,N,N]) * 2 / 3**0.5
    return rpos, cellmat

def facecenter(N):
    rpos = np.array([[x,y,z] for x,y,z in it.product(range(N),range(N),range(N)) if (x+y+z)%2 == 0]) / N
    cellmat = np.diag([N,N,N]) / 2**0.5
    return rpos, cellmat

def diamond(N):
    rpos, cellmat = facecenter(N)
    rpos *= N
    rpos = np.vstack([rpos, rpos+0.5]) / N
    cellmat *= (8/3)**0.5
    return rpos, cellmat

def find_cycle(graph, history):
    """
    Recursively find a cycle.
    """
    last = history[-1]
    if len(history) > 1:
        last2 = history[-2]
    else:
        last2 = -1
    nexts = list(graph.neighbors(last))
    for next in random.sample(nexts, len(nexts)):
        if next != last2:
            # going forward
            if next in history:
                # a shortcut loop
                head = history.index(next)
                yield history[head:]
            yield from find_cycle(graph, history + [next])


def delete_cycle(g, cy):
    for i in range(len(cy)):
        g.remove_edge(cy[i-1], cy[i])
    for node in cy:
        if len(list(g.neighbors(node))) == 0:
            g.remove_node(node)


def icerule(g):
    # とりあえず、transversing cycleは気にしない。
    # gは破壊される
    dg = nx.DiGraph()
    while len(g) > 0:
        start = random.choice(list(g.nodes()))
        print(start)
        for cycle in find_cycle(g, [start]):
            break
        print(cycle, len(list(g.neighbors(start))))
        nx.add_cycle(dg, cycle)
        delete_cycle(g, cycle)
    return dg

if __name__ == "__main__":
    rpos, cellmat = simplelattice(8)
    g = nx.Graph([[i,j] for i,j in pairs_iter(rpos, 1.1, cellmat, distance=False)])
    dg = icerule(g)
    print(dg.edges())
