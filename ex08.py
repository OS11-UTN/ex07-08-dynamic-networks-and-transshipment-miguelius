#!/usr/bin/env python
import numpy as np
from pprint import pprint
from scipy.optimize import linprog

"""
#EX09: direct LP approach
The shortest path between node s and t has to be found. Each arc has a distance 
but also the time that it takes to travel between its corresponding vertices.

1) If a person had to travel between s and t in less than 9 hours (T). 
   What’s the shortest path? Try to solve the problem with a simple LP model.
2) What if the maximum available time that this person has drops to 8 hours? 
   What’s the new shortest path? Understand the LP model outputs.
3) What’s the first solution that comes to your mind in order to solve point 2 issues? 
   Is it feasible in reality?
"""
def nn2na(NN):
  idx = np.argwhere(NN)
  NA = np.zeros([NN.shape[0], idx.shape[0]]).astype(int)
  for i, arc in enumerate(idx):
    NA[arc[0], i] = 1
    NA[arc[1], i] = -1

  arc_idx = [ (arc[0], arc[1]) for arc in idx]
  return NA, arc_idx

def get_usage(arc_idxs, use, max_flow):
  return [f"{x} -> {np.round(use[i])} / {max_flow[i]}" for i, x in enumerate(arc_idxs)]

def min_cut(arc_idxs, use, max_flow):
  return list(filter(lambda x: x is not None,
                     [x if max_flow[i] != None and np.isclose(use[i], max_flow[i]) == [True] else None for i, x in
                      enumerate(arc_idxs)]))

def get_selected_arcs(arc_idxs, selected_arcs):
  arc = []
  for idx, i in enumerate(selected_arcs):
      if round(i) == 1:
          arc.append(arc_idxs[idx])
  return arc

def get_arcs_as_tuple_list(NN):
    return [ tuple(x) for x in (np.transpose(np.nonzero(NN)).tolist())]

stocks = 2

Plantas_Stock= np.array([
    [100, 200],
    [150, 150],
    [200, 100]
])

Stocks_Puntos= np.array([
    [100, 200],
    [150, 150],
    [200, 100]
])

Production = np.array([
    [20, 30],
    [10, 40],
    [30, 10]
])

Demand = np.array([
    [30, 40],
    [10, 20],
    [20, 20]
])

#vertices = ['P1a', 'P1b', 'P2a', 'P2b','P3a', 'P3b','S1a', 'S1b','S2a', 'S2b','V1a', 'V1b', 'V2a', 'V2b',' V3a', 'V3b']
vertices = []
arcos = []
ps = []
eses = []
vs = []
for p in range(Plantas_Stock.shape[0]):
        for prod in range(Plantas_Stock.shape[1]):
            e = f"P{p+1}{chr(97+prod)}"
            vertices.append(e)
            ps.append(e)

for s in range(stocks):
    for prod in range(Plantas_Stock.shape[1]):
        e = f"S{s + 1}{chr(97 + prod)}"
        vertices.append(e)
        eses.append(e)

for d in range(Plantas_Stock.shape[0]):
    for prod in range(Plantas_Stock.shape[1]):
        e = f"V{d + 1}{chr(97 + prod)}"
        vertices.append(e)
        vs.append(e)

for i in range(Plantas_Stock.shape[0]):
    for j in range(Plantas_Stock.shape[1]):
        for k in range(stocks):
            e = f"P{i + 1}{chr(97 + j)}"
            d = f"S{k + 1}{chr(97 + j)}"
            arcos.append((e,d,Plantas_Stock[i][j],0,Production[i][j]))

for i in range(Stocks_Puntos.shape[0]):
    for j in range(Stocks_Puntos.shape[1]):
        for k in range(stocks):
            s = f"S{k + 1}{chr(97 + j)}"
            v = f"V{i + 1}{chr(97 + j)}"
            arcos.append((s,v,Stocks_Puntos[i][j],0,-Demand[i][j]))

#pprint(vertices)
#pprint(arcos)

dims = len(vertices)
#print(dims)
NN = np.zeros((dims,dims),int)                     # la matriz se inicializa en 0s
NNeq = np.zeros((dims,dims),int)                   # la matriz se inicializa en 0s
C = np.zeros(len(arcos), int)
b = np.array([20, 30, 10, 40, 30, 10])
beq = np.array([-30, -40, -10, -20, -20, -20])

for i,n in enumerate(arcos):
    # print(n[0][0])
    if n[0][0] == 'P':
      NN[vertices.index(n[0]), vertices.index(n[1])] = 1  # pongo en 1 los indices de la matriz donde hay nodo
    if n[1][0] == 'V':
      NNeq[vertices.index(n[0]), vertices.index(n[1])] = 1  # pongo en 1 los indices de la matriz donde hay nodo
    C[i] = n[2]
#print(c1)
#exit(0)
if __name__ == '__main__':
    A, arcs = nn2na(NN)
    # Aprima va a tener los arcos de Stock
    Aprima = np.delete(A,slice(0,6), 0)
    Aprima = np.pad(Aprima, [(6,0), (0,12)]) #
    print(Aprima)
    A = np.delete(A,slice(6,10), 0)
    A = np.pad(A, [(0,4), (0,12)], mode="constant", constant_values=0)
    #print(A)
    #print(A.shape)
    b = np.pad(b, (0,10))
    #print(b.shape)
    print("**** A = b ****")
    for i, a in enumerate(A):
       print(f"{a} = {b[i]}")
    print("FIN\n\n")
    Aeq, arcs_eq = nn2na(NNeq)
    Aeq = np.pad(Aeq, [(0,0), (12,0)] , mode="constant", constant_values=0)
    # Junto los vertices de stock en la matriz Aeq
    Aeq = np.add(Aeq, Aprima)
    beq = np.pad(beq, (10, 0))

    print("\n\n**** Aeq = beq ****")
    for i, a in enumerate(Aeq):
       print(f"{a} = {beq[i]}")
    print("FIN\n\n")

    l = np.zeros(A.shape[1])
    u = np.full(A.shape[1], None)
    bounds = tuple(zip(l, u))

    res = linprog(C, A_eq=Aeq, b_eq=beq, bounds=bounds, A_ub=A, b_ub=b, method='revised simplex')

    usage = get_usage([(x[0], x[1]) for x in arcos ], res.x.astype(float), [x[2] for x in arcos ])
    #min_cut = min_cut(arcos, res.x, np.array(capacidades_nueavas))
    total_cost = res.fun
    pprint(res)
    print("## Results ##")
    print(f"Usage: (from,to) -> used/max: {usage}")
    print(f"Min cost: {total_cost:0.2f}")
