#!/usr/bin/env python
import numpy as np
from pprint import pprint
from scipy.optimize import linprog

from ex07 import nn2na, get_usage

"""
#EX8
A technology company has 3 production plants where it produces 2 chip models. It has 2 stock plants and 3 
sale points that requires different quantities of each model. The costs, production and demand related 
data:

 1) Find the optimum distribution for the period using the scipy.linprog package
"""

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

# GENERAMOS LOS VERTICES
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

dims = len(vertices)

NN = np.zeros((dims,dims),int)                     # la matriz se inicializa en 0s
C = np.zeros(len(arcos), int)                      # inicializo C en zero
b = Production.flatten()                           # Aplano la producción
beq = Demand.flatten() * -1                        # Aplano la demanda y la hago negativa

for i,n in enumerate(arcos):
    NN[vertices.index(n[0]), vertices.index(n[1])] = 1   # agrego arco
    C[i] = n[2]                                          # agrego el costo

if __name__ == '__main__':
    A, arcs = nn2na(NN)                # arma nodo array
    Aeq = A.copy()                     # Me armo una copia para eq

    A[6:] = np.zeros(A.shape[1])       # Solo dejamos los arcos de P
    Aeq[:6] = np.zeros(A.shape[1])     # Sólo dejamos los arcos de S y V

    b = np.pad(b, (0,10))              # expando b

    print("**** A = b ****")
    for i, a in enumerate(A):
       print(f"{a} = {b[i]}")
    print("FIN\n\n")

    beq = np.pad(beq, (10, 0))         # expando beq

    print("\n\n**** Aeq = beq ****")
    for i, a in enumerate(Aeq):
       print(f"{a} = {beq[i]}")
    print("FIN\n\n")

    # definimos las cotas entre 0 e infinito
    l = np.zeros(A.shape[1])
    u = np.full(A.shape[1], None)
    bounds = tuple(zip(l, u))

    res = linprog(C, A_eq=Aeq, b_eq=beq, bounds=bounds, A_ub=A, b_ub=b, method='revised simplex')

    usage = get_usage([(x[0], x[1]) for x in arcos ], res.x.astype(float), [x[2] for x in arcos ])

    total_cost = res.fun
    pprint(res)

    print("## Results ##")
    print(f"Usage: (from,to) -> used/max: {usage}")
    print(f"Min cost: {total_cost:0.2f}")
