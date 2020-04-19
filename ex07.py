from pprint import pprint

import numpy as np
from scipy.optimize import linprog

"""
#EX0: Build a function that transforms Node-Node matrix (incidence matrix) to Node-Arc matrix
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


vertices = ['s', 'a', 'b', 't']
arcos = [(0,1), (0,2), (1,2), (1,3), (2,3)]
capacidades = [5, 10, 6, 3, 3]
tiempos = [1, 3, 2, 1, 1]

#
tiempo = 6

# The strategy is to convert the graph in a bigger one where each
# Then we can use the linprog or ford - fulkerson

dims = tiempo*len(vertices) + 2          # cada vertice por tiempo + los dos extremos
NN = np.zeros((dims,dims),int)           # la matriz se inicializa en 0s

# Registro cada nodo en una lista de nodos y sus capacidades
capacidades_nueavas = []
nodos = []
todos_los_nodos = ['s']                  # node inicial s
for t in range(tiempo):
  # de s vamos a s(i)
  NN[0][t + 1] = 1
  capacidades_nueavas.append(None)       # los iniciales tienen capacidad infinita

for j,v in enumerate(vertices):
  for t in range(tiempo):
      todos_los_nodos.append(v+str(t+1))                          # voy agregando a mi lista de labels de nodos
      for k, a in enumerate(arcos):
              if a[0] == j:
                  if t + tiempos[k] < tiempo:
                      ta = t+1
                      v1 = f"{vertices[a[0]]}{ta}"                # construyo el vertice nuevo
                      v2 = f"{vertices[a[1]]}{ta+tiempos[k]}"     # el vertice destino corrido en tiempo de ese arco
                      nodos.append((v1, v2))
                      capacidades_nueavas.append(capacidades[k])  # registro la capacidad. es importante conservar el orden

todos_los_nodos.append('t')                                       # agrego el ultimo label de nodo

for n in nodos:
    NN[todos_los_nodos.index(n[0]),todos_los_nodos.index(n[1])] = 1  # pongo en 1 los indices de la matriz donde hay nodo

for t in range(tiempo):
  # de t(i) a t
  NN[dims - tiempo - 1 + t][dims - 1] = 1
  capacidades_nueavas.append(None)                                   # capacidad infinita

# cerramos el ciclo de t a s
NN[dims-1][0] = 1
capacidades_nueavas.append(None)                                     # ya dije capacidad infinta?

if __name__ == '__main__':
    Aeq, arc_idxs = nn2na(NN)
    bounds = tuple([(0, None) for arcs in range(0, Aeq.shape[1])])
    beq = np.zeros(Aeq.shape[0])
    C = np.zeros(Aeq.shape[1])
    C[-1] = -1
    len(arc_idxs)
    l = np.zeros(len(arc_idxs))
    u = capacidades_nueavas
    bounds = tuple(list(zip(l, u)))

    res = linprog(C, A_eq=Aeq, b_eq=beq, bounds=bounds)

    usage = get_usage(arc_idxs, res.x.astype(float), capacidades_nueavas)
    min_cut = min_cut(arc_idxs, res.x, np.array(capacidades_nueavas))
    acum_flow = - res.fun
    pprint(res)
    print("## Results ##")
    print(f"Usage: (from,to) -> used/max: {usage}")
    print(f"Min arcs to be cut to cut: (from,to) {min_cut}" )
    print(f"Max flow: {acum_flow:0.2f}")


