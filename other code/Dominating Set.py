import networkx as nx
import pandas as pd
import networkx.algorithms.approximation as nxaa

# build up a graph

edge=pd.read_csv("forode\\network2.csv")

edge=edge.apply(tuple, axis=1).tolist()

print(edge)

G = nx.Graph()

G.add_edges_from(edge)

nx.draw_networkx(G)

# output
print(nx.dominating_set(G),"\n" ) # return a dominating sets


print("min_dominating_set:",nxaa.min_weighted_dominating_set(G, weight=None),"\n")


