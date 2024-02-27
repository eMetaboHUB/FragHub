from anytree.exporter import JsonExporter
from anytree.importer import JsonImporter
from anytree import Node, RenderTree
import pandas as pd
import json
import networkx as nx
import matplotlib.pyplot as plt

df = pd.read_excel(r"C:\Users\Axel\Documents\PYTHON\FragHub\FragHub\TOOLS\generate_decision_tree\GRAPH_CONSTRUCT_SOLUTIONS_FINAL.xlsx", sheet_name="Shimadzu")

# Copier la colonne 'SOLUTION'
saved_solution = df['SOLUTION'].copy()

# Convertir toutes les colonnes à lower case excepté 'SOLUTION'
df = df.applymap(lambda s: s.lower() if type(s) == str else s)

# Restaurer la colonne 'SOLUTION'
df['SOLUTION'] = saved_solution

# Initialiser le noeud racine
root = Node("root")

# création du dictionnaire pour stocker les références des nœuds
nodes = {"root": root}

for _, row in df.iterrows():
    parent_name = "root"
    for val in row:
        # Concaténer le nom du parent et la valeur actuelle pour créer un identifiant unique
        node_name = f"{parent_name}-{val}"
        if node_name not in nodes:
            nodes[node_name] = Node(node_name, parent=nodes[parent_name])
        parent_name = node_name

for pre, fill, node in RenderTree(root):
    if node.is_leaf:
        # Découpe le nom du noeud en liste à chaque occurrence de "-"
        name_parts = node.name.split('-')
        # Prend seulement les parts après la troisième occurrence de "-"
        leaf_name = '-'.join(name_parts[6:])
        print(f"{pre}{leaf_name}")
    else:
        # Sinon, affiche uniquement la dernière partie du nom du nœud
        print(f"{pre}{node.name.split('-')[-1]}")

exporter = JsonExporter(indent=2, sort_keys=True)

# Export the tree to a JSON string
jsonstring = exporter.export(root)

# Save to a JSON file
with open(r'C:\Users\Axel\Documents\PYTHON\FragHub\FragHub\TOOLS\generate_decision_tree\tree.json', 'w') as f:
    f.write(jsonstring)

with open(r'C:\Users\Axel\Documents\PYTHON\FragHub\FragHub\TOOLS\generate_decision_tree\tree.json') as f:
    data = json.load(f)
root = data
G = nx.DiGraph()

# Note that we must add the root node separately here
G.add_node(root['name'], layer=0)

def add_edges(parent, depth=0):
    for i, child in enumerate(parent.get('children', [])):
        child_name = child['name']
        G.add_node(child_name, layer=depth)
        G.add_edge(parent['name'], child_name)
        add_edges(child, depth=depth + 1)

# Ajout des arêtes à G
add_edges(root)

# Création de la disposition multipartite
pos = nx.multipartite_layout(G, subset_key="layer")

# Inversion des valeurs y de la disposition
pos = {node: (x, -y) for node, (x, y) in pos.items()}

# Dessin du graphe
plt.figure(figsize=(15, 15))

# Dessin des nœuds et des arêtes
nx.draw(G, pos, with_labels=False, arrows=False)

# Préparation des labels
labels = {}
# Pour chaque nœud dans le graphe
for node in G.nodes():
    if node is not None:  # Assurez-vous que le nœud n'est pas None
        if G.out_degree(node) == 0:  # Si le nœud est une feuille
            # Découpez le nom du nœud en liste à chaque occurrence de "-"
            name_parts = node.split('-')
            # Prenez seulement les parts après la troisième occurrence de "-"
            leaf_name = '-'.join(name_parts[6:])
            labels[node] = leaf_name
        else:  # Si le nœud n'est pas une feuille
            labels[node] = node.split('-')[-1]

# Dessin des labels
nx.draw_networkx_labels(G, pos, labels=labels, horizontalalignment='left')

plt.show()