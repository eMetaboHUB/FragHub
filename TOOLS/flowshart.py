from graphviz import Digraph

# Créez un objet Digraph
dot = Digraph()

# Ajoutez des nœuds au graphe
dot.node('A', 'Nœud A')
dot.node('B', 'Nœud B')
dot.node('C', 'Nœud C')

# Ajoutez des arêtes entre les nœuds
dot.edges(['AB', 'AC'])

# Vous pouvez également personnaliser l'apparence du graphe
dot.attr(rankdir='LR')  # Définir l'orientation de gauche à droite

# Enregistrez le graphe au format PDF (ou tout autre format souhaité)
dot.render('mon_graphe', format='pdf')