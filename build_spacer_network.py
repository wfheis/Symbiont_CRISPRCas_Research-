"""
Will Heisler
Script to build a network of strain genomes based on shared CRISPR array spacers from BLAST results using NetworkX.
11/3/2025
"""

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib

# Load BLAST
blast_tsv = "bb_bh_outgroupRC_blast.tsv"

array_blast = pd.read_csv(blast_tsv, sep="\t", header=None)
array_blast.columns = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

# Filter for hits that are not self-hits and are highly similar
blast_filtered = array_blast[
    (array_blast["qseqid"] != array_blast["sseqid"]) &
    (array_blast["pident"] >= 98) & 
    #(array_blast["pident"] >= 90) & 
    (array_blast["length"] >= 20)
]

# Function to extract strain ID from the sequence ID
def get_strain_id(seq_id):
    return seq_id.split("_")[0]

# Count shared spacers between strains
#strain_pairs = blast_filtered.copy()
#strain_pairs["strain1"] = strain_pairs["qseqid"].apply(get_strain_id)
#strain_pairs["strain2"] = strain_pairs["sseqid"].apply(get_strain_id)

# Group by strain pairs and count shared spacers
#edges = strain_pairs.groupby(["strain1", "strain2"]).size().reset_index(name="weight")

# Create graph
#G = nx.Graph()

all_strains = set(array_blast["qseqid"].apply(get_strain_id)) | \
              set(array_blast["sseqid"].apply(get_strain_id))

# Count shared spacers between strains
strain_pairs = blast_filtered.copy()
strain_pairs["strain1"] = strain_pairs["qseqid"].apply(get_strain_id)
strain_pairs["strain2"] = strain_pairs["sseqid"].apply(get_strain_id)

# Group by strain pairs and count shared spacers
edges = strain_pairs.groupby(["strain1", "strain2"]).size().reset_index(name="weight")
print(edges)
# Save edges to CSV
edges.to_csv("outgroupRC_spacer_edges.csv", index=False)

# Create graph AND ensure isolated nodes are included
G = nx.Graph()
G.add_nodes_from(all_strains)  # ensure single-node strains are present

# Add edges with weights
for _, row in edges.iterrows():
    G.add_edge(row["strain1"], row["strain2"], weight=row["weight"])



# Set node colors based on strain prefix
node_colors = []
for node in G.nodes():
    if node.startswith("bb"):
        node_colors.append("tomato")  ## P.bonniea
    elif node.startswith("bh"): ## P.hayleyella
        node_colors.append("indianred")
    elif node.startswith("bg"): ## B.gladioli
        node_colors.append("pink")
    elif node.startswith("bp"): ## B.plantarii
        node_colors.append("violet")
    elif node.startswith("ccal"): ## C.calidae
        node_colors.append("orange")
    elif node.startswith("cconc"): ## C.concitans
        node_colors.append("darkgoldenrod")
    else:
        node_colors.append("gray")  # fallback for unknown strains

# Normalize edge weights for color mapping
edge_weights = [d["weight"] for (_, _, d) in G.edges(data=True)]
norm = matplotlib.colors.Normalize(vmin=min(edge_weights), vmax=max(edge_weights))
#norm = matplotlib.colors.Normalize(vmin=2, vmax=20)  # fixed scale for better comparison
cmap = plt.cm.viridis
edge_colors = [cmap(norm(w)) for w in edge_weights]

# Draw the graph
pos = nx.spring_layout(G, k=2.2, seed=10)  
plt.figure(figsize=(10, 9))
nx.draw_networkx_nodes(
    G, pos,
    node_color=node_colors,
    node_size=600, 
    #edgecolors="black", 
    #linewidths=0.8,
    alpha=0.9, 
)
nx.draw_networkx_edges(
    G, pos,
    #width=[w/5 for w in edge_weights],
    width= 2,
    edge_color=edge_colors,  
)
nx.draw_networkx_labels(
    G, pos,
    font_size=5.5,
    font_color="black",
    font_family="Arial", 
    font_weight="regular",
)

legend_handles = [
    mpatches.Patch(color="tomato", label="P.bonniea (bb)"),
    mpatches.Patch(color="indianred", label="P.hayleyella (bh)"),
    mpatches.Patch(color="violet", label="B.plantarii (bp)"),
    mpatches.Patch(color="pink", label="B.gladioli (bg)"),
    mpatches.Patch(color="orange", label="C.calidae (ccal)"),
    mpatches.Patch(color="darkgoldenrod", label="C.concitans (cconc)"),
]
plt.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=8)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cbar = plt.colorbar(sm, ax=plt.gca(), shrink=0.7)
cbar.set_label("Number of Shared Spacers", fontsize=10, labelpad=10)

plt.title("Array Spacer Sharing Network")
plt.tight_layout()
#plt.savefig("arrayA_bb+bh+outgroup_network.png", dpi=600, bbox_inches='tight', transparent=True)
plt.show()



