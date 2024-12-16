import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
data = pd.read_csv("clusters.csv")

# Define a color palette (noise is assigned a neutral color)
num_clusters = data['cluster'].max()
palette = sns.color_palette("hsv", num_clusters)
palette.insert(0, (0.5, 0.5, 0.5))  # Gray for noise (cluster 0)

# Plot the clusters
plt.figure(figsize=(8, 6))
sns.scatterplot(
    x="x",
    y="y",
    hue="cluster",
    palette=palette,
    data=data,
    legend="full",
    s=100
)

# Add plot details
plt.title("Cluster Visualization")
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
plt.legend(title="Cluster", loc="upper right", bbox_to_anchor=(1.15, 1))
plt.grid(True)
plt.tight_layout()

# Show the plot
plt.show()
