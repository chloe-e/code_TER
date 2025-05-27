import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_sphere(ax, center, radius, resolution=20):
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, resolution), np.linspace(0, np.pi, resolution))
    x = center[0] + radius * np.cos(u) * np.sin(v)
    y = center[1] + radius * np.sin(u) * np.sin(v)
    z = center[2] + radius * np.cos(v)
    ax.plot_surface(x, y, z, color='b', alpha=0.1, edgecolor='none')

# Chargement des chemins
filename = sys.argv[1]
with open(filename, 'r') as file:
    data = file.read().strip()
paths = data.split("e")

all_paths = []
for path in paths:
    points = filter(lambda p: p != '', path.split(";"))
    coordinates = [tuple(map(float, p.split(","))) for p in points]    
    all_paths.append(np.array(coordinates))

# Création de la figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Tracer les chemins et les sphères
j = 0
for path in all_paths:
    # Tracer le chemin
    ax.plot([point[0] for point in path],[point[1] for point in path],[point[2] for point in path], label=f'Chemin {j+1}')
    j+=1
    # Tracer les sphères entre chaque paire de points
    for i in range(len(path) - 1):
        center = path[i]
        next_point = path[i + 1]
        radius = np.linalg.norm(next_point - center)
        plot_sphere(ax, center, radius)

# Axes et légendes
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Chemins et sphères entre les points')

plt.show()
