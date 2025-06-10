
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]
file = open(filename,'r')
data = file.read().strip()
paths = data.split("e")


all_paths = []
for path in paths:
    points = (path.split(";"))
    points = filter(lambda p : p != '',points)
    coordinates = [tuple(map(float,p.split(","))) for p in points]    
    all_paths.append(np.array(coordinates))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
print(all_paths)

# Tracer chaque chemin
for i, path in enumerate(all_paths):
    ax.plot([point[0] for point in path],[point[2] for point in path], label=f'Chemin {i+1}')

# Activer l'interactivité
ax.set_xlabel('X')
ax.set_zlabel('Z')
ax.set_title('Représentation des chemins')

plt.legend()
plt.show()
file.close()
