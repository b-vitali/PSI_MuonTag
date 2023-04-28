import numpy as np

def create_helix(height, radius, num_turns):
    z = np.linspace(0, height, num=1000)
    theta = np.linspace(0, num_turns*2*np.pi, num=1000)
    x = radius*np.cos(theta)
    y = radius*np.sin(theta)
    return np.column_stack((x,y,z))

def create_helix_opposite(height, radius, num_turns):
    z = np.linspace(0, height, num=1000)
    theta = np.linspace(0, num_turns*2*np.pi, num=1000)
    x = radius*np.cos(-theta)
    y = radius*np.sin(-theta)
    return np.column_stack((x,y,z))

import math

def create_cylinder(height, radius, num_turns, num_fibers):
    cylinder = []
    for i in range(num_fibers):
        angle = i * 2*math.pi / num_fibers
        fiber = create_helix(height, radius, num_turns)
        rotation_matrix = np.array([[np.cos(angle), -np.sin(angle), 0],
                                    [np.sin(angle), np.cos(angle), 0],
                                    [0, 0, 1]])
        fiber = np.dot(fiber, rotation_matrix)
        cylinder.append(fiber)
    return np.concatenate(cylinder, axis=0)

def create_cylinder_opposite(height, radius, num_turns, num_fibers):
    cylinder = []
    for i in range(num_fibers):
        angle = i * 2*math.pi / num_fibers
        fiber = create_helix_opposite(height, radius, num_turns)
        rotation_matrix = np.array([[np.cos(angle), -np.sin(angle), 0],
                                    [np.sin(angle), np.cos(angle), 0],
                                    [0, 0, 1]])
        fiber = np.dot(fiber, rotation_matrix)
        cylinder.append(fiber)
    return np.concatenate(cylinder, axis=0)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Create cylinder
height = 10
radius = 2
thikness = 0.5
num_turns = 0.5

num_fibers = int(2*radius*np.pi / thikness)

cylinder_in = create_cylinder(height, radius, num_turns, num_fibers)
cylinder_out = create_cylinder_opposite(height, radius+0.2, num_turns, num_fibers)

# Create 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(num_fibers):
    ax.plot(cylinder_in[i*1000:(i+1)*1000, 0], cylinder_in[i*1000:(i+1)*1000, 1], cylinder_in[i*1000:(i+1)*1000, 2], color='blue')
    ax.plot(cylinder_out[i*1000:(i+1)*1000, 0], cylinder_out[i*1000:(i+1)*1000, 1], cylinder_out[i*1000:(i+1)*1000, 2], color='green')

# Set axis labels and limits
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim(-radius, radius)
ax.set_ylim(-radius, radius)
ax.set_zlim(0, height)

# Show plot
plt.show()
