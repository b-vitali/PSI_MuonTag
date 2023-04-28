import numpy as np
import matplotlib.pyplot as plt

def create_mesh(nz, nt):
    z = np.linspace(0, 10, nz)
    t = np.linspace(0, 360, nt)
    zv, tv = np.meshgrid(z, t)
    print(zv)
    print(tv)
    return zv, tv

def bold(h):
    for i in range(len(h)):
        skip = False
        for j in range(len(h[0])):
            if skip:
                if j< len(h[0])-1:
                    j = j+1
            skip = False
            #neighbour = h[i][j] + h[i][j] + h[i][j] + h[i][j] + h[i][j] + h[i][j] + h[i][j] + h[i][j]
            if h[i][j]==1:
                h[i][j-1]=1
                if j< len(h[0])-1:
                    h[i][j+1]=1
                    skip = True
    return h

def line(z,turns,hit):
    m = 360/10 * turns
    t = []
    z_max = 10
    t_min = 0
    t_max = 360
    
    print(hit[0],hit[1])
    q = hit[0] - hit[1]*m
    print(q,m)
    for zi in z:
        ti = zi*m + q
        #print(ti)
        if ti < t_min:
            #print("low")
            ti=ti+t_max
        if ti > t_max:
            #print("high")
            ti=ti-t_max
        #print(ti)
        t.append(ti)
    return t
hits = [[70,2],[120,5],[70,8]]
print(hits)

nt = 361
nz = 1001
turns = 1.5
tax = np.linspace(0, 360, nt)
z = np.linspace(0, 10, nz)

t_list = [] # t_list[0] = [t_plus , t_minus]
h_list = [] # h_list[0] = [h_plus , h_minus]
for hit in hits:
    t_plus = line(z, turns, hit)
    t_minus = line(z, -turns, hit)

    h_plus, x, y, im = plt.hist2d(z, t_plus, bins = (z, tax), cmap=plt.cm.Greys)
    h_minus, x, y, im = plt.hist2d(z, t_minus, bins = (z, tax), cmap=plt.cm.Greys)

    h_plus = bold(h_plus)
    h_minus = bold(h_minus)
    t_list.append([t_plus,t_minus])
    h_list.append([h_plus,h_minus])

for i in range(len(hits)):
    plt.figure(figsize=(1, 3))
    plt.subplot(131)
    plt.imshow(h_list[i][0], cmap='gray', interpolation='none', extent=[0, 360, 1000, 0])
    plt.subplot(132)
    plt.imshow(h_list[i][1], cmap='gray', interpolation='none', extent=[0, 360, 1000, 0])
    plt.subplot(133)
    h = h_list[i][0]*h_list[i][1]
    plt.imshow(h, cmap='gray', interpolation='nearest', extent=[0, 360, 1000, 0])
    plt.show()

'''
print(h1)
plt.figure()
plt.imshow(h1)
plt.show()

print(h2)
plt.figure()
plt.imshow(h2)
plt.show()

print(h2*h1)
h = h1 + h2
plt.figure()
plt.imshow(h)
plt.show()
'''