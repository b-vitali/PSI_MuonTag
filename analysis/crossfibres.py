import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

#colours = ['blue','green','C2','C3','C4','C5','C6','C7','C8','C9']
custom_colors = ['green', 'blue', 'green', 'purple', 'orange']
custom_markerstyles = ['o', 's', '^', 'x', 'd']
custom_linestyles = ['-', '--', '-.', ':']

step5 = False
dx = 1
dy = 10
limit = 2

f_dy = []
f_dx = []
def YError(phi):
    #return (dx / np.sin((phi)*np.pi/180))
    Dyx = (dx / np.sin((phi)*np.pi/180))                                                                                                                 
    Dyy = (dy / np.cos((phi)*np.pi/180)) 
    for i in range(len(phi)):
        values = [Dyx[i],Dyy[i]]
        f_dy.append(min(values))        
    return f_dy    

def XError(phi):
    #return (dx / np.cos((phi)*np.pi/180))
    Dxx= (dx / np.cos((phi)*np.pi/180))                                                                                                                 
    Dxy= (dy / np.sin((phi)*np.pi/180)) 
    for i in range(len(phi)):
        values = [Dxx[i],Dxy[i]]
        f_dx.append(min(values))                                                                                                                           
    return f_dx     

t = np.arange(0.1, 89.9, 0.1)
Dy = YError(t)
Dx = XError(t)
s = np.sin(2*np.pi*t)

Dy = np.array(Dy)
Dx = np.array(Dx)

fig, ax = plt.subplots()
plt.title('Projected uncertainties for rotated fibres', fontsize=15)
plt.ylabel('$dx; dy$ [mm]', fontsize=15)
plt.xlabel(r'$\vartheta$ [deg]', fontsize=15)
plt.xticks(np.arange(0, 95, step=5))
plt.yticks(np.arange(1, 11, step=0.5))
plt.yscale('log')
plt.grid(True, which="both")
plt.plot(t, Dx, color=custom_colors[0], linestyle=custom_linestyles[0], label=r"$dx = min (\delta x\ \sec\phi;\delta y\ \csc\phi)$")
plt.plot(t, Dy, color=custom_colors[1], linestyle=custom_linestyles[0], label=r"$dy = min (\delta x\ \csc\phi;\delta y\ \sec\phi)$")
#ax.fill_between(t, min(Dx), max(Dx), where=(Dy < limit)&(Dx < limit), facecolor='green', alpha=0.3, label = r'$\delta x='+str(dx)+r'; \delta y='+str(dy)+r'\rightarrow dx;dy<'+str(limit)+'$')
plt.legend(loc="upper center", fontsize = 15)
plt.show()

#reset
f_dy = []
f_dx = []
DT = 1
z_max = 15
r = 3.5
def anglesD1Turn(a):
    t = np.tan(a*np.pi/180)*z_max/(2*np.pi*r)
    b = np.arctan(2*np.pi*r* (t+DT)/z_max)*180/np.pi
    return b, t

if step5:
    a = np.arange(0, 60, 5)
    b, t = anglesD1Turn(a)
    b = np.floor(b/5)*5
    tb = np.tan(b*np.pi/180)*z_max/(2*np.pi*r)
    print(tb-t)

else:
    a = np.arange(0, 60, 0.1)
    b, t = anglesD1Turn(a)
    tb = np.tan(b*np.pi/180)*z_max/(2*np.pi*r)
    print(tb-t)

fig0, ax0 = plt.subplots()
plt.title('Turns vs Angles', fontsize=15)
plt.ylabel("Turns", fontsize=15)
plt.xlabel(r'$\vartheta$ [deg]', fontsize=15)
plt.xticks(np.arange(0, 95, step=5))
plt.yscale('log')
plt.grid(True, which="both")
plt.plot(a, t, color=custom_colors[0], linestyle=custom_linestyles[0])
#ax.fill_between(a, min(b), max(b), where=(a < 45), facecolor='green', alpha=0.3)
plt.show()

fig, ax = plt.subplots()
plt.title('Angles: 1 turn difference', fontsize=15)
plt.ylabel(r'$\vartheta$ [deg]', fontsize=15)
plt.xlabel(r'$\vartheta_{out}$ [deg]', fontsize=15)
plt.xticks(np.arange(0, 95, step=5))
plt.yticks(np.arange(0, 95, step=5))
plt.grid(True, which="both")
if step5:
    plt.plot(a, b, color=custom_colors[1],  marker=custom_markerstyles[0], linestyle='None', label=r'$\vartheta_{in}$')
    plt.plot(a, abs(a-b), color=custom_colors[0],  marker=custom_markerstyles[0], linestyle='None', label=r'$\vartheta_{in}-\vartheta_{out}$')
else:
    plt.plot(a, b, color=custom_colors[1], linestyle=custom_linestyles[0], label=r'$\vartheta_{in}$')
    plt.plot(a, abs(a-b), color=custom_colors[0], linestyle=custom_linestyles[0], label=r'$\vartheta_{in}-\vartheta_{out}$')
#ax.fill_between(a, min(b), max(b), where=(a < 45), facecolor='green', alpha=0.3)
plt.legend(loc="best", fontsize = 15)
plt.show()

dxa = XError(a)
dyb = YError(b)
dxa = np.array(dxa)
dyb = np.array(dyb)

fig2, ax2 = plt.subplots()
plt.title('Projected uncertainties: 1 turn difference', fontsize=15)
plt.xlabel(r'$\vartheta_{out}$ [deg]', fontsize=15)
plt.xticks(np.arange(0, 95, step=5))
#plt.yticks(np.arange(1, 11, step=0.5))
plt.yticks(np.arange(1, 11, step=0.05))
#plt.yscale('log')
plt.grid(True, which="both")
plt.ylabel('$dx; dy$ [mm]', fontsize=15)
if step5:
    plt.plot(a, dxa, color=custom_colors[0],  marker=custom_markerstyles[0], linestyle='None', label=r'$dx_{out}$')
    plt.plot(a, dyb, color=custom_colors[1],  marker=custom_markerstyles[0], linestyle='None', label=r'$dy_{in}$')
else:
    plt.plot(a, dxa, color=custom_colors[0], linestyle=custom_linestyles[0], label=r'$dx_{out}$')
    plt.plot(a, dyb, color=custom_colors[1], linestyle=custom_linestyles[0], label=r'$dy_{in}$')  
plt.legend(loc="best", fontsize = 15)
#ax2.fill_between(a, min(dxa), max(dxa), where=(dxa < 2)&(dyb<2), facecolor='green', alpha=0.3)
plt.show()

def Length(phi):
    t = np.tan(phi*np.pi/180)*z_max/(2*np.pi*r)
    base = t*2*np.pi*r
    height = z_max
    l = np.sqrt(base*base+height*height)
    return l

La = Length(a)
Lb = Length(b)

fig3, ax3 = plt.subplots()
plt.title("Fibres' length: 1 turn difference", fontsize=15)
plt.xlabel(r'$\vartheta_{out}$ [deg]', fontsize=15)
plt.xticks(np.arange(0, 95, step=5))
#plt.yticks(np.arange(1, 11, step=0.5))
#plt.yscale('log')
plt.grid(True, which="both")
plt.ylabel('Fibre length [cm]', fontsize=15)
if step5:
    plt.plot(a, La, color=custom_colors[0],  marker=custom_markerstyles[0], linestyle='None', label=r'$L_{out}$')
    plt.plot(a, Lb, color=custom_colors[1],  marker=custom_markerstyles[0], linestyle='None', label=r'$L_{in}$')
else:
    plt.plot(a, La, color=custom_colors[0], linestyle=custom_linestyles[0], label=r'$L_{out}$')
    plt.plot(a, Lb, color=custom_colors[1], linestyle=custom_linestyles[0], label=r'$L_{in}$')
plt.legend(loc="best", fontsize = 15)
#ax2.fill_between(a, min(dxa), max(dxa), where=(dxa < 2)&(dyb<2), facecolor='green', alpha=0.3)
plt.show()
