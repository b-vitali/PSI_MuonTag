import numpy as np
import matplotlib.pyplot as plt
colours = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']

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
plt.plot(t, Dx, 'C0-', label=r"$dx = min (\delta x\ \sec\phi;\delta y\ \csc\phi)$")
plt.plot(t, Dy, 'C1-', label=r"$dy = min (\delta x\ \csc\phi;\delta y\ \sec\phi)$")
ax.fill_between(t, min(Dx), max(Dx), where=(Dy < limit)&(Dx < limit), facecolor='green', alpha=0.3, label = r'$\delta x='+str(dx)+r'; \delta y='+str(dy)+r'\rightarrow dx;dy<'+str(limit)+'$')
plt.legend(loc="upper left", fontsize = 15)
plt.show()

