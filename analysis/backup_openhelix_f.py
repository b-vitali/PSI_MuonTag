import numpy as np
import matplotlib.pyplot as plt

z_min = 0
z_max = 15
turns = 0.8 # 0.22 = 15*2 = 30; 0.8 = 45*2 = 90
thikness = 10 #mm
r = 3
#hits = [[0.5,10],[3.7,10],[7.3,10], [0.5,175],[3.7,175],[7.3,175], [0.5,88],[3.7,88],[7.3,88]]
hits = [[7.3,10],[9,175],[12,88]]

colours = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']

# Connect this to the thickness of the fibres
def round_to(n, precision): 
    correction = 0.5 if n >= 0 else -0.5
    return int( n/precision+correction ) * precision



def line(z,turns,hit):
    m = 360/z_max * turns
    #print(m)
    t = []
    
    #print(hit[0],hit[1])
    q = hit[1] - hit[0]*m
    #print(q)
    q = round_to(q,0.05)
    #print("rounded q ", q)
    #print(q,m)
    t = z*m + q
    #print("raw", t)
    for i in range(len(t)):
        while t[i] < 0:
            #print("low")
            t[i]=t[i]+360
        while t[i] > 360:
            #print("high")
            t[i]=t[i]-360
    #print("wrapped", t)
    #print(np.abs(np.diff(t)))
    #print(len(np.abs(np.diff(t))))
    #print(len(t))
    t[:-1][np.abs(np.diff(t)) >= 200] = np.nan
    return t

x = np.arange(0, z_max, 0.0001)

angle = np.arctan(2*np.pi*r*turns/z_max)*180/np.pi
print(angle)
plt.title('Study for ghost hits: #turns = ' + str(turns) +'; angle = '+str(round(angle,1))+' [deg]')
plt.ylabel('#theta [deg]')
plt.xlabel('z [cm]')
plt.xlim(z_min, z_max)
plt.ylim(0, 360)

for h in range(len(hits)):
    hit = hits[h]
    t1 = line(x, turns, hit)
    t2 = line(x, -turns, hit)

    plt.plot(x, t1, colours[h]+'--')
    plt.plot(x, t2, colours[h]+'-.')

    Dt = []
    for i in range(len(t1)):
        Dti = t1[i]-t2[i]
        Dt.append(Dti)
    
    #print(np.sign(Dt))
    #idx = np.argwhere(np.sign(Dt)==0).flatten()
    Dt_array = np.array(Dt)
    #print(Dt_array)
    idx = np.argwhere(np.abs(Dt_array) < 1).flatten()

    start = 0
    stop = 1  
    idx_short = []
    for i in range(len(idx)):
        if i == len(idx)-1:        
            #print(i)
            idx_short.append(idx[int((stop+start)/2)])
        elif idx[i+1] == idx[i]+1:
            stop = stop+1
        else:
            #print(idx[int((stop+start)/2)] , int((stop+start)/2))
            idx_short.append(idx[int((stop+start)/2)])
            start = stop - 1
        #print(i,start,stop)

    #print(idx)
    #print(idx_short)

    t1_array = np.array(t1)
    print("Hits: ", x[idx_short], t1_array[idx_short])
    #print(t1_array[idx])

    plt.plot(hit[0], hit[1], 'k*', markersize=25)
    plt.plot(x[idx_short], t1_array[idx_short], colours[h]+'o', markersize=10)
plt.show()
