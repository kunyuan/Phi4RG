import numpy as np
import matplotlib.pyplot as plt

UV=8.0

Bin=64

Bare=16

Scale=[(i*1.0+1.0)*UV/Bin for i in range(Bin)]
# Scale=[1.0/(Bin-i)*UV for i in range(Bin)]

dScale=[Scale[i]-Scale[i-1] for i in range(1, Bin)]

print Scale
print dScale

EffVer=[Bare for i in range(Bin)]

d=-3.0/16/np.pi

for i in range(Bin-1):
    start=Bin-1-i
    end=Bin-2-i
    print start, end, dScale[end], Scale[start]
    # EffVer[end]=(EffVer[start]*Scale[start]-3.0/16.0/np.pi*EffVer[start]**2*dScale[end])/Scale[end]
    dg=(-EffVer[start]-d*EffVer[start]**2)*dScale[end]/Scale[start]
    EffVer[end]=EffVer[start]-dg

plt.plot(Scale, EffVer, 'o-')
plt.show()
