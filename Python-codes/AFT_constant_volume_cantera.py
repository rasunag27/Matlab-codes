import cantera as ct
import matplotlib.pyplot as plt
import numpy as np


gas = ct.Solution('gri30.xml')
phi = np.arange(0.5,1.5,0.1)
n = len(phi)

for i in range(n):
	gas.TPX = 298.15, 101325, {'CH4':1,'O2':2,'N2':7.25}
	gas.set_equivalence_ratio(phi[i],'CH4','O2:1.0,N2:3.76')
	gas.equilibrate('UV','auto')
	plt.plot(phi[i],gas.T,'b*--',markersize=8)

plt.xlabel('Equivalence ratio (phi)')
plt.ylabel('Adiabatic flame temparature (K)')
plt.title('Equivalence ratio vs AFT')
plt.grid(True)
plt.show()
hold(True)


