"""
Evaluation of Adiabatic flame temparature for hydrocarbons in Constant pressure reactor with heat loss

The general hydrocarbon stiochiometric equation is:

CxHy + (x+y/4)*(O2+3.76*N2) = xCO2 + (y/2)H2O + 3.76*(x+y/4)*N2
Root problem for constant pressure: H_products - H_reactants + LHV*heat_loss
LHV = Low heat value of the gas
Heat_loss = heat loss factor = 0.35

"""

import matplotlib.pyplot as plt
import math
import numpy as np

def h(T,co_effs):

	""" 
	Entering the coefficient values and temp for polynomial
	function.

	"""
	R =8.314 #J/mol.K
	a1 = co_effs[0];
	a2 = co_effs[1];
	a3 = co_effs[2];
	a4 = co_effs[3];
	a5 = co_effs[4];
	a6 = co_effs[5];

	return(a1 + a2*T/2 + a3*pow(T,2)/3 + a4*pow(T,3)/4 + a5*pow(T,4)/5 + a6/T )*R*T

"""
Let the number of carbon atoms be 2. The comparison is done for three hydrocarbons
1. C2H2 , LHV = 1250e3
2. C2H4 , LHV = 1480e3
3. C2H6 , LHV = 1435e3

"""

"""
Reactants: C2H2, C2H4, C2H6, O2, N2 - Low temp
Products: CO2, H2O, N2 - High temp

"""
# Reactants
c2h2_coeffs_l = [8.08681094E-01, 2.33615629E-02, -3.55171815E-05, 2.80152437E-08, -8.50072974E-12, 2.64289807E+04, 1.39397051E+01]
c2h4_coeffs_l = [3.95920148E+00, -7.57052247E-03, 5.70990292E-05, -6.91588753E-08, 2.69884373E-11, 5.08977593E+03, 4.09733096E+00]
c2h6_coeffs_l = [4.29142492E+00, -5.50154270E-03, 5.99438288E-05, -7.08466285E-08, 2.68685771E-11, -1.15222055E+04, 2.66682316E+00]
o2_coeffs_l = [3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03, 3.65767573E+00]
n2_coeffs_l = [0.03298677E+02, 0.14082404E-02, -0.03963222E-04, 0.05641515E-07, -0.02444854E-10, -0.10208999E+04, 0.03950372E+02]

# Products
co2_coeffs_h = [3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00]
h2o_coeffs_h = [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00]
n2_coeffs_h = [0.02926640E+02, 0.14879768E-02, -0.05684760E-05, 0.10097038E-09, -0.06753351E-13, -0.09227977E+04, 0.05980528E+02]

def f(T,n):

	""" 
	Root finding problem
	"""
	#Products enthaply

	h_co2_p = h(T,co2_coeffs_h)
	h_h2o_p = h(T,h2o_coeffs_h)
	h_n2_p = h(T,n2_coeffs_h)

	# Reactants enthaply

	Tstd = 298.15
	R = 8.314
	h_c2h2_r = h(Tstd,c2h2_coeffs_l)
	h_c2h4_r = h(Tstd,c2h4_coeffs_l)
	h_c2h6_r = h(Tstd,c2h6_coeffs_l)
	h_o2_r = h(Tstd,o2_coeffs_l)
	h_n2_r = h(Tstd,n2_coeffs_l)


	x = 2;
	
	if n==1:
		y = 2
		H_products = x*h_co2_p + (y/2)*h_h2o_p + 3.76*(x+y/4)*h_n2_p
		H_reactants = h_c2h2_r + (x+y/4)*(h_o2_r + 3.76*h_n2_r)
		LHV = 1250*pow(10,3)
	elif n==2:
		y = 4
		H_products = x*h_co2_p + (y/2)*h_h2o_p + 3.76*(x+y/4)*h_n2_p
		H_reactants = h_c2h4_r + (x+y/4)*(h_o2_r + 3.76*h_n2_r)
		LHV = 1480*pow(10,3)
	elif n==3:
		y = 6
		H_products = x*h_co2_p + (y/2)*h_h2o_p + 3.76*(x+y/4)*h_n2_p
		H_reactants = h_c2h6_r + (x+y/4)*(h_o2_r + 3.76*h_n2_r)
		LHV = 1435*pow(10,3)

	return H_products - H_reactants + LHV*heat_loss

def fprime(T,n):
	dT = 1e-6
	return (f(T+dT,n) - f(T,n))/dT

n = [1,2,3]

T_guess = np.ones(len(n))

for i in range(len(n)):
	T_guess[i] = 1500;

ct = 0;
alpha = 0.5;
tol = 1e-5;
heat_loss = 0.35;

"""
Newton-Raphson method
"""
for i in range(len(n)):
	while(abs(f(T_guess[i],n[i])) > tol):
		T_guess[i] = T_guess[i] - alpha*(f(T_guess[i],n[i])/fprime(T_guess[i],n[i]))
		ct = ct+1;
print(ct)
print(T_guess)
plt.plot(n,T_guess,'r*--',markersize=8)
plt.title('AFT for C2H2, C2H4, C2H6')
plt.xlabel('Hydrocarbons(CxHy)')
plt.ylabel('Adiabatic flame temparature (K)')
plt.grid(True)
plt.show()