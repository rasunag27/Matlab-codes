"""
Program to evaluate Adiabatic flame temparature for methane with different Equivalence ratio(phi).

Combustion reaction: CH4 + 2(O2+3.76*N2) = CO2 + 2*H20+ 7.52*N2
Lean mixture (phi<1): CH4 + (2/phi)*(O2+3.76N2) = CO2+2H2O+((2-2phi)/phi)*O2+(7.52/phi)*N2
Rich mixture (phi>1): CH4 + (2/phi)*(O2+3.76N2) = (4-4/phi)*CO+(4/phi - 3)*CO2+2H2O+(7.52/phi)*N2

"""
import numpy as np
import matplotlib.pyplot as plt
import math

def h(T,co_effs):

	"""
	Entering the co-efficient values and temp for polynomial function.

	"""

	R = 8.314 #J/mol.K
	a1 = co_effs[0];
	a2 = co_effs[1];
	a3 = co_effs[2];
	a4 = co_effs[3];
	a5 = co_effs[4];
	a6 = co_effs[5];

	return(a1 + a2*T/2 + a3*pow(T,2)/3 + a4*pow(T,3)/4 + a5*pow(T,4)/5 + a6/T )*R*T

"""
Reactants: CH4, O2, N2 - Low temp (T<1000)
Products: CO2, H2O, CO, O2, N2 - High temp (T>1000)

"""

#Reactants
ch4_coeffs_l = [5.14987613E+00, -1.36709788E-02, 4.91800599E-05, -4.84743026E-08, 1.66693956E-11, -1.02466476E+04, -4.64130376E+00]
o2_coeffs_l = [3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03, 3.65767573E+00]
n2_coeffs_l = [0.03298677E+02, 0.14082404E-02, -0.03963222E-04, 0.05641515E-07, -0.02444854E-10, -0.10208999E+04, 0.03950372E+02]

#Products
co2_coeffs_h = [3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00]
h2o_coeffs_h = [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00]
n2_coeffs_h = [0.02926640E+02, 0.14879768E-02, -0.05684760E-05, 0.10097038E-09, -0.06753351E-13, -0.09227977E+04, 0.05980528E+02]
o2_coeffs_h = [3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00]
co_coeffs_h = [2.71518561E+00, 2.06252743E-03, -9.98825771E-07, 2.30053008E-10, -2.03647716E-14, 1.41518724E+04, 7.81868772E+00]

def f(T,phi):

	"""
	Root finding problem with phi as input

	"""
	#Products enthalpy

	h_co2_p = h(T,co2_coeffs_h)
	h_h2o_p = h(T,h2o_coeffs_h)
	h_n2_p = h(T,n2_coeffs_h)
	h_co_p = h(T,co_coeffs_h)
	h_o2_p = h(T,o2_coeffs_h)

	#Reactants enthalphy

	Tstd = 298.15
	R = 8.314
	h_ch4_r = h(Tstd,ch4_coeffs_l)
	h_o2_r = h(Tstd,o2_coeffs_l)
	h_n2_r = h(Tstd,n2_coeffs_l)


	if phi==1:
		H_products = h_co2_p + 2*h_h2o_p + 7.52*h_n2_p
		H_reactants = h_ch4_r + 2*h_o2_r + 7.52*h_n2_r
		N_p = N_r = 10.52

	elif phi<1:
		H_products = h_co2_p + 2*h_h2o_p + (2/phi - 2) * h_o2_p + (7.52/phi)*h_n2_p
		N_p = 3 + (2/phi - 2) + (7.52/phi)
	
	else:
		H_products = (4-4/phi)*h_co_p + (4/phi - 3)*h_co2_p + 2*h_h2o_p + (7.52/phi)*h_n2_p
		N_p = (4-4/phi) + (4/phi - 3) + 2 + (7.52/phi)


	H_reactants = h_ch4_r + (2/phi)*h_o2_r + (7.52/phi)*h_n2_r
	N_r = 1 + (2/phi) + (7.52/phi)

	return H_reactants - H_products - R*(N_r*Tstd - N_p*T)

def fprime(T,phi):
	return (f(T+1e-6,phi) - f(T,phi))/1e-6


ct = 0; # Iteration
alpha = 0.2; # Relaxation parameter

phi = np.arange(0.5,1.5,0.1)
n = len(phi);

T_guess = np.ones(n);

for i in range(n):
	T_guess[i] = 1500;

tol = 1e-5; 

"""
Using Newton-Raphson method

"""
for i in range(n):
	while(abs(f(T_guess[i],phi[i])) > tol):
		T_guess[i] = T_guess[i] - alpha*(f(T_guess[i],phi[i])/fprime(T_guess[i],phi[i]))
		ct = ct+1


print(ct)
print(max(T_guess))
plt.plot(phi,T_guess,'r*--')
plt.xlabel('Equivalence ratio (phi)')
plt.ylabel('Adiabatic flame temparature')
plt.title('Equivalence ratio vs AFT')
plt.grid(True)
plt.show()

