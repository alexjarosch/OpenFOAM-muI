import numpy as np

""" This is a quick Python script to calculate IN1 as well as A_minus in the
    same way the muIreg viscosity module is doing. I use a coding style very
    close to the C++ syntax of OpenFOAM.
"""

# parameters from the tutorial case "damBreak_muIreg"
mus = 0.342
mud = 0.557
muInf = 0.05
alphaReg = 1.9
I0 = 0.069

# I use Z instead of the lowercase zeta in Barker & Gray 2017
tanZs = np.tan((np.arctan(mus) + np.arctan(mud))/2)
# lowest I value, close to zero
Ilower = 1e-9
# upper I value
Iupper = (I0*(tanZs - mus)/(mud - tanZs))
# I for now, I call all values reg here
Ireg = (Ilower + Iupper)/2

# solve the equation in a while loop
while(Iupper - Ilower > 1e-9):

    Ireg = (Ilower + Iupper)/2

    muIreg = (mus + Ireg*(mud - mus)/(Ireg + I0))
    # muPrime in Barker & Gray 2017, derivative of mu with respect to I
    muPrime = ((I0*(mud - mus))/(Ireg + I0)/(Ireg + I0))
    # Inupnu is the fraction used in the C eq 3.9 in Barker & Gray
    Inupnu = muPrime*Ireg/muIreg
    # C equation 3.9 in Barker & Gray 2017
    C = 4*Inupnu*Inupnu - 4*Inupnu + muIreg*muIreg*(1 - Inupnu/2)*(1 - Inupnu/2)

    if C < 0:
        Iupper = Ireg
    else:
        Ilower = Ireg

# the estimated I becomes I^N_1
IN1 = Ireg
# A minus according to eq 6.4
A_minus = IN1*np.exp(alphaReg*np.power(I0 + IN1, 2)/np.power(mus*I0 + mud*IN1 + muInf*IN1*IN1, 2))

print("IN1    :", IN1)
print("A_minus:", A_minus)
