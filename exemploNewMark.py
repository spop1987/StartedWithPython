from ctypes import sizeof
import math as m
import numpy as np
import matplotlib.pyplot as plt

def Newmark_coefficients(dt):
    alpha = 0.25
    delta = 0.5

    a0 = 1/(alpha*(dt**2))
    a1 = 1/(alpha*dt)
    a2 = (1/2*alpha) - 1
    a3 = (1 - delta)*dt
    a4 = delta*dt
    a5 = a4*a0
    a6 = a4*a1 - 1
    a7 = a4*a2 - a3

    return a0, a1, a2, a3, a4, a5, a6, a7

def MechanicalProperties(mass, k, zeta, nof, m):
    
    if nof == 1:
        M = mass
        K = k
        Wn = m.sqrt(K/M)
        C = 2*M*zeta*Wn
    else:
        M = np.ones((nof,nof), float)
        K = np.ones((nof,nof), float)
        C = np.ones((nof,nof), float)
    
    return M,K,C

def CalculateConst(vetorConstantes, MatrizM, MatrizK, MatrizC):
    return (vetorConstantes[0]*MatrizM + vetorConstantes[5]*MatrizC + MatrizK)**(-1)

constantes = np.zeros(8, float)
dt = 0.01
NGDL = 1

constantes = Newmark_coefficients(dt)

tempo = np.arange(0, 20, dt)
print(tempo)
M, K, C = MechanicalProperties(10, 40, 0.1, NGDL, m)

comprimentoTempo = int(len(tempo))

Desloc = np.zeros(NGDL, comprimentoTempo)
Desloc[:,0] = 0.5
print(Desloc)
Veloc = np.zeros(NGDL, comprimentoTempo)
F = np.zeros(shape=(NGDL, comprimentoTempo))

A_0 = (M**(-1))*F[:,0] - C*Veloc[:,0] - K*Desloc[:,0]
Acel = np.zeros(NGDL, comprimentoTempo)
Acel[:,0] = A_0

constante = CalculateConst(constantes, M, K, C)

i = 0
for t in tempo:
    c1 = M*(constantes[0]*Desloc[:,i] + constantes[1]*Veloc[:,i] + constantes[2]*Acel[:,i])
    c2 = C*(constantes[5]*Desloc[:,i] + constantes[6]*Veloc[:,i] + constantes[7]*Acel[:,i])
    Desloc[:,i+1] = constante*( F[:,i+1] + c1 + c2 )
    Veloc[:,i+1] = constantes[5]*( Desloc[:,i+1] - Desloc[:,i] ) - constantes[6]*Veloc[:,i] - constantes[7]*Acel[:,i]
    Acel[:,i+1] = constantes[0]*( Desloc[:,i+1] - Desloc[:,i] ) - constantes[1]*Veloc[:,i] - constantes[2]*Acel[:,i]


# plt.plot(tempo, Desloc)




