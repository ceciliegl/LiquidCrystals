import numpy as np

N = 1000
cLL = 0.1

f = np.ones(N)/(2*np.pi)

dang = 2*np.pi/N

E = np.sum(f*np.log(2*np.pi*f))*dang

for i in range(N):
    for j in range(N):
        E += 0.5*cLL*f[i]*f[j]*np.abs(np.sin(dang*abs(i-j)))*dang*dang

print(E)
