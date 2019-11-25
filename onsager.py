import numpy as np
import scipy.integrate as int

def func1(phi, alpha):
    return np.cosh(alpha*np.cos(phi))

def normalize(func, alpha):
    phis = np.linspace(0, 2*np.pi, 1000)
    N = np.trapz(func(phis, alpha), phis)
    return N

def f(func, phi, alpha):
    A = normalize(func, alpha)
    return func(phi, alpha)/float(A)

def F(func, cLL, alpha):
    phis = np.linspace(0, 2*np.pi, 1000)
    phiss = np.linspace(0, 2*np.pi, 1000)
    int1 = int.quad(f(func, x, alpha)*np.log(2*np.pi*f(func, x, alpha)), 0, 2*pi)
    int2 = 0.5*cLL*np.traps(np.trapz(f(func, phis, alpha)*f(func, phiss, alpha)*np.abs(sin(phis-phiss)), phis), phiss)
    return (int1 + int2)
