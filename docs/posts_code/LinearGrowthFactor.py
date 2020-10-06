#%%
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
plt.style.use('./Modified5308.mplstyle') 
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
figtransparent = 1
#%%
def g_flat(ai, Omi):
    Ol = 1-Omi # flat Universe of only mater and lambda
    hubble = lambda a:  np.sqrt(Omi/a**3 + Ol)
    rhs_g = lambda a: 1/a**3/hubble(a)**3
    gint = quad(rhs_g,0,ai)[0]
    g_f = 2.5*Omi*hubble(ai)/ai
    return g_f*gint

def g_fit_flat(ai, Omi):
# Fitting formula for flat Universes Carroll, Press & Turner
    Ol = 1- Omi
    Om = Omi/ai**3
    return (5*Om/(2*(Om**(4/7)-Ol+(1+Om/2)*(1+Ol/70))))*ai

#%%
g_flat_vec = np.vectorize(g_flat)

a = np.logspace(-2,0,50)
plt.figure(figsize=(10,8))
plt.plot(1/a, g_flat_vec(a,1)*a, label=r"$\Omega_m=1, \, \Omega_\Lambda=0$")
plt.plot(1/a, g_flat_vec(a,0.4)*a, label=r"$\Omega_m=0.4, \, \Omega_\Lambda=0.6$")
plt.plot(1/a, g_flat_vec(a,0.3)*a, label=r"$\Omega_m=0.3, \, \Omega_\Lambda=0.7$")
#plt.plot(1/a, g_fit_flat(a,0.3)*a, '-', label=r"fit $\Omega_m=0.3, \, \Omega_\Lambda=0.7$", lw=2)
plt.plot(1/a, g_flat_vec(a,0.2)*a, label=r"$\Omega_m=0.2, \, \Omega_\Lambda=0.8$")
plt.plot(1/a, g_flat_vec(a,0.1)*a, label=r"$\Omega_m=0.1, \, \Omega_\Lambda=0.9$")
plt.plot(1/a, g_flat_vec(a,0.01)*a, label=r"$\Omega_m=0.01, \, \Omega_\Lambda=0.99$")
plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.xlim(10,1)
plt.ylim(1e-1,1)
plt.xlabel(r"$1+z$")
#plt.xticks([7,5,3,1], ["7","5","3","1"])
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
plt.gca().xaxis.set_minor_formatter(FormatStrFormatter('%d'))
plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%d'))
plt.gca().yaxis.set_minor_formatter(FormatStrFormatter('%.1f'))
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.title("growth factor for flat cosmologies")
#plt.show()
plt.savefig("../assets/growthfactor-inhomogeneous-as-function-of-a.png", dpi=150,transparent=figtransparent)
# %%
