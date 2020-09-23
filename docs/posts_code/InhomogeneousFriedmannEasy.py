# %% 
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
plt.style.use('./Modified5308.mplstyle') 
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

figtransparent = 0

#%%
# Relative Change of Volume PLOT
def rhs_del(a,t, Om, zi, delta_opzi):
    opzi = 1+ zi
    ai = 1/opzi
    Ok0 = (delta_opzi/opzi*(ai**3*(-6 + delta_opzi/opzi)*(-1 + Om) - (-15 + delta_opzi/opzi)*Om))/(9.*ai)
    OL0 = 1-Om-Ok0 
    return np.sqrt(Om/a*(1+delta_opzi/opzi)- Ok0 + OL0*a**2) 

def a_from_t_for_delta(delta, lrhs):
    a0 = [1e-5] # at a=1 we want da/dt=1 we can later get units from H0.
    zi = 1/a0[0]-1
    t = np.logspace(-8,1, 190000)
    ts=t
    sol = odeint(lrhs, a0, ts, args= (Om, zi, delta) )
    ts = ts[np.isfinite(sol.flatten())]
    sol = sol[np.isfinite(sol.flatten())].flatten()
    return ts, interp1d(ts,sol,bounds_error=None, fill_value="extrapolate")

def rhs(a,t, Om, OL0):
    return np.sqrt(Om/a-(1-Om-OL0) + OL0*a**2)

def rhs_sph(a,t, Om, zi, delta):
    Ok0 = (1+zi)*Om*delta *5/3
    OL0 = 1-Om-Ok0 
    return np.sqrt(Om/a*(1+delta)-Ok0 + OL0*a**2) #/bH

def rhsq(a,t, Om, zi, delta_opzi):
    abdot = rhs_del(a,t, Om, zi, 0)
    addot = rhs_del(a,t, Om, zi, delta_opzi)
    return 



Om = 0.3
ts, a_from_t = a_from_t_for_delta(0., rhs_del)

d1 = 0.5 # three different deltas and their negative values
d2 = 1
d3 = 7
ts, ad1 = a_from_t_for_delta(d1, rhs_del)
ts, mad1 = a_from_t_for_delta(-d1, rhs_del)
ts, ad2 = a_from_t_for_delta(d2, rhs_del)
ts, mad2 = a_from_t_for_delta(-d2, rhs_del)
ts, ad3 = a_from_t_for_delta(d3, rhs_del)
ts, mad3 = a_from_t_for_delta(-d3, rhs_del)

plt.plot(1/a_from_t(ts), (mad3(ts)/a_from_t(ts))**(3)-1,'--',c=cycle[2],label=r"$\delta/a_i=-{:.1g}$".format(d3))
plt.plot(1/a_from_t(ts), (mad2(ts)/a_from_t(ts))**(3)-1,'--',c=cycle[1],label=r"$\delta/a_i=-{:.1g}$".format(d2))
plt.plot(1/a_from_t(ts), (mad1(ts)/a_from_t(ts))**(3)-1,'--',c=cycle[0],label=r"$\delta/a_i=-{:.1g}$".format(d1))
plt.plot(1/a_from_t(ts), (ad1(ts)/a_from_t(ts))**(3)-1,c=cycle[0],label=r"$\delta/a_i={:.1g}$".format(d1))
plt.plot(1/a_from_t(ts), -1*((mad1(ts)/a_from_t(ts))**(3)-1),'--',c=cycle[0],alpha=.9,lw=1)

plt.plot(1/a_from_t(ts), (ad2(ts)/a_from_t(ts))**(3)-1,c=cycle[1],label=r"$\delta/a_i={:.1g}$".format(d2))
plt.plot(1/a_from_t(ts), -1*((mad2(ts)/a_from_t(ts))**(3)-1),'--',c=cycle[1],alpha=.9,lw=1)

plt.plot(1/a_from_t(ts), (ad3(ts)/a_from_t(ts))**(3)-1,c=cycle[2],label=r"$\delta/a_i={:.1g}$".format(d3))
plt.plot(1/a_from_t(ts), -1*((mad3(ts)/a_from_t(ts))**(3)-1),'--',c=cycle[2],alpha=.9,lw=1)

plt.plot(1/a_from_t(ts), -((d1)*a_from_t(ts)),'-.',c=cycle[0],alpha=.9,lw=2,label=r"$\propto a$")
plt.plot(1/a_from_t(ts), (d1*a_from_t(ts)),'-.',c=cycle[0],alpha=.9,lw=2)
plt.plot(1/a_from_t(ts), -(d2*a_from_t(ts)),'-.',c=cycle[1],alpha=.9,lw=2)
plt.plot(1/a_from_t(ts), (d2*a_from_t(ts)),'-.',c=cycle[1],alpha=.9,lw=2)
plt.plot(1/a_from_t(ts), -(d3*a_from_t(ts)),'-.',c=cycle[2],alpha=.9,lw=2)
plt.plot(1/a_from_t(ts), (d3*a_from_t(ts)),'-.',c=cycle[2],alpha=.9,lw=2)

plt.xscale("log")
plt.xlim(1000,1)
plt.ylim(-1.25,2.25)
plt.xlabel(r"$1+z$")
plt.ylabel(r"${\delta V/V=(a_\delta}/{a_b})^3-1$")
plt.title(r"change of density vs background $\Omega_m={:.1f}$".format(Om))
plt.legend()
plt.savefig("../assets/deltaVoV-evolution-inhomogeneous-as-function-of-opz.png", dpi=150,transparent=figtransparent)


# %%
## Growth Factor comparison
a0 = [1e-3] # at a=1 we want da/dt=1 we can later get units from H0.
zi = 999
t = np.linspace(1e-6,5, 190000)
ts=t

def g_flat(ai, Omi):
    """Linear growth factor for a flat Universe"""
    Ol = 1-Omi # flat Universe of only matter and lambda
    hubble = lambda a:  np.sqrt(Omi/a**3 + Ol)
    rhs_g = lambda a: 1/a**3/hubble(a)**3
    gint = quad(rhs_g,0,ai)[0]
    g_f = 2.5*Omi*hubble(ai)/ai
    return g_f*gint

def rhsf(a,t,Om):
    return np.sqrt(Om/a+(1-Om)*a**2)
sol = odeint(rhsf, a0, ts, args= (Om, ) )
ts = ts[np.isfinite(sol.flatten())]
sol = sol[np.isfinite(sol.flatten())].flatten()
t_from_a = interp1d(sol,ts,bounds_error=None,fill_value="extrapolate")


ctime = t_from_a(1)
ft = np.array([1e-9,t_from_a(0.1),t_from_a(0.5), ctime])
delt = -1*np.logspace(-3,1,100)[::-1]
delt = np.append(delt,-delt[::-1])

fdel = np.array([])
fdelz1= np.array([])
fdelz9= np.array([])
for delta in delt:
    sol = odeint(rhs_del, 1e-19, ft, args= (Om, 999, delta))
    fdelz9 = np.append(fdelz9,sol[1]/0.1)
    fdelz1 = np.append(fdelz1,sol[2]/0.5)
    fdel = np.append(fdel,sol[3]/1)

plt.plot(np.abs(delt),fdel**(3)-1,c=cycle[0],label=r"$z=0$")
plt.plot(np.abs(delt),fdelz1**(3)-1,c=cycle[1],label=r"$z=1$")
plt.plot(np.abs(delt),fdelz9**(3)-1,c=cycle[2],label=r"$z=9$")

plt.plot(np.abs(delt), 1/(1*delt*g_flat(1 , Om)+1)-1,'--', lw=2, c=cycle[0],alpha=.7,label=r"$aD_1(z)\delta_0$")
plt.plot(np.abs(delt), 1/(1+.5*delt*g_flat(.5 , Om))-1,'--', lw=2, c=cycle[1],alpha=.7)
plt.plot(np.abs(delt), 1/(1+.1*delt*g_flat(.1 , Om))-1,'--', lw=2, c=cycle[2],alpha=.7)

plt.xlabel(r"$|\delta_0|=|\delta/a_i|$")
plt.ylabel(r"$\delta V/V=(a_\delta/a_b)^3-1$")
plt.title("relative change of volume vs initial density")
plt.xscale("log")
plt.ylim(-.75,.75)
plt.xlim(9e-3,1e1)
plt.legend()
#plt.savefig("../assets/deltaVoV-evolution-inhomogeneous-as-function-of-delta0.png", dpi=150,transparent=figtransparent)
# %%

## Growth Factor comparison
a0 = [1e-3] # at a=1 we want da/dt=1 we can later get units from H0.
zi = 999
t = np.linspace(1e-6,5, 190000)
ts=t

def g_flat(ai, Omi):
    """Linear growth factor for a flat Universe"""
    Ol = 1-Omi # flat Universe of only matter and lambda
    hubble = lambda a:  np.sqrt(Omi/a**3 + Ol)
    rhs_g = lambda a: 1/a**3/hubble(a)**3
    gint = quad(rhs_g,0,ai)[0]
    g_f = 2.5*Omi*hubble(ai)/ai
    return g_f*gint

def rhsf(a,t,Om):
    return np.sqrt(Om/a+(1-Om)*a**2)
sol = odeint(rhsf, a0, ts, args= (Om, ) )
ts = ts[np.isfinite(sol.flatten())]
sol = sol[np.isfinite(sol.flatten())].flatten()
t_from_a = interp1d(sol,ts,bounds_error=None,fill_value="extrapolate")

aini = 1e-4
ctime = t_from_a(1)
ft = np.array([aini,t_from_a(0.1),t_from_a(0.5), ctime])
delt = np.linspace(-1,1,500)
#delt = -1*np.logspace(-1.5,1,100)[::-1]
#delt = np.append(delt,-delt[::-1])

fdel = np.array([])
fdelz1= np.array([])
fdelz9= np.array([])
for delta in delt:
    fadi = (1+delta*aini)**(-1/3)
    sol = odeint(rhs_del, 1e-11, ft, args= (Om, 1/aini-1, delta))
    fdelz9 = np.append(fdelz9,sol[1]/0.1)*fadi
    fdelz1 = np.append(fdelz1,sol[2]/0.5)*fadi
    fdel = np.append(fdel,sol[3]/1)*fadi

deltp = delt


plt.plot(deltp,fdel**(3)-1,c=cycle[0],label=r"$z=0$")
plt.plot(deltp,fdelz1**(3)-1,c=cycle[1],label=r"$z=1$")
plt.plot(deltp,fdelz9**(3)-1,c=cycle[2],label=r"$z=9$")

plt.plot(deltp, 1/(1*delt*g_flat(1 , Om)+1)-1,'--', lw=2, c=cycle[0],alpha=.7,label=r"$(1+aD_1(z)\delta_0)^{-1}-1$")
plt.plot(deltp, 1/(.5*delt*g_flat(.5 , Om)+1)-1,'--', lw=2, c=cycle[1],alpha=.7)
plt.plot(deltp, 1/(.1*delt*g_flat(.1 , Om)+1)-1,'--', lw=2, c=cycle[2],alpha=.7)

plt.xlabel(r"$|\delta_0|=|\delta/a_i|$")
plt.ylabel(r"$\delta V/V=(a_\delta/a_b)^3-1$")
plt.title("relative change of volume vs initial density")
#plt.xscale("log")
plt.ylim(-.75,.75)
plt.xlim(-.75,.75)
plt.legend()
plt.savefig("../assets/deltaVoV-evolution-inhomogeneous-as-function-of-delta0.png", dpi=150,transparent=figtransparent)

# %%
# Now scale out the redshift
deltp = delt
l = (1/(1*delt*g_flat(1 , Om)+1))**(1/3)
lz1 = (1/(.5*delt*g_flat(.5 , Om)+1))**(1/3)
lz9 = (1/(.1*delt*g_flat(.1 , Om)+1))**(1/3)

fi1 = (1+delt/1e5)**(-1/3)
plt.plot(deltp,(fdel**(1))/l-1,c=cycle[0],label=r"$z=0$")
plt.plot(deltp,(fdelz1**(1))/lz1-1,c=cycle[1],label=r"$z=1$")
plt.plot(deltp,(fdelz9**(1))/lz9-1,c=cycle[2],label=r"$z=9$")

plt.xlabel(r"$|\delta_0|=|\delta/a_i|$")
plt.ylabel(r"$\delta V/V=(a_\delta/a_b)^3-1$")
plt.title("relative change of volume vs initial density")
#plt.xscale("log")
plt.ylim(-.0025,.0025)
plt.xlim(-.075,.075)
plt.legend()
#plt.savefig("../assets/deltaVoV-evolution-inhomogeneous-as-function-of-delta0.png", dpi=150,transparent=figtransparent)
# %%
plt.loglog(np.abs(deltp),np.abs(fdel**(3)-1),np.abs(deltp),np.abs(deltp))
# %%
plt.loglog(np.abs(deltp),np.abs(((fdel**(1))/l-1)),np.abs(deltp),np.abs(deltp))
# %%
