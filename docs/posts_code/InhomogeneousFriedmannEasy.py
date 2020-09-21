# %% 
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
plt.style.use('./Modified5308.mplstyle') 
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

figtransparent = 0

class cosmo:
    """Keep parameters for a Universe with curvature
    H0 = Hubble constant today: in km/s/Mpc
    Om0: Omega_matter at current time
    Ok0: Omega_curvature today set to 1-Omega_matter-Omega_Lambda
    OL0: Omega_lambda 
     """
    def __init__(self,params):
        self.H0 = params['H0']
        self.Om0 = params['Om0']
        self.OL0 = params['OL0']
        self.Ok0 = 1.0 - self.Om0 - self.OL0
        
    def hubble(self,a):
        return np.sqrt(self.H0**2*(self.Om0/a**3-self.Ok0/a**2+(1-self.Om0)))

    def hubble_from_z(self,z):
        return hubble(self,1/(1+z))

params = {'H0':67.11, 'Om0':0.3175, 'OL0':0.6825}

u1 = cosmo( {'H0':67.11, 'Om0':0.3175, 'OL0':0.6825})
u2 = cosmo( {'H0':70, 'Om0':0.3175, 'OL0':0.6825})
u3 = cosmo( {'H0':67.11, 'Om0':0.28, 'OL0':0.6825})
u4 = cosmo( {'H0':67.11, 'Om0':0.4, 'OL0':0.6825})
universe = cosmo(params)

ar = np.logspace(-2,.4,400)

plt.figure(figsize=(10,8))
plt.loglog(ar,ar*u1.hubble(ar), label=r"$H_0={:.1f},\ \Omega_m$={:.2f}".format(u1.H0,u1.Om0),lw=4)
plt.loglog(ar,ar*u2.hubble(ar), '--',label=r"$H_0={:.1f},\ \Omega_m$={:.2f}".format(u2.H0,u2.Om0),alpha=.7)
plt.loglog(ar,ar*u3.hubble(ar), '--',label=r"$H_0={:.1f},\ \Omega_m$={:.2f}".format(u3.H0,u3.Om0),alpha=.7)
plt.loglog(ar,ar*u4.hubble(ar), '--',label=r"$H_0={:.1f},\ \Omega_m$={:.2f}".format(u4.H0,u4.Om0),alpha=.7)

plt.plot([1],[u1.H0],'o',ms=8,alpha=.75, label="today")

plt.xlim(1e-2,3)
plt.ylim(1e1,5e2)

plt.legend()
plt.ylabel(r"$\dot{a}=a\,H(a)$",)
plt.xlabel(r"$a$")
plt.title("spatially flat cosmologies")
#plt.savefig("../assets/inhomogeneous-Hubble-as-function-of-a.png", dpi=150,transparent=figtransparent)
# %% [markdown]
# <h2>Integrate</h2>
#

# %%
def rhs_del(a,t, Om, zi, delta_opzi):
    opzi = 1+ zi
    ai = 1/opzi
    Ok0 = Om*delta_opzi
    Ok0 = (delta_opzi/opzi*(ai**3*(-6 + delta_opzi/opzi)*(-1 + Om) - (-15 + delta_opzi/opzi)*Om))/(9.*ai)
    OL0 = 1-Om*(1+delta_opzi/opzi)-Ok0 # *(1+zi)**3
 #   bH = rhsf(1/(1+zi), 0, Om)*
    return np.sqrt(Om/a*(1+delta_opzi/opzi)- Ok0 + OL0*a**2) 

def a_from_t_for_delta(delta, lrhs):
    a0 = [1e-7] # at a=1 we want da/dt=1 we can later get units from H0.
    zi = 1/a0[0]-1
    t = np.linspace(1e-6,5, 190000)

    ts=t
    sol = odeint(lrhs, a0, ts, args= (Om, zi, delta) )

    ts = ts[np.isfinite(sol.flatten())]
    sol = sol[np.isfinite(sol.flatten())].flatten()
    return ts, interp1d(ts,sol,bounds_error=None, fill_value="extrapolate")

def rhsf(a,t,Om):
    return np.sqrt(Om/a+(1-Om)*a**2)
def rhs(a,t, Om, OL0):
    return np.sqrt(Om/a-(1-Om-OL0) + OL0*a**2)

def rhs_sph(a,t, Om, zi, delta):
    Ok0 = (1+zi)*Om*delta *5/3
    OL0 = 1-Om-Ok0 
    return np.sqrt(Om/a*(1+delta)-Ok0 + OL0*a**2) #/bH

Om = .3
ts, a_from_t = a_from_t_for_delta(0., rhs_del)
d1 = 10
ts, ad1 = a_from_t_for_delta(d1, rhs_del)
ts, mad1 = a_from_t_for_delta(-d1, rhs_del)
#ts, as1 = a_from_t_for_delta(d1, rhs_sph)
#ts, mas1 = a_from_t_for_delta(-d1, rhs_sph)

plt.plot(1/a_from_t(ts), (mad1(ts)/a_from_t(ts))**3-1,'--',c=cycle[0],label=r"$a_\delta/a_b \delta/a_i=-{:.1g}$".format(d1),alpha=.5)
plt.plot(1/a_from_t(ts), -((ad1(ts)/a_from_t(ts))**3-1),c=cycle[0],label=r"$a_b/a_\delta\, \delta/a_i={:.1g}$".format(d1),alpha=.5)
#plt.plot(1/a_from_t(ts), mas1(ts)/a_from_t(ts)-1,'--',c=cycle[1],label=r"$\delta=-{:.1g}$".format(d1),alpha=.5)
#plt.plot(1/a_from_t(ts), as1(ts)/a_from_t(ts)-1,c=cycle[1],label=r"$\delta={:.1g}$".format(d1),alpha=.5)

plt.plot(1/a_from_t(ts),d1* a_from_t(ts),'-',c=cycle[2], alpha=.5, label=r"$\propto a$")

plt.xscale("log")
plt.yscale("log")
plt.xlim(100,1)
#plt.ylim(-1.3,1.3)
plt.ylim(1e-2,10)
plt.xlabel(r"$1+z$")
plt.ylabel(r"${a_\delta}/{a_b}$")
plt.title(r"change of scale factor vs background $\Omega_m={:.3f}$".format(Om))
#plt.loglog(1/solb,sol/solb,label=r"$\delta(z=99)=0.01$")
plt.legend()
# %%
def g_flat(ai, Omi):
    Ol = 1-Omi # flat Universe of only mater and lambda
    hubble = lambda a:  np.sqrt(Omi/a**3 + Ol)
    rhs_g = lambda a: 1/a**3/hubble(a)**3
    gint = quad(rhs_g,0,ai)[0]
    g_f = 2.5*Omi*hubble(ai)/ai
    return g_f*gint
# Fitting formula for flat Universes
#    Om = Omi/a**3
#    return (5*Om/(2*Om**(4/7)-Ol+(1+Om/2)*(1+Ol/70)))

g_flat_vec = np.vectorize(g_flat)

a = np.logspace(-2,0,50)
plt.figure(figsize=(10,8))
plt.plot(1/a, g_flat_vec(a,1)*a, label=r"$\Omega_m=1, \, \Omega_\Lambda=0$")
plt.plot(1/a, g_flat_vec(a,0.4)*a, label=r"$\Omega_m=0.4, \, \Omega_\Lambda=0.6$")
plt.plot(1/a, g_flat_vec(a,0.3)*a, label=r"$\Omega_m=0.3, \, \Omega_\Lambda=0.7$")
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
plt.show()
plt.savefig("../assets/growthfactor-inhomogeneous-as-function-of-a.png", dpi=150,transparent=figtransparent)



#%%

d1 = 0.5
d2 = 1
d3 = 7
ts, ad1 = a_from_t_for_delta(d1, rhs_del)
ts, mad1 = a_from_t_for_delta(-d1, rhs_del)
ts, ad2 = a_from_t_for_delta(d2, rhs_del)
ts, mad2 = a_from_t_for_delta(-d2, rhs_del)
ts, ad3 = a_from_t_for_delta(d3, rhs_del)
ts, mad3 = a_from_t_for_delta(-d3, rhs_del)

plt.plot(1/a_from_t(ts), (mad3(ts)/a_from_t(ts))**3-1,'--',c=cycle[2],label=r"$\delta/a_i=-{:.1g}$".format(d3))
plt.plot(1/a_from_t(ts), (mad2(ts)/a_from_t(ts))**3-1,'--',c=cycle[1],label=r"$\delta/a_i=-{:.1g}$".format(d2))
plt.plot(1/a_from_t(ts), (mad1(ts)/a_from_t(ts))**3-1,'--',c=cycle[0],label=r"$\delta/a_i=-{:.1g}$".format(d1))
plt.plot(1/a_from_t(ts), (ad1(ts)/a_from_t(ts))**3-1,c=cycle[0],label=r"$\delta/a_i={:.1g}$".format(d1))
plt.plot(1/a_from_t(ts), -1*((mad1(ts)/a_from_t(ts))**3-1),'--',c=cycle[0],alpha=.9,lw=1)

plt.plot(1/a_from_t(ts), (ad2(ts)/a_from_t(ts))**3-1,c=cycle[1],label=r"$\delta/a_i={:.1g}$".format(d2))
plt.plot(1/a_from_t(ts), -1*((mad2(ts)/a_from_t(ts))**3-1),'--',c=cycle[1],alpha=.9,lw=1)

plt.plot(1/a_from_t(ts), (ad3(ts)/a_from_t(ts))**3-1,c=cycle[2],label=r"$\delta/a_i={:.1g}$".format(d3))
plt.plot(1/a_from_t(ts), -1*((mad3(ts)/a_from_t(ts))**3-1),'--',c=cycle[2],alpha=.9,lw=1)

plt.plot(1/a_from_t(ts), -(d1*a_from_t(ts)),'-.',c=cycle[0],alpha=.9,lw=2,label=r"$\propto a$")
plt.plot(1/a_from_t(ts), (d1*a_from_t(ts)),'-.',c=cycle[0],alpha=.9,lw=2)
plt.plot(1/a_from_t(ts), -(d2*a_from_t(ts)),'-.',c=cycle[1],alpha=.9,lw=2)
plt.plot(1/a_from_t(ts), (d2*a_from_t(ts)),'-.',c=cycle[1],alpha=.9,lw=2)
plt.plot(1/a_from_t(ts), -(d3*a_from_t(ts)),'-.',c=cycle[2],alpha=.9,lw=2)
plt.plot(1/a_from_t(ts), (d3*a_from_t(ts)),'-.',c=cycle[2],alpha=.9,lw=2)


plt.xscale("log")
#plt.yscale("log")
plt.xlim(1000,1)
plt.ylim(-1.25,1.25)
#plt.ylim(1e-3,.3)
plt.xlabel(r"$1+z$")
plt.ylabel(r"${\delta=(a_\delta}/{a_b})^3-1$")
plt.title(r"change of density vs background $\Omega_m={:.1f}$".format(Om))
#plt.loglog(1/solb,sol/solb,label=r"$\delta(z=99)=0.01$")
plt.legend()

# %%

d1 = 0.5
d2 = 1
d3 = 7
ts, ad1 = a_from_t_for_delta(d1, rhs_del)
ts, mad1 = a_from_t_for_delta(-d1, rhs_del)
ts, ad2 = a_from_t_for_delta(d2, rhs_del)
ts, mad2 = a_from_t_for_delta(-d2, rhs_del)
ts, ad3 = a_from_t_for_delta(d3, rhs_del)
ts, mad3 = a_from_t_for_delta(-d3, rhs_del)

plt.plot(1/a_from_t(ts), ((mad3(ts)/a_from_t(ts))**3-1)/(d3*a_from_t(ts)),'--',c=cycle[2],label=r"$\delta/a_i=-{:.1g}$".format(d3))
plt.plot(1/a_from_t(ts), ((mad2(ts)/a_from_t(ts))**3-1)/(d2*a_from_t(ts)),'--',c=cycle[1],label=r"$\delta/a_i=-{:.1g}$".format(d2))
plt.plot(1/a_from_t(ts), ((mad1(ts)/a_from_t(ts))**3-1)/(d1*a_from_t(ts)),'--',c=cycle[0],label=r"$\delta/a_i=-{:.1g}$".format(d1))
plt.plot(1/a_from_t(ts), ((ad1(ts)/a_from_t(ts))**3-1)/(d1*a_from_t(ts)),c=cycle[0],label=r"$\delta/a_i={:.1g}$".format(d1))
plt.plot(1/a_from_t(ts), 1*((mad1(ts)/a_from_t(ts))**3-1)/(d1*a_from_t(ts)),'--',c=cycle[0],alpha=.9,lw=1)

plt.plot(1/a_from_t(ts), ((ad2(ts)/a_from_t(ts))**3-1)/(d2*a_from_t(ts)),c=cycle[1],label=r"$\delta/a_i={:.1g}$".format(d2))
plt.plot(1/a_from_t(ts), 1*((mad2(ts)/a_from_t(ts))**3-1)/(d2*a_from_t(ts)),'--',c=cycle[1],alpha=.9,lw=1)

plt.plot(1/a_from_t(ts), ((ad3(ts)/a_from_t(ts))**3+1)/(1+d3*a_from_t(ts)),c=cycle[2],label=r"$\delta/a_i={:.1g}$".format(d3))
plt.plot(1/a_from_t(ts), 1*((mad3(ts)/a_from_t(ts))**3+1)/(1+d3*a_from_t(ts)),'--',c=cycle[2],alpha=.9,lw=1)

#plt.plot(1/a_from_t(ts), -(d1*a_from_t(ts)),'-.',c=cycle[0],alpha=.9,lw=2,label=r"$\propto a$")
#plt.plot(1/a_from_t(ts), (d1*a_from_t(ts)),'-.',c=cycle[0],alpha=.9,lw=2)
#plt.plot(1/a_from_t(ts), -(d2*a_from_t(ts)),'-.',c=cycle[1],alpha=.9,lw=2)
#plt.plot(1/a_from_t(ts), (d2*a_from_t(ts)),'-.',c=cycle[1],alpha=.9,lw=2)
#plt.plot(1/a_from_t(ts), -(d3*a_from_t(ts)),'-.',c=cycle[2],alpha=.9,lw=2)
#plt.plot(1/a_from_t(ts), (d3*a_from_t(ts)),'-.',c=cycle[2],alpha=.9,lw=2)


plt.xscale("log")
#plt.yscale("log")
plt.xlim(1000,1)
plt.ylim(-1.25,2.25)
#plt.ylim(1e-3,.3)
plt.xlabel(r"$1+z$")
plt.ylabel(r"${\delta=(a_\delta}/{a_b})^3-1$")
plt.title(r"change of density vs background $\Omega_m={:.1f}$".format(Om))
#plt.loglog(1/solb,sol/solb,label=r"$\delta(z=99)=0.01$")
plt.legend()

#%%
# ## 
d1 = 5e-4
d2 = 1e-3
d3 = 7e-3
ts, ad1 = a_from_t_for_delta(d1, rhs_del)
ts, mad1 = a_from_t_for_delta(-d1, rhs_del)
ts, ad2 = a_from_t_for_delta(d2, rhs_del)
ts, mad2 = a_from_t_for_delta(-d2, rhs_del)
ts, ad3 = a_from_t_for_delta(d3, rhs_del)
ts, mad3 = a_from_t_for_delta(-d3, rhs_del)

plt.plot(1/a_from_t(ts), (mad3(ts)/a_from_t(ts)-1)/d3/1e3,'--',c=cycle[2],label=r"$\delta(z=999)=-{:.1g}$".format(d3))
plt.plot(1/a_from_t(ts), (mad2(ts)/a_from_t(ts)-1)/d2/1e3,'--',c=cycle[1],label=r"$\delta=-{:.1g}$".format(d2))
plt.plot(1/a_from_t(ts), (mad1(ts)/a_from_t(ts)/d1-1)/d1/1e3,'--',c=cycle[0],label=r"$\delta=-{:.1g}$".format(d1))
plt.plot(1/a_from_t(ts), (ad1(ts)/a_from_t(ts)-1)/d1/1e3,c=cycle[0],label=r"$\delta={:.1g}$".format(d1))
plt.plot(1/a_from_t(ts), -1*(mad1(ts)/a_from_t(ts)-1)/d1/1e3,'--',c=cycle[0],alpha=.9,lw=1)

plt.plot(1/a_from_t(ts), (ad2(ts)/a_from_t(ts)-1)/d2/1e3,c=cycle[1],label=r"$\delta={:.1g}$".format(d2))
plt.plot(1/a_from_t(ts), -1*(mad2(ts)/a_from_t(ts)-1)/d2/1e3,'--',c=cycle[1],alpha=.9,lw=1)

plt.plot(1/a_from_t(ts), (ad3(ts)/a_from_t(ts)-1)/d3/1e3,c=cycle[2],label=r"$\delta={:.1g}$".format(d3))
plt.plot(1/a_from_t(ts), -1*(mad3(ts)/a_from_t(ts)-1)/d3/1e3,'--',c=cycle[2],alpha=.9,lw=1)

da = ((1+d1*a_from_t(ts)/a_from_t(0))**(-1/3)-1)/d1/1e3
plt.plot(1/a_from_t(ts),da*a_from_t(ts), "o", c="grey", markevery=10, ms=6, label=r"growth factor $g(a)$")


plt.xscale("log")
plt.xlim(1000,1)
plt.ylim(-.5,.5)
plt.xlabel(r"$1+z$")
plt.ylabel(r"$({a_\delta}/{a_b}-1)/\delta$")
plt.title(r"change of scale factor vs background $\Omega_m={:.1f}$".format(Om))#plt.loglog(1/solb,sol/solb,label=r"$\delta(z=99)=0.01$")
plt.legend()

# %%


# %%

delt = np.logspace(-5,-2,20)
d1 = 5e-4
d2 = 1e-3
d3 = 7e-3
ts, ad1 = a_from_t_for_delta(d1, rhs_del)
ts, mad1 = a_from_t_for_delta(-d1, rhs_del)
ts, ad2 = a_from_t_for_delta(d2, rhs_del)
ts, mad2 = a_from_t_for_delta(-d2, rhs_del)
ts, ad3 = a_from_t_for_delta(d3, rhs_del)
ts, mad3 = a_from_t_for_delta(-d3, rhs_del)

plt.plot(1/a_from_t(ts), mad3(ts),'--',c=cycle[2],label=r"$\delta(z=999)=-{:.1g}$".format(d3))
plt.plot(1/a_from_t(ts), ad3(ts),'--',c=cycle[2],label=r"$\delta(z=999)={:.1g}$".format(d3))

plt.plot(1/a_from_t(ts), mad2(ts),'--',c=cycle[1],label=r"$\delta(z=999)=-{:.1g}$".format(d2))
plt.plot(1/a_from_t(ts), ad2(ts),'--',c=cycle[1],label=r"$\delta(z=999)={:.1g}$".format(d2))

plt.plot(1/a_from_t(ts), mad2(ts),'--',c=cycle[1],label=r"$\delta(z=999)=-{:.1g}$".format(d2))
plt.plot(1/a_from_t(ts), ad2(ts),'--',c=cycle[1],label=r"$\delta(z=999)={:.1g}$".format(d2))

plt.xscale("log")
plt.yscale("log")
plt.xlim(1000,1)
plt.ylim(1e-3,2)
plt.xlabel(r"$1+z$")
plt.ylabel(r"$({a_\delta}/{a_b}-1)/\delta$")
plt.title(r"change of scale factor vs background $\Omega_m={:.1f}$".format(Om))#plt.loglog(1/solb,sol/solb,label=r"$\delta(z=99)=0.01$")
plt.legend()
# %%
## Growth Factor comparison


a0 = [1e-3] # at a=1 we want da/dt=1 we can later get units from H0.
zi = 999
t = np.linspace(1e-6,5, 190000)
ts=t
sol = odeint(rhsf, a0, ts, args= (Om, ) )
ts = ts[np.isfinite(sol.flatten())]
sol = sol[np.isfinite(sol.flatten())].flatten()
t_from_a = interp1d(sol,ts,bounds_error=None,fill_value="extrapolate")


ctime = t_from_a(1)
ft = np.array([1e-6,ctime])
delt = -1*np.logspace(-4,1.6,100)


fdel = np.array([])
for delta in delt:
    sol = odeint(rhs_del, 1e-3, ft, args= (Om, 999, 1e-3*delta))
    fdel = np.append(fdel,sol[1])

#plt.plot(-delt,-(fdel**3-1)/delt,'o')
#plt.plot(-delt,-(fsphdel**3-1)/delt,'o')
plt.plot(-delt,(fdel-1),'o')
#plt.plot(-delt,fsphdel-1,'o')
plt.xlabel(r"$\delta/a_i$")
plt.xscale("log")
plt.ylim(1e-4,1e2)
plt.yscale("log")
# %%

# %%
import veusz.embed as veusz

g = veusz.Embedded('new win')
g.EnableToolbar(enable=True)
g.To( g.Add('page') )
g.To( g.Add('graph') )
g.SetData('x', np.arange(20))
g.SetData('y', np.arange(20)**2)
g.Add('xy')
g.Zoom(0.5)

g.WaitForClose()# %%

# %%

