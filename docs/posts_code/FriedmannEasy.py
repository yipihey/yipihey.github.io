# %% 
from definitions import *
h = 0.7
Mpc = 3.086e24
H100 = 100*1e5/Mpc
H0 = H100*h 
G = 6.6743e-8
c = 2.99792458e10
rhocrit = 3*H0**2/(8*np.pi*G)
k = 0
R_c = 10e3 * Mpc  # 10 Gpc

class cosmo:
    """Keep parameters for a flat Universe
    H0 = Hubble constant today: in km/s/Mpc
    Om0: Omega_matter at current time
    OL0: Omega_lambda will be set to 1-Omega_matter for a flat Universe
     """
    def __init__(self,params):
        self.H0 = params['H0']
        self.Om0 = params['Om0']
        self.OL0 = 1.0 - self.Om0
        
#        print ('H0', self.H0)
#        print ('Om0', self.Om0)
#        print ('OL0', self.OL0)
        
    def hubble(self,a):
        return np.sqrt(self.H0**2*(self.Om0/a**3+(1-self.Om0)))

    def hubble_from_z(self,z):
        return hubble(self,1/(1+z))

params = {'H0':67.11, 'Om0':0.3175}

u1 = cosmo( {'H0':67.11, 'Om0':0.3175})
u2 = cosmo( {'H0':74, 'Om0':0.3175})
u3 = cosmo( {'H0':67.11, 'Om0':0.1})
u4 = cosmo( {'H0':67.11, 'Om0':0.5})
universe = cosmo(params)

ar = np.logspace(-2,.4,400)
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
plt.savefig("../assets/Hubble-as-function-of-a.png", dpi=150,transparent=figtransparent)
# %% [markdown]
# <h2>Integrate</h2>
#

# %%
from scipy.integrate import odeint

def rhs(a,t,Om):
    return np.sqrt(Om/a+(1-Om)*a**2)
a0 = [1] # at a=1 we want da/dt=1 we can later get units from H0.
Om = 1
t = np.linspace(0, 1.5, 19001)
tm = -1*t
sol = odeint(rhs, a0, t, args= (Om, ) )
solb = odeint(rhs, a0, tm, args= (Om, ) )

plt.figure(figsize=(10,6))

p1 = plt.plot(t, (odeint(rhs, a0, t , args=(u1.Om0, )))[:,0], label=r"$\Omega_m={:.2f}$".format(u1.Om0), alpha=.7)
plt.plot(tm,    (odeint(rhs, a0, tm, args=(u1.Om0, )))[:,0], alpha=.7, color=p1[0].get_color())

p3 = plt.plot(t, (odeint(rhs, a0, t , args=(u3.Om0, )))[:,0], label=r"$\Omega_m={:.2f}$".format(u3.Om0), alpha=.7)
plt.plot(tm,    (odeint(rhs, a0, tm, args=(u3.Om0, )))[:,0], alpha=.7, color=p3[0].get_color())

p4 = plt.plot(t, (odeint(rhs, a0, t , args=(u4.Om0, )))[:,0], label=r"$\Omega_m={:.2f}$".format(u4.Om0), alpha=.7)
plt.plot(tm,    (odeint(rhs, a0, tm, args=(u4.Om0, )))[:,0], alpha=.7, color=p4[0].get_color())

p2 = plt.plot(t, (odeint(rhs, a0, t , args=(1, )))[:,0], label=r"$\Omega_m=1$", alpha=.7)
plt.plot(tm,    (odeint(rhs, a0, tm, args=(1, )))[:,0], alpha=.7, color=p2[0].get_color())


plt.plot([-2/3,-2/3],[1e-5,1e4],'--',lw=2,alpha=.5 ,color="grey",label=r"$\Omega_m=1$")
plt.plot([0],[1],'o',ms=8,alpha=.75, color="grey", label="today")


#plt.xscale("log")
plt.yscale("log")
plt.xlim(-1.3,.49)
plt.ylim(1e-2,2)
plt.xlabel(r"$(t-t_{today}) H_0$")
plt.ylabel(r"$a(t)$")
plt.legend();
plt.title("spatially flat cosmologies")

plt.savefig("../assets/flat-a-of-t.png", dpi=150,transparent=figtransparent)
# %%
# Zoomes in version of the same plot from before
plt.figure(figsize=(10,6))

p1 = plt.plot(t, (odeint(rhs, a0, t , args=(u1.Om0, )))[:,0], label=r"$\Omega_m={:.2f}$".format(u1.Om0), alpha=.7)
plt.plot(tm,    (odeint(rhs, a0, tm, args=(u1.Om0, )))[:,0], alpha=.7, color=p1[0].get_color())

p3 = plt.plot(t, (odeint(rhs, a0, t , args=(u3.Om0, )))[:,0], label=r"$\Omega_m={:.2f}$".format(u3.Om0), alpha=.7)
plt.plot(tm,    (odeint(rhs, a0, tm, args=(u3.Om0, )))[:,0], alpha=.7, color=p3[0].get_color())

p4 = plt.plot(t, (odeint(rhs, a0, t , args=(u4.Om0, )))[:,0], label=r"$\Omega_m={:.2f}$".format(u4.Om0), alpha=.7)
plt.plot(tm,    (odeint(rhs, a0, tm, args=(u4.Om0, )))[:,0], alpha=.7, color=p4[0].get_color())

p2 = plt.plot(t, (odeint(rhs, a0, t , args=(1, )))[:,0], label=r"$\Omega_m=1$", alpha=.7)
plt.plot(tm,    (odeint(rhs, a0, tm, args=(1, )))[:,0], alpha=.7, color=p2[0].get_color())


plt.plot([-2/3,-2/3],[1e-5,1e4],'--',lw=2,alpha=.5 ,color="grey",label=r"$\Omega_m=1$")
plt.plot([0],[1],'o',ms=8,alpha=.75, color="grey", label="today")


#plt.xscale("log")
#plt.yscale("log")
plt.xlim(-.66,.1)
plt.ylim(.3,1.1)
plt.xlabel(r"$(t-t_{today}) H_0$")
plt.ylabel(r"$a(t)$")
plt.legend();
plt.title("spatially flat cosmologies")

plt.savefig("../assets/flat-a-of-t-zoom.png", dpi=150,transparent=figtransparent)
# %%

