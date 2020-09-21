---
mathjax: true
layout: post
title:  "The Inhomogeneous Friedmann Universe"
date:   2020-09-10
categories: fragments
author: Tom
---

## The Inhomogeneous Friedmann Universe

To get insight on how structure like galaxies and the cosmic web they inhabit format we continue the previous discussion of [the Friedmann Universe]({% post_url 2020-09-10-Friedmann %}) and see what happens to initial inhomogeneities in this expanding background.

We saw that the Friedmann equation is:

$$\frac{\dot{a}^2 + kc^2/R_c^2}{a^2} = \frac{8 \pi G \rho + \Lambda c^2}{3}, $$

and that using $\rho_{crit,0}= {3H_0^2 \over 8\pi G}$ and $H=\dot{a}/a$ we could rewrite this as

$$H(a)^2= \Omega_{m,0}\,a^{-3} - \Omega_{k,0}\,a^{-2} + \Omega_\Lambda,
$$

where $\Omega_k= k c^2/R_c^2/\rho_{crit,0}$ and $\Omega_\Lambda = \Lambda c^2/3/\rho_{crit,0}$. One important thing to note is that all the $\Omega$s have to add up to 1. This is true at all times but we wrote the equation here so that we give all the $\Omega$s only at the time when $a=1$ and keep their (general versions) $a$ dependence explicit in the equation. This is indicated by using the subscript $0$ in $\Omega_{m,0}$. The zero refers to the cosmological redshift $1+z\equiv 1/a$ which is often used as an  way to denote which epoch of the Universe one is referring to. For a fixed cosmological model either of $a$, $z$, or time are used for this depending on what one is interesting to highlight. How much the light redshifts from that period to now, how old the Universe was or how much smaller a specific distance would have been at that time.

All $\Omega$s have to add to 1 by definition. In an analogy to the Newtonian equivalent equation the Friedmann equation is derived from the Energy condition. The $k$ denotes there whether the total energy is exactly zero, negative for a gravitational bound objects or positive for which the total energy would be positive.


At the current time when $a=1$ and $\dot{a}=H_0$, $\Omega_{k,0} = 1 - \Omega_{m,0} - \Omega_{\Lambda,0}$, can be thought of as the amount of matter missing that would have made the Universe closed. If $\Omega_{k,0}<0$ then we have a Universe that is closed and will collapse eventually. 

If we consider these density parameters at an earlier time, say at $z=99$ where $a=0.01$ they would take very different numerical values because the $\rho_{crit}$ which they are measured against depends on $H(z)^2$ 

Consider a Universe with matter density $\rho$ that has a matter density slightly different than a flat universe with the same cosmological constant and matter density $\rho_b$. Write this difference as a density contrast $\delta=\rho/\rho_b-1$.
This Universe now has a curvature that depends on $\delta$. So we would have 

$$\frac{\dot{a}^2}{a^2} = \frac{8 \pi G \rho}{3} (1+\delta) -  \frac{kc^2}{a^2R_c^2} + \frac{\Lambda c^2}{3}, $$

and if we assume that for a very small delta at very high redshift $z_i$ (small $a$) it expands as fast as the background Universe for which,

$$
\frac{\dot{a}^2}{a^2} = \frac{8 \pi G \rho}{3} + \frac{\Lambda c^2}{3}, 
$$

we can set these two equations equal and see that,

$$\frac{8 \pi G \rho}{3} \delta =  \frac{kc^2}{a^2R_c^2}. $$ 

In terms of density parameters that tells us that the curvature corresponding to a slight density perturbation $\delta$ at redshift $z_i$ is given by:

$$\Omega_{k,0} = (1+z_i)\Omega_{m,0}\delta .$$

To derive this we assumed that the expansion rate in our perturbed Universe is exactly the same as the one of the unperturbed. If one includes that the Hubble rate at $z_i$ where we specify $\delta$ is perturbed so that $H(z_i)=H_b(z_i)(1-\delta)$ one finds the equivalent $\Omega_{k,0}$. 

$$\Omega_{k,0} = {\delta \left[ a_i^3(\delta-6)(\Omega_{m,0}-1)- (\delta-15)\Omega_{m,0} \right] \over 9 a_i}
  \approx \frac{5}{3}(1+z_i)\Omega_{m,0}\delta .$$

The last approximation here is very accurate as all other terms have at least order $a_i^2$ or $\delta^2$ or even smaller. 
So we can now compare now how different the scale factors evolve in the flat background Universe as compared to the one corresponding to the one with the slight early difference in density. 

We can also define the parameter describing the evolution as $\delta_0\equiv \delta(z)/a= \delta(z) (1+z)$ so that we can write:

$$\Omega_{k,0} = {\delta_0 \left[ a_i^3(\delta_0/a_i-6)(\Omega_{m,0}-1)- (\delta_0/a_i-15)\Omega_{m,0} \right] \over 9}
  \approx \frac{5}{3}\Omega_{m,0}\delta_0.$$


