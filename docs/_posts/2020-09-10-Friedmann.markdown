---
mathjax: true
layout: post
title:  "The Friedmann Universe"
date:   2020-09-10
categories: fragments
author: Tom
---

# The Friedmann Universe

Assuming an isotropic and homogneous Universe allowed Friedmann to drastically simplify Einsteins equations and get quite straigh-forward equations that describe how distances and time evolve for different assumptions about the Universe. E.g. the more mass you have in that Universe will make a big difference in how quickly an initial expansion slows down just as your intuition would expect balls thrown on heavier planets landing sooner than what we are used to.
 
## Introduction to the Friedmann Equation
The [Friedmann equation](https://en.wikipedia.org/wiki/Friedmann_equations) describes the evolution of a homogeneous isotropic Universe. Distances change according to a scale factor $a(t)$. We use the typical convention in which $a(t=today)=1$ so that any scale $R_i$ at the present day can be thought to actually have the size $a(t)R_i$ at any other time. $R_i$ is typically called the co-moving scale or radius. So in an assumed homogeneous and isotropic Universe all distances scale by the same factor. Whether you think of a 1cm cherry or a $1.2\times 10^{11}$cm star. In fact one is usually concerned with distance measured in Mega [Parsec](https://en.wikipedia.org/wiki/Parsec) ($1\mathrm{Mpc}\approx 3.086\times 10^{24}$cm). So we see $a(t)$ is unit-less. It gives just the ratio of length-scales at different times. 

One way of writing the Friedmann equation that describes how $a$ changes over time involves its derivative with respect to time, $\dot{a}\equiv {da \over dt}$. It is given by this combination,

$$\frac{\dot{a}^2 + kc^2/R_c^2}{a^2} = \frac{8 \pi G \rho + \Lambda c^2}{3}. $$

It involves a few more constants. $G=6.6743 \times 10^{-8} \mathrm{cm^3 g^{-1} s^{-2}}$ is Newton's gravitational constant, $c=2.99792458 \times 10^{10} {\mathrm{ cm\, s}}^{-1}$ the speed of light in vacuum. $\Lambda$ is Einstein's famous cosmological constant that is now thought to make an enormous difference to the future of the Universe. Last but not least, $k$, is a constant that takes one of the three possible values -1, 0, 1 and describes the overall geometry of the Universe. $k/a^2$ describes the curvature of the Universe and for any Universe $k$ remains fixed. The parameter $R_c$ is a length scale specifying the curvature which is important when $k\neq 0$. The remaining symbol not yet explained is $\rho$. This is the mass energy density in units of $\mathrm{g/cm^3}$. Assuming a spherical blob of mass, $M_i$ in a sphere with radius, $R_i$, at the present time that is co-moving with the Universe. Its density would then be $\rho(t)=M_i/[4\pi/3\, (a\,R_i)^3 ]$. This scaling with $a^{-3}$ is correct for most normal matter in which most of the energy is in form of its rest mass. For radiation or relativistically moving particles one obtains an extra $1/a$ factor accounting for cosmological redshifting of the energy on top of the simple volume dilution. For the following the use of $\rho\propto a^{-3}$ is an excellent approximation for all times since $a$ was about $0.01$. 

Now this is a first order differential equation and to integrate it will require a suitable boundary condition to fix the associated integration constant. If the age of the Universe were given, together with the choice of $a=1$ today, that would be all that is needed. However, this equation is used to infer the age of the Universe from measurements astronomers make! Alternatively, one can specify the timescale of expansion using the ratio $a/\dot{a}$ which has units of time. For historical reasons its inverse is what is usually talked about, the [Hubble constant](https://en.wikipedia.org/wiki/Hubble%27s_law), $H(t)\equiv \dot{a}/a$, most often quoted with units of $\mathrm{km/s/Mpc}$. So for something $300 \mathrm{Mpc}$ away from us we expect $300 \times H_0$ $\mathrm{km/s}$ large speed of expansion. So the Hubble constant relates the speed of expansion to distance. A word of caution here. The Hubble constant as defined here ($\dot{a}/{a}$) is in fact not constant. It evolves by the Friedmann equation. However, its value today which we denote $H_0\approx 70 \mathrm{km/s/Mpc}\approx 2.27\times 10^{-18}\mathrm{s^{-1}}$ is a constant and it will serve as boundary condition to integrate the Friedmann equation. 

So specifying size and speed ($a$, $\dot{a}$) at the current time allows then to integrate the Friedmann equation both backwards and forwards in time. 

## Density Parameters
One more thing of note of jargon you will see everywhere, are many types of $\Omega$'s. $\Omega_\Lambda$, $\Omega_B$, $\Omega_m$, $\Omega_0$, $\Omega_{\nu}$, $\Omega_{BH}$, $\Omega_R$ and so on. To see what these are about just see that *at the present time* $\dot{a}/a=H_0$ and when one sets $k=0$ and $\Lambda=0$ in the Friedmann equation it simplifies to $H^2=\dot{a}^2/a^2 = H_0^2=8\pi \rho G/3$. So it defines a particular density value that makes both sides equal which is called the critical density. 

$$\rho_{crit}=\frac{3H_0^2}{8\pi G}=1.878\times 10^{-29} h^2 \mathrm{g/cm^3}\approx 9.2\times 10^{-30}\mathrm{g/cm^3}$$

This then is used as a unit to describe average densities in the Universe of certain components (baryons, neutrinos, photons, black holes, stars, planets, hydrogen, etc.). For baryons, for example, these are all the particles like protons, neutrons etc. that make up all the familiar things like us, our planet and all the stars. Their average density is $\Omega_B h^2\sim 2.2\times 10^{-2}$. Little $h$ here is the value of the Hubble constant expressed in $100\mathrm{km/s/Mpc}$[^1]. The $\Omega$ means it is expressed in the units of the critical density. I.e. $\Omega_B\equiv \rho_B/\rho_{crit}$. So $\rho_B=\Omega_B \rho_{crit} \approx \Omega_B 1.878\times 10^{-29} h^2 \mathrm{g/cm^3}\approx 4.13\times 10^{-31}{g/cm^3}$ which is equivalent to about one proton every 4 cubic meters(!). So yes space is quite empty if you consider taking everything in the Universe and spreading it out uniformly. Remarkably the Big Bang theory does indeed suggest that the entire observable Universe existed for a very long time as an almost perfectly uniform expanding region before it eventually started to allow for more variety. 

## Interpreting the Equation


[^1]: so $h=0.7$ for a Hubble constant of $70\mathrm{km/s/Mpc}$