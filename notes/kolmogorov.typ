#set heading(numbering:"9.1")
#set math.equation(
  block:true,
  numbering: (..nums) => {
  let section = counter(heading).get().first()
  numbering("(1.1)", 9, ..nums)
})
#set document(title:[Chapter 9: Propagation through Atmospheric Turbulence])

#title()
= Kolmogorov's Theory of Turbulence
Notes based on J.Schmidt's Numerical Simulation of Optical Propagation.
#v(1em)

Turbulence, in the sense of optical propagation in the atmosphere, is caused by a variation of refractive index. Kolmogorov developed a theory with turbulent flow in which, the kinetic energy is transferred from large eddies to smaller eddies.

$L_0$ = Outer scale (Average size of *largest* eddies) $[m]$

$l_0$ = Inner scale (Average size of *smallest* eddies) $[m]$

The range of eddie sizes is called the *inertial subrange*. Kolmogorov states that eddies within this range are statistically *homogeneous* (spacially invariant) and *isotropic* (directionally invariant) within small regions of space. Properties like eddie velocity $v$ and refractive index $n$ have stationary increments. Hence the use of structure function rather than covariance.

Note: Structure function is a function of *distance* between 2 points rather than the properties of absolute positions.

The relationship between the average eddie velocity $v$ and the average scale size of eddies $r$.

#math.equation[$v prop r^(1/3)$]
$v$ = Average velocity of eddie

$r$ = Average size of eddie

Structure function is defined as the square of the speeds so the structure function of eddie velocity is defined as:

#math.equation[$D_v (r) = C_v^2r^(2/3)$]
$C_v^2$ = Velocity structure parameter

(the square is a result of the structure function definition)

For the rest of the document, $D$ denotes a structure function and the subscript represents the kind of structure function.

For laminar flow (occurs at small scales: $0 <= r << l_0$), the properties differ so the velocity structure function is defined as:

#math.equation[$D_v (r) = C_v^2l_0^(-4/3)r^2$]

For larger scales, the flow is highly anisotropic. The potential temperature $theta$ can follow a similar analysis since potential temperature and ordinary temperature are linearly related. The potential temperature structure function follows:

#math.equation[$D_(theta)(r)= cases(
  C_(theta)^2l_0^(-4/3)r^2 \,& 0 <= r << l_0,
  C_(theta)^2r^(2/3) \,& l_0<< r << L_0,
)$]<pottempstrfn>

$C_(theta)^2$ = Potential temperature structure parameter

To be able to produce a refractive index statistical model we need a few more considerations:

The refractive index $n$ at a point in space $bold(r)$ can be written as:
#math.equation[$n(bold(r))=mu_n (bold(r))+n_1(bold(r))$]<refidxpoint>
where, $mu_n approx 1$ and is slowly varying mean value of refractive index and 
$n_1(bold(r))$ is the deviation from the mean value. This means that $n_1(bold(r))$ is modelled with a 0 mean random process.

For *optical wavelengths*, the refractive index at a point in space can approximated by:

#math.equation[$n(bold(r)) approx 1 + 77.6 times 10^(-5)(P(bold(r)))/(T(bold(r))) "for" lambda = 0.5mu m$]<refidxoptwvl>
Where:

$P$ = Pressure [millibars]

$T$ = Temperatrue [Kelvin]

Variation in refractive index is given as:

#math.equation[$d n = 7.99 times 10^(-5)(d P-(d T)/(T^2))$]

Given that pressure is relatively *uniform* and potential temperature $theta$ is linearly related to ordinary temperature $T$, refractive index variation can then be written as:

#math.equation[$d n = 7.99 times 10^(-5)(d theta)/(T^2)$]

Variation in refractive index $d n$ is directly proportional to the variation in potential temperature $d theta$ so the refractive index structure function $D_n(r)$ follows the same power law as $D_(theta)(r)$:

#math.equation[$D_n (r) = cases(
  C_n^2l_0^(-4/3)r^2 \, & 0<=r<<l_0,
  C_n^2r^(-2/3) \, & l_0<<r<<L_0,
)$]<refidxstrfn>

$C_n^2$ = Refractive index structure parameter

- $C_n^2$ is measured in $m^(-2/3)$
- $C_n^2$ ranges from $10^(-17)$ to $10^(-13)m^(-2/3)$ 
- Small values at high altitude
- Large values near ground

It is often necessary to have a spectral description of the refractive index fluctations. The refractive index power spectral density $bold(Phi)_n (kappa)$ can be computed from @refidxstrfn. The Kolmogorov refractive index power spectral density is defined as:

#math.equation[$bold(Phi)_n^K (kappa) = 0.033C_n^2 kappa^(-11/3)$ for $1/L_0<<kappa<<1/l_0$]<PSDkolmogorov>

Where $kappa = 2pi (f_x bold(hat(i))+f_y bold(hat(j)))$ is the angular spacial frequency in rad/m.
@PSDkolmogorov is valid for random fields that are:
- Locally homogeneous
- Locally isotropic

There are other models such as the *von Karman* power spectral density (PSD):

#math.equation[$Phi_n^(v K) = (0.033C_n^2)/((kappa^2+kappa_0^2)^(11/6))$ for $0<=kappa<=1/l_0$]

And *modified von Karman* PSD:

#math.equation[$bold(Phi)_n^(m v K) = 0.033C_n^2(exp(-kappa^2"/"kappa_m^2))/((kappa^2+kappa_0^2)^(11/6))$ for $0<=kappa<infinity$]
$kappa_m = 5.92"/"l_0$

$kappa_0 = 2pi"/"L_0$

When dealing with electromagnetic wave propagation through the atmosphere, the refractive index can be considered as time independent because the speed of light is so fast, the times it takes for light to propagate through a turbulent eddie is much lower than the time it takes for the turbulent properties to change. Given the knowledge of wind speed, you can convert spatial statistics to temporal one.

#math.equation[$phi.alt(x,y,t) = phi.alt(x-v_x t,y-v_y t,0)$]
$v_x$ and $v_y$ are x,y components of the wind velocity

= Optical Propagation through Turbulence
#v(1em)
The wave equation that covers any of the six field components is:
#math.equation[$[nabla^2+k^2n^2(bold(r))]U(bold(r))=0$]<waveeqopt>

To solve @waveeqopt, we assume weak fluctations which means $|n_1(bold(r))|<<1$. Then to rewrite $n^2(bold(r))$, we can use @refidxpoint and the assumption to derive an approximation for the square of the refractive index:

#math.equation[$n^2(bold(r)) approx 1 + 2n_1(bold(r))$]<n2r>

Substituting @n2r into @waveeqopt, results in the wave equation shown below:

#math.equation[${nabla^2+k^2[1+2n_1(bold(r))]}U(bold(r)) = 0$]

When a medium has a constant refractive index, @waveeqopt is solved with Fourier Optics and Green's functions. When medium is randomly inhomogeneous, perturbative method with Green's function is used.

In Rytov's method, optical field is written as:

#math.equation[$U(bold(r)) = U_0(bold(r))exp(psi (bold(r)))$]

Where $U_0(bold(r))$ is the vacuum solution when $n_1 = 0$, and $psi(bold(r))$ is the complex phase perturbation of the form:

#math.equation[$psi(bold(r)) = psi_1(bold(r)) + psi(bold(r)) + ...$]
And 
#math.equation[$psi(bold(r)) = chi + i phi.alt$]
Where $chi$ is log-amplitude perturbation, and $phi.alt$ is phase perturbuation.

= Optical parameters of atmosphere
#v(1em)
The derivation of these parameters have been ommited, akin to the original book. Nevertheless, here are some
useful information from Rytov's theory:

- Mean value of optical field
#math.equation[$chevron.l U(bold(r) chevron.r = U_0(bold(r)) chevron.l exp psi(bold(r)) chevron.r$]

- Mutual coherence function #footnote[Not sure what $bold(r)'$ is referring to here]
#math.equation[$Gamma(bold(r),bold(r'),z) = chevron.l U(bold(r))U^*(bold(r)') chevron.r = U_0(bold(r))U^*_0(bold(r)') chevron.l exp[psi(bold(r))
psi^*(bold(r)')]chevron.r$]

- Modulus of complex coherence factor (coherence factor)
#math.equation[$mu(bold(r),bold(r)',z) = (|(Gamma(bold(r),bold(r)',z))|)/[Gamma(bold(r),bold(r)',z)Gamma(bold(r)',bold(r)',z)]^(1/2)$]

- Wave structure function #footnote[log-amplitude perturbation $chi$ structure function plus phase perturbation $phi.alt$ structure function]
#math.equation[$D(bold(r),bold(r)',z) = -2ln mu(bold(r),bold(r)',z) = D_chi (bold(r),bold(r)',z) + D_phi.alt (bold(r),bold(r)',z)$]

- Phase power spectral density
#math.equation[$Phi_phi.alt(kappa) = 1/(4pi^2kappa^2) integral_0^infinity (sin(kappa r))/(kappa r)d/(d r)[r^2d/(d r)D_phi.alt (r)]d r$]

- Mean MTF#footnote[Not sure what MTF is] of the turbulent path
#math.equation[$cal(H)(f) = exp[-1/2D(lambda f_l f)]$]
Where $f_l$ is the system focal length