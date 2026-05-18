#set math.equation(block:true,numbering: "(9.1)")
#set heading(numbering:"9.1")
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
$C_v$ = Velocity structure parameter

For the rest of the document, $D$ denotes a structure function and the subscript represents the kind of structure function.

For laminar flow (occurs at small scales: $0 <= r << l_0$), the properties differ so the velocity structure function is defined as:

#math.equation[$D_v (r) = C_v^2l_0^(-4/3)r^2$]

For larger scales, the flow is highly anisotropic. The potential temperature $theta$ can follow a similar analysis since potential temperature and ordinary temperature are linearly related.

#math.equation[$D_(phi)$$]