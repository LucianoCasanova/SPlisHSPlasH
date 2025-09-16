0_0_X -> van der Waals nonPressureForce + IISPH

0_0_0   r: 0.05385  rho0: 1000.61   a: 2.0  R: 1.5
0_0_1   r: 0.05385  rho0: 1000.61   a: 3.0  R: 1.5
0_0_2   r: 0.05385  rho0: 1000.61   a: 1.0  R: 1.5
0_0_3   r: 0.05385  rho0: 1000.61   a: 2.0  R: 1.0
0_0_4   r: 0.05385  rho0: 1000.61   a: 2.0  R: 2.0

m = 0.8 * 8 * r^3 * rho0
a_bar = a/m^2
To copy paper, I wanted m=1 with rho0=1000 -> r=0.05385

To study the effect of discretization, what should I change?
|-> rho0 should stay the same.
|-> r will change.
|-> Then m changes.
|-> a_bar scales with 1/m^2 so the relative strength stays the same. 

0_0_5   r: 0.04274  rho0: 1000.61   a: 2.0  R: 1.5 (m: 0.5)
0_0_6   r: 0.06164  rho0: 1000.61   a: 2.0  R: 1.5 (m: 1.5)

Now I want to study the influence of H (before it was 2h).
|-> Changed in SPHKernels.h and SurfaceTension_vdW.cpp

0_0_7   r: 0.05385  rho0: 1000.61   a: 2.0  R: 1.5  H: 3h
0_0_8   r: 0.05385  rho0: 1000.61   a: 2.0  R: 1.5  H: h    (breaks)
0_0_9   r: 0.05385  rho0: 1000.61   a: 2.0  R: 1.5  H: 4h

No habia static_cast<Real> en 7 y 0... preguntar a gpt

0_1_X -> vdWSPH (TimeStep) + van der Waals nonPressureForce

Not working