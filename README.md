Vascular and interstitial flow solver for discrete microvascular networks

A C++ programme using the Armadillo library (Sanderson & Curtin, 2016), which loads discrete microvascular networks and solves for vascular blood flow (blood pressure, vessel wall shear stress etc) and interstitial fluid transport (interstitial fluid pressure and velocity fields). This software package has been applied to realistic, whole vascular tumours in d'Esposito et al. 2018. Code was written by Paul W. Sweeney (with biconjugate gradient method and contour files used and altered from Timothy Secomb’s O2 transport code: https://physiology.arizona.edu/people/secomb/greens), under the supervision of Rebecca J. Shipley and Simon Walker-Samuel based at University College London in Mechanical Engineering and the Centre for Advanced Biomedical Imaging.

The vascular component uses Poiseuille's Law for flow, empirical laws to calculate blood viscosity (Pries et al. 2005) and (if selected) uses an iterative scheme to solve for vascular haematocrit distributions (Pries et al. 1989).

The interstitial component solves Darcy's Law in an axisymmetric, spherical domain where point sources of flux are distributed along the vasculature and Green's functions to link the vascular and interstitial domains.

References:

d'Esposito A*, Sweeney P*, Ali M, Saleh M, Ramasawmy R, Roberts T, Agliardi G, Desjardins A, Lythgoe M, Pedley R, Shipley R*, Walker-Samuel S*. (2018). Computational fluid dynamics with imaging of cleared tissue and of in vivo perfusion predicts drug uptake and treatment responses in tumours. Nature Biomedical Engineering.
*Joint first/last author.

C. Sanderson and R. Curtin. Armadillo: a template-based C++ library for linear algebra. Journal of Open Source Software, Vol. 1, pp. 26, 2016.

Pries, A. R. and Secomb, T. W. 2005. Microvascular blood viscosity in vivo and the endothelial surface layer. American journal of physiology. Heart and circulatory physiology 289:H2657–H2664.

Pries, A. R., Ley, K., Claassen, M., and Gaehtgens, P. 1989. Red cell distribution at microvascular bifurcations. Microvascular research 38:81–101.
