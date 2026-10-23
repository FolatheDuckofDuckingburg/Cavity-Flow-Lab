# Vorticity–Strain Alignment: Novelty Audit and Derivation Log

**Status:** exploratory research note; not peer reviewed.  
**Branch:** `research/alignment-novelty-audit`  
**Purpose:** record exact identities, assumptions, failed novelty claims, and concrete proof obligations for the Navier–Stokes manuscript.

> **Important:** This note does not claim a new regularity criterion or a solution of the 3D Navier–Stokes global-regularity problem. The identities below are derivations to check and compare against the literature. Novelty remains unestablished.

## 1. Equations and conventions

Consider a smooth incompressible solution on the three-dimensional periodic torus:
[
\partial_t u+(u\cdot\nabla)u=-\nabla p+\nu\Delta u,
\qquad \nabla\cdot u=0.
]
Let (D_t=\partial_t+u\cdot\nabla), (A=\nabla u), (S=(A+A^\top)/2), (W=(A-A^\top)/2), and (\omega=\nabla\times u). Pressure is normalized so that (\operatorname{tr} S=0) and (\Delta p=-\operatorname{tr}(A^2)).

The vorticity equation is
[
D_t\omega=S\omega+\nu\Delta\omega.
]
Where (\omega\ne0), set (\xi=\omega/|\omega|) and (\sigma=\xi^\top S\xi). Projecting the vorticity equation orthogonally to \(\xi\) gives
[
D_t\xi=(I-\xi\xi^\top)S\xi+
\frac{\nu}{|\omega|}(I-\xi\xi^\top)\Delta\omega.
]
This formula is only pointwise defined where (\omega\ne0).

The material strain equation is
[
D_tS=-S^2-W^2-\nabla^2p+\nu\Delta S.
]
Since (W\xi=0) in three dimensions when \(\xi\) is parallel to vorticity, differentiating \(\sigma=\xi^\top S\xi\) yields
[
\begin{aligned}
D_t\sigma={}&|S\xi|^2-2\sigma^2
-\xi^\top(\nabla^2p)\xi+\nu\xi^\top(\Delta S)\xi\\
&+\frac{2\nu}{|\omega|}
\big[(I-\xi\xi^\top)\Delta\omega\big]^\top S\xi.
\end{aligned}
]
This is an exact identity under the smoothness and nonzero-vorticity assumptions above. It is not a closed scalar evolution equation: the pressure Hessian and viscous terms remain.

## 2. Principal-strain eigenframe

Let \(\lambda_1\ge\lambda_2\ge\lambda_3\) be the eigenvalues of \(S\). At points where \(\lambda_1\) is simple, choose a local orthonormal eigenframe \(e_1,e_2,e_3\), and write \(\alpha_i=(\xi\cdot e_i)^2\). Then
[
\sum_i\alpha_i=1,\qquad \alpha_i\ge0,
\]
and
[
\sigma=\lambda_1\alpha_1+\lambda_2\alpha_2+\lambda_3\alpha_3.
]
Define the alignment deficit
[
d=|\omega|^2(\lambda_1-\sigma).
]
The spectral identity is
[
d=|\omega|^2\big[(\lambda_1-\lambda_2)\alpha_2+
(\lambda_1-\lambda_3)\alpha_3\big]\ge0.
]
This is algebraic and follows from the spectral decomposition. By itself it is not a new theorem.

If \(\lambda_1-\lambda_2\ge\delta>0\) on a region \(\Omega\), then
[
\int_\Omega d\,dx\ge\delta\int_\Omega
|\omega|^2(1-\alpha_1)\,dx.
]
The assumption must be stated explicitly; a uniform gap cannot be presumed in general.

For a simple eigenvalue, the eigenvector derivative has the formal expression
[
D_t e_1=\sum_{j=2}^3 e_j
\frac{e_j^\top(D_tS)e_1}{\lambda_1-\lambda_j}.
]
This requires a differentiable eigenframe and nonzero eigenvalue gaps. Small gaps can make the eigenframe coordinates ill-conditioned without proving that the fluid itself is singular. A spectral projector removes eigenvector sign ambiguity but does not remove the spectral-gap issue.

## 3. Enstrophy balance and what the deficit does not prove

Define
[
\Omega(t)=\frac12\int_{\mathbb T^3}|\omega|^2\,dx,
\qquad P(t)=\int_{\mathbb T^3}|\nabla\omega|^2\,dx.
]
For a smooth periodic solution,
[
\Omega'(t)+\nu P(t)=\int_{\mathbb T^3}\omega\cdot S\omega\,dx.
]
Since \(\omega\cdot S\omega=\sigma|\omega|^2\), the pointwise decomposition gives
[
\Omega'+\nu P=\int\lambda_1|\omega|^2\,dx-\int d\,dx,
]
provided the integrals and definitions are interpreted consistently (for example, away from degeneracies or through a carefully justified measurable spectral formulation).

This is an exact rearrangement, not a new estimate. The condition
[
\int\lambda_1|\omega|^2-\int d\le\nu P
]
implies nonincreasing enstrophy, but it merely restates the sign condition on the enstrophy derivative. It is not an independent regularity criterion.

Likewise, writing \(S=S_+-S_-\) with positive and negative spectral parts gives
[
\Omega'+\nu P=\mathcal P_+-\mathcal P_-,
\]
where \(\mathcal P_+=\int\omega\cdot S_+\omega\ge0\) and \(\mathcal P_-=\int\omega\cdot S_-\omega\ge0\). Assuming \(\mathcal P_+\le\mathcal P_-+\nu P\) proves \(\Omega'\le0\), but the assumption is the entire difficulty and is not established here.

## 4. Novelty assessment

### Established in this note
- Exact vorticity-direction and directional-stretching identities under explicit smoothness/nonzero-vorticity assumptions.
- A strain-eigenvalue decomposition of the deficit.
- A conditional lower bound under a uniform spectral gap.
- Identification of pressure-Hessian, viscous, and eigenframe-rotation terms that obstruct a simple closed estimate.

### Not established
- Any new theorem compared with the literature.
- A dynamical lower bound for the deficit.
- A closed evolution inequality for the deficit.
- Propagation of a depletion condition.
- A new conditional or unconditional regularity criterion.
- Any proof of global regularity for 3D Navier–Stokes.

### Proof obligations for a potentially new result
1. State a precise theorem with function spaces, domain, time interval, and assumptions.
2. Derive all localization/cutoff terms if using a local formulation.
3. Control pressure-Hessian and eigenframe-rotation terms without assuming the desired conclusion.
4. Test scale invariance and dimensions.
5. Search for counterexamples among exact flows and known special solutions.
6. Compare the theorem line by line with existing strain-eigenvalue and vorticity-direction criteria.
7. Do not claim novelty until the literature comparison and proof are independently checked.

## 5. Literature leads to verify

These are starting points for a proper novelty review, not evidence that the current approach is novel:

- B. Galanti, J. Gibbon, and M. Heritage, “Vorticity alignment results for the three-dimensional Euler and Navier–Stokes equations,” *Nonlinearity* 10 (1997), DOI: [10.1088/0951-7715/10/6/013](https://doi.org/10.1088/0951-7715/10/6/013).
- E. Dresselhaus and M. Tabor, “The kinematics of stretching and alignment of material elements in general flow fields,” *Journal of Fluid Mechanics* (1992), DOI: [10.1017/S0022112092001460](https://doi.org/10.1017/S0022112092001460).
- A. P. Encinas-Bartos and G. Haller, “Vorticity alignment with Lyapunov vectors and rate-of-strain eigenvectors,” *European Journal of Mechanics B/Fluids* (2024), DOI: [10.1016/j.euromechflu.2024.02.003](https://doi.org/10.1016/j.euromechflu.2024.02.003).
- Constantin–Fefferman geometric regularity criteria and subsequent vorticity-direction literature.
- Modern criteria based on strain eigenvalues, including work on the middle strain eigenvalue.

Each reference should be checked against the publisher/preprint version and read in full before the manuscript makes a novelty claim.

## 6. Immediate next experiment

Do not begin with a global regularity claim. First investigate a falsifiable, narrower question:

> Can one derive an independently useful upper bound for the positive stretching integral from a scale-localized alignment deficit plus quantities controlled by the energy class, with all pressure, cutoff, and eigenframe terms accounted for?

Start by trying to prove the estimate, then test its scaling and exact-flow cases. If it reduces to the enstrophy identity or fails for known solutions, record the failure and pivot. A failed candidate is useful research evidence; it is not a theorem.
