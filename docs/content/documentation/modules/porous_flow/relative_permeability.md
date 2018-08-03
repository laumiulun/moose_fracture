# Relative permeability

The relative permeability of a phase is a function of its effective
saturation:
\begin{equation}
S_{\mathrm{eff}}(S) = \frac{S - S_{\mathrm{res}}^{\beta}}{1 -
  \sum_{\beta'}S_{\mathrm{res}}^{\beta'}}
\end{equation}
In this equation $S_{\mathrm{res}}^{\beta}$ is the residual
saturation for phase $\beta$.  If $S_{\mathrm{eff}} < 0$ then the
relative permeability is zero, while if $S_{\mathrm{eff}}>1$ then the
relative permeability is unity.  Otherwise, the relative permeability
is given by the expressions below.

## Constant
[`PorousFlowRelativePermeabilityConst`](/porous_flow/PorousFlowRelativePermeabilityConst.md)

The relative permeability of the phase is constant
\begin{equation}
k_{\mathrm{r}} = C.
\end{equation}
This is not recommended because there is nothing to discourage phase
disappearance, which manifests itself in poor convergence.  Usually
$k_{\mathrm{r}}(S) \rightarrow 0$ as $S\rightarrow 0$ is a much better
choice. However, it is included as it is useful for testing purposes.

## Corey
[`PorousFlowRelativePermeabilityCorey`](/porous_flow/PorousFlowRelativePermeabilityCorey.md)

The relative permeability of the phase is given by \citet{corey1954}
\begin{equation}
k_{\mathrm{r}} = S_{\mathrm{eff}}^{n},
\end{equation}
where $n$ is a user-defined quantity. Originally, \citet{corey1954} used $n = 4$ for the wetting phase, but the PorousFlow module allows an arbitrary exponent to be used.

!media media/porous_flow/relperm_corey.png width=80% margin-left=10px caption=Corey relative permeability id=relperm_corey


## van Genuchten
[`PorousFlowRelativePermeabilityVG`](/porous_flow/PorousFlowRelativePermeabilityVG.md)

The relative permeability of the phase is given by \citet{vangenuchten1980}
\begin{equation}
k_{\mathrm{r}} = \sqrt{S_{\mathrm{eff}}} \left(1 - (1 -
S_{\mathrm{eff}}^{1/m})^{m} \right)^{2}.
\end{equation}
This has the numerical disadvantage that its derivative as
$S_{\mathrm{eff}}\rightarrow 1$ tends to infinity.  This means that
simulations where the saturation oscillates around
$S_{\mathrm{eff}}=1$ do not converge well.

Therefore, a *cut* version of the van Genuchten expression is also offered, which is
almost definitely indistinguishable experimentally from the original expression:
\begin{equation}
k_{\mathrm{r}} =
\begin{cases}
\textrm{van Genuchten} & \textrm{for } S_{\mathrm{eff}} < S_{c} \\
\textrm{cubic} & \textrm{for } S_{\mathrm{eff}} \geq S_{c}.
\end{cases}
\end{equation}
Here the cubic is chosen so that its value and derivative match the
van Genuchten expression at $S=S_{c}$, and so that it is unity at
$S_{\mathrm{eff}}=1$.

!media media/porous_flow/relperm_vg.png width=80% margin-left=10px caption=van Genuchten relative permeability id=relperm_vg

## Brooks-Corey
[`PorousFlowRelativePermeabilityBC`](/porous_flow/PorousFlowRelativePermeabilityBC.md)

The \citet{brookscorey1966} relative permeability model is an extension of the previous  \citet{corey1954} formulation where the relative permeability of the wetting phase is given by
\begin{equation}
k_{\mathrm{r, w}} = \left(S_{\mathrm{eff}}\right)^{(2 + 3 \lambda)/\lambda},
\end{equation}
and the relative permeability of the non-wetting phase is
\begin{equation}
k_{\mathrm{r, nw}} = (1 - S_{\mathrm{eff}})^2 \left[1 - \left(S_{\mathrm{eff}}\right)^{(2 + \lambda)/\lambda}\right],
\end{equation}
where $\lambda$ is a user-defined exponent. When $\lambda = 2$, this formulation reduces
to the original \citet{corey1954} form.

!media media/porous_flow/relperm_bc.png width=80% margin-left=10px caption=Brooks-Corey relative permeability id=relperm_bc

## Broadbridge-White
[`PorousFlowRelativePermeabilityBW`](/porous_flow/PorousFlowRelativePermeabilityBW.md)

The relative permeability of a phase given by \citet{broadbridge1988} is
\begin{equation}
k_{\mathrm{r}} = K_{n} + \frac{K_{s} - K_{n}}{(c - 1)(c -
  S_{\mathrm{eff}})}S_{\mathrm{eff}}^{2}.
\end{equation}

## FLAC
[`PorousFlowRelativePermeabilityFLAC`](/porous_flow/PorousFlowRelativePermeabilityFLAC.md)

A form of relative permeability used in FLAC, where
\begin{equation}
k_{\mathrm{r}} = (m + 1)S_{\mathrm{eff}}^{m} - m S_{\mathrm{eff}}^{m + 1}.
\end{equation}
This has the distinct advantage over the Corey formulation that its
derivative is continuous at $S_{\mathrm{eff}}=1$.

##References
\bibliographystyle{unsrt}
\bibliography{porous_flow.bib}
