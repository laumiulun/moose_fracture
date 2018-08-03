# PorousFlowCapillaryPressureBW
!syntax description /UserObjects/PorousFlowCapillaryPressureBW

The Broadbridge-White capillarity relationship valid for small $K_{n}$ is \citep{broadbridge1988}
\begin{equation}
S_{\mathrm{eff}} = S_{n} + (S_{s} - S_{n}) \frac{c}{1 + L(x)} \ .
\end{equation}
where
\begin{equation}
x = (c - 1) e^{c - 1 - c P/\lambda} \ ,
\end{equation}
and $L(x)$ is the Lambert W-function that satisfies $L(z)e^{L(z)}=z$.
This is of limited use in real simulations, and is only used in the Porous
Flow module for comparison with the analytical solutions of \citet{broadbridge1988} and
\citet{warrick1990} for multi-phase infiltration and drainage problems.

!!! note
    Only effective saturation as a function of capillary pressure is available

!syntax parameters /UserObjects/PorousFlowCapillaryPressureBW

!syntax inputs /UserObjects/PorousFlowCapillaryPressureBW

!syntax children /UserObjects/PorousFlowCapillaryPressureBW

##References
\bibliographystyle{unsrt}
\bibliography{porous_flow.bib}
