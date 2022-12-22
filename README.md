# Angular-Correlation

Overview of Angular Correlations programs used for the analysis of ILL Dataset.

Note many of these codes are set to angles respective to the detector configuration of 155mm of the IFIN clover to Target distance and 134mm of the FIPPS clover to Target distance.

## FitCode (Attenuation Factor Extraction)

This Program is utilized to Evaluate the Attenuation factors given known/pure Transitions.

Two Methods are utilized in the program. The first labelled as *Sambu* utilizes the truncated Angular Correlation equation: 
```math
W(\theta) = 1 + a_2P_2(cos(\theta)) + a_4P_4(cos(\theta))
```

```math
a_2 = Q_2B_2(J_2)R_2(J_2J_4)
```
```math
a_4 = Q_4B_4(J_2)R_4(J_2J_3) 
```

Where Q is the solid-angle correction factor related to the finite size of detector size and the Coefficient of B and R factors are adopted from [H. J. Rose and D. M. Brink. Angular distribution of γ rays in terms of phase-defined matrix
elements. Rev. Mod. Phys., 39, (1967).] The function sets free the Q values and performs a fit on the experimental angular correlations and yields the fitted Q factors.

The Second Method labelled as *Norm C* imploys the same Angular equation as above, however, the a2 and a4 factors are evaluated using the online GRIFFIN TOOL Angular Correlation calculator (griffincollaboration.github.io/AngularCorrelationUtility/)
however, the a2 and a4 were multiplyed by their Q factors respectively, which once again are set free in the program and are fitted for.


## 
