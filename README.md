The _Estimation and Projection Package_ (EPP) is a simple epidemiological model for inferring HIV incidence trends from population HIV surveillance data. The EPP model has been developed for over fifteen years by a large number of collaborators to include a complex representation of demography, HIV natural history, and intricacies of HIV treatment scale-up required to model the HIV epidemic with the level of realism required for national planning purposes. The model has principally been developed for fitting to survey data about HIV prevalence in the general population and among pregnant women attending ANC, though in principle could be calibrated to any variety of data sources with suitable extensions to appropriately model the data generating process, such as HIV deaths, HIV diagnoses, AIDS cases, or recent infections.

The objective of `simpleepp` is to create a stylized version of the EPP model that strips away unnecessary complexities in order to focus on the underlying inference problem. `simpleepp` will provide a conceptual introduction to EPP for new users, and be a useful tool for stylized exploration of how alternative modelling choices may affect epidemic inference, such as:

* Alternative specifications and functional forms for HIV transmission or incidence rate.
* Various configurations of potentially available data sources for inference.
* Incorporation of hierarchical models, spatial autocorrelation, and covariates predictors for the HIV transmission rate.
* Influence of considering uncertainty about various model processes, such as HIV natural history or effects of ART on mortality.
* Mis-specification of the data model, for example unaccounted for biases in particular data sources.


## Current functionality

The current implementation represents a susceptible population and progression through four stages of infection. Effects of antiretroviral treatment on survival and transmission are not yet implemented (ART compartments are present but ART initiation is not yet implemented). The 

The model is implemented as discrete difference approximation to the differential equations, with timestep `dt=0.1`. The model is implemented in Stan and returns a matrix with three columns consisting of the transmission rate, incidence rate (force of infection), and prevalence at each time step. Presently, the example script `simpleepp.R` illustrates maximum-likelihood optimization of parameter estimates. Note that optimization will tend to oversmooth penalized semiparametric components.

The following variants of the EPP model are implemented:

* EPP classic
* r-logistic 
* random walk (r-flex) with first order or second order penalty
* r-spline with first order or second order penalty

In the implementation of the r-spline model, the spline function is specified for $\log\kappa(t)$, rather than directly specified for $\kappa(t)$ as in the implemented EPP model developed by Hogan and colleagues. This has implications for the second-order smoothness because the constraint $\kappa(t) > 0$ forces a bend in the spline to avoid negative values when specifying the spline directly for $\kappa(t)$. This constraint is not in place when modelling $\log\kappa(t)$ allowing somewhat stronger penalization.

## To do:

The following features are still to be implemented:

* Inference in Stan (for optimization and MCMC inference).
* ART initiation and effects on transmission.
* Semiparametric specification of the force of infection ($\lambda(t)$).
* The following models: r-spline with equilibrium prior, r-trend model, logistic/random-walk, logistic-AR1.
* Documentation and text about all of it.
