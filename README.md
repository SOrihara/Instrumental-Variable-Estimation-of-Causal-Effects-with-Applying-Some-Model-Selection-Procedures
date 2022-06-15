# Instrumental-Variable-Estimation-of-Causal-Effects-with-Applying-Some-Model-Selection-Procedures
This repository stores a R code to estimate a treatment effect by 2SRI and LIML.
This repository stores a R code to estimate a treatment effect by 2SRI and LIML with model selection procedures discussed by Orihara et al. (2022). Also, a toy example is prepared in the R code.

The R code needs to prepare the following variables:
- TT: treatment variable, XX1: matrix of expectation variables for the treatment model
- YY: outcome variable,   XX2: matrix of expectation variables for the outcome model

Unfortunately, the proposed method has some limitations:
- Treatment variables are only continuous, and outcome variables are only dichotomous.
- Unmeasured covariates are only the standard binormal distribution.
- Treatment & outcome models are only linear models.
- The homogeneous effect is assumed for the outcome model.
- All variables are handled as continuous variables except for the outcome.

Also please note the following when use the R code:
- The programs has not validated correctly such as double programming (there is only self check).
- The programs are freely available; however, there are no responsibility.
