# Interventions_SES_PROMs
Code to perform the analyses in Lindmark, A. & Darehed, D. (2024). "Investigating hypothetical interventions to mitigate socioeconomic differences in Patient Reported Outcome Measures 3 months post-stroke among working-age patients: A nationwide register-based study"

Contains the following code files:
- "estimation_func.R" with the code for the function to estimate the interventional direct and indirect effects as described in Section 1.2 in the Supplementary Materials.
- "call.R" the function call used to estimate the interventional direct and indirect effects and bootstrap standard errors.

Note that the code is written for the specific example in the article. It can be adapted to other examples, but if you only have single mediators (i.e. no mediator groups that you wish to intervene on) a more general implementation can be found at https://github.com/moreno-betancur/medRCT.
