{smcl}
{* *! version 0.9, 25 April 2019}{...}
{cmd:help for rwrmed}{right:Geoffrey T. Wodtke}
{hline}

{title:Title}

{p2colset 5 18 18 2}{...}
{p2col : {cmd:rwrmed} {hline 2}}causal mediation analysis using regression-with-residuals{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 18 2}
{cmd:rwrmed} {varname}{cmd:,} avar({varname}) mvar({varname}) a({it:real}) astar({it:real}) m({it:real}) 
mreg({it:string}) [cvars({varlist}) lvars({varlist}) {opt nointer:action} 
{reps({it:integer 1000}) strata({varname}) cluster({varname}) level(cilevel) seed({it:passthru})]

{phang}{opt varname} - this specifies the outcome variable.

{phang}{opt avar(varname)} - this specifies the treatment (exposure) variable.

{phang}{opt mvar(varname)} - this specifies the mediator variable.

{phang}{opt a(real)} - this specifies the reference level of treatment.

{phang}{opt astar(real)} - this specifies the alternative level of treatment. Together, (astar - a) defines
the treatment contrast of interest.

{phang}{opt m(real)} - this specifies the level of the mediator at which the controlled direct effect 
is evaluated. If there is no treatment-mediator interaction, then the controlled direct effect
is the same at all levels of the mediator and thus an arbitary value can be chosen.

{phang}{opt mreg(string)} - this specifies the form of regression model to be estimated for the mediator variable. This 
model can be either {it:linear}, {it:logistic}, or {it:poisson}.


{title:Options}

{phang}{opt cvars(varlist)} - this option specifies the list of pre-treatment covariates to be included in the analysis. Categorical 
variables need to be coded as a series of dummy variables before being entered as covariates.

{phang}{opt lvars(varlist)} - this option specifies the list of post-treatment covariates to be included in the analysis. Categorical 
variables need to be coded as a series of dummy variables before being entered as covariates.

{phang}{opt nointer:action} - this option specifies whether a treatment-mediator interaction is not to be
included in the outcome model (the default assumes an interaction is present).

{phang}{opt reps(integer 200)} - this option specifies the number of replications for bootstrap resampling (the default is 200).

{phang}{opt strata(varname)} - this option specifies a variable that identifies resampling strata. If this option is specified, 
then bootstrap samples are taken independently within each stratum.

{phang}{opt cluster(varname)} - this option specifies a variable that identifies resampling clusters. If this option is specified,
then the sample drawn during each replication is a bootstrap sample of clusters.

{phang}{opt level(cilevel)} - this option specifies the confidence level for constructing bootstrap confidence intervals. If this 
option is omitted, then the default level of 95% is used.

{phang}{opt seed(passthru)} - this option specifies the seed for bootstrap resampling. If this option is omitted, then a random 
seed is used and the results cannot be replicated. {p_end}


{title:Description}

{pstd}{cmd:rwrmed} performs causal mediation analysis using regression-with-residuals. Two models are estimated: a 
model for the mediator conditional on treatment and the pre-treatment covariates (if specified) after centering them 
around their sample means, and a model for the outcome conditional on treatment, the mediator, the pre-treatment covariates 
(if specified) after centering them around their sample means, and the post-treatment covariates (if specified) after 
centering them around their estimated conditional means given all prior variables. It extends causal mediation 
analysis -- for example, that performed by {cmd:paramed} -- to allow for the presence of treatment-induced confounders, 
which are post-treatment covariates that confound the mediator-outcome relationship.

{pstd}{cmd:rwrmed} allows linear models for the outcome, and linear, logistic, or poisson models for the mediator. It requires the user
to specify an appropriate form for the mediator model. Linear models are used to residualize the post-treatment covariates with
respect to treatment and the pre-treatment covariates.

{pstd}{cmd:rwrmed} provides estimates of the controlled direct effect and then the randomized intervention analogues to 
the natural direct effect, the natural indirect effect, and the total effect when a set of measured post-treatment covariates
are included in {opt lvars(varlist)}. When post-treatment covariates are not included in {opt lvars(varlist)}, the command instead 
provides estimates of the conventional natural direct and indirect effects and the conventional total effect. Standard errors 
and confidence intervals are derived using the nonparametric bootstrap. See references for further details. {p_end}


{title:Assumptions}

{pstd}Let C be the measured pre-treatment covariates included in {opt cvars(varlist)}, and let L be the measured post-treatment covariates
included in {opt lvars(varlist)}. Obtaining consistent estimates of the controlled direct effect requires two main assumptions: {p_end}

{phang2}(1) There are no unmeasured treatment-outcome confounders given C {p_end}
{phang2}(2) There are no unmeasured mediator-outcome confounders given C and L {p_end}

{pstd}Obtaining consistent estimates of the randomized intervention analogues to natural direct and indirect effects requires assumptions (1) 
and (2) and then an additional assumption: {p_end}

{phang2}(3) There are no unmeasured treatment-mediator confounders given C {p_end}

{pstd}Note that assumptions (1) and (3) are satisified by random allocation of the treatment variable. See references for further details. {p_end}


{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. use rwrmed_example.dta} {p_end}

 
{pstd}Continuous outcome, continuous mediator, binary treatment coded 0 and 1, two pre-treatment covariates, two post-treatment covariates, 
no interaction between treatment and mediator, bootstrap standard errors with default settings: {p_end}
 
{phang2}{cmd:. rwrmed y_cont, avar(treat_bin) mvar(m_cont) cvars(cvar1 cvar2) lvars(lvar1 lvar2) a(0) astar(1) m(0) mreg(linear) nointer} {p_end}

 
{pstd}Continuous outcome, binary mediator, binary treatment coded 0 and 1, two pre-treatment covariates, two post-treatment covariates, 
include an interaction between treatment and mediator, bootstrap standard errors with 500 replications fixing the seed to 1234: {p_end}
 
{phang2}{cmd:. rwrmed y_cont, avar(treat_bin) mvar(m_bin) cvars(cvar1 cvar2) lvars(lvar1 lvar2) a(0) astar(1) m(0) mreg(logistic) reps(500) seed(1234)} {p_end}

 
{pstd}Continuous outcome, count mediator, continuous treatment, two pre-treatment covariates, two post-treatment covariates, include an interaction 
between treatment and mediator, bootstrap standard errors with default settings: {p_end}
 
{phang2}{cmd:. rwrmed y_cont, avar(treat_cont) mvar(m_num) cvars(cvar1 cvar2) lvars(lvar1 lvar2) a(-1) astar(1) m(0) mreg(poisson)} {p_end}


{title:Saved results}

{pstd}{cmd:rwrmed} saves the following results in {cmd:e()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}matrix containing direct, indirect and total effect estimates{p_end}
{synopt:{cmd:e(V)}}matrix containing variance of the effect estimates{p_end}


{title:Author}

{pstd}Geoffrey T. Wodtke {break}
Department of Sociology{break}
University of Toronto{p_end}

{phang}Email: geoffrey.wodtke@utoronto.ca


{title:References}

{pstd}Wodtke GT and Zhou X. (2019). Effect Decomposition in the Presence of Treatment-induced Confounding: A Regression-with-Residuals Approach. In preparation. {p_end}

{pstd}Zhou X and Wodtke GT. (2018). A Regression-with-residuals Method for Estimating Controlled Direct Effects. Political Analysis. Advance online publication. https://doi.org/10.1017/pan.2018.53 {p_end}


{title:Acknowledgments}

{pstd}This command is based on the {cmd:paramed} command by Hanhua Liu, Richard Emsley, Graham Dunn, Tyler VanderWeele, and Linda Valeri. {p_end}


{title:Also see}

{psee}
Help: {manhelp regress R}, {manhelp logit R}, {manhelp poisson R}, {manhelp bootstrap R}, {manhelp paramed R}
{p_end}
