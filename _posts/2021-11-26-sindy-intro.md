---
layout: post
title:  "Intro to SINDy"
date:   2021-11-26 17:21:40 +0200
is_series: true
series_title: "SIR-SINDy"
---


<em>Sparse Identification of Nonlinear Dynamics</em> (SINDy) is a systems identification method introduced in [Brunton et al., 2016].
Here the authors postulate that most dynamical systems given by a (autonomous) ordinary differential equation (ODE)

$$ {d\over dt}\,s(t) = \dot s(t)= f(s(t)) \qquad (1) $$

(where \\( s(t) \in \mathbb R \\) is assumed to be at least once differentiable with respect to time \\( t \\))
typically have relatively few terms on the right-hand side, given by the function \\( f \\). And thus, given time-series observations of \\( s(t) \\), written here as \\( s_{t_i},\; i \in \{1,2,\ldots,n-1,n\} \\),
one might be able to use sparse linear regression techniques to obtain a pointwise approximation of the derivative \\(\dot s(t) \\) in (1).<br>

This results in a sparse linear regression estimator \\(\xi \\) (vector of coefficients with few non-zero entries),
by which one can approximate the function \\( f \\).<br>
Finding a sparse estimator \\(\xi \\) corresponds to finding a parsimonious model for the process \\( s(t) \\).

If the derivatives \\( \dot s_{t_i},\; i \in \{1,2,\ldots,n-1,n\} \\) are not part of the observations (which they seldom are), then these need to be approximated using, for example, finite difference, or the derivatives of a fitted smoothing spline, because these make out the target in the sparse linear regression.

Let's elaborate a little on this:<br>
Say you have \\(n\\) observations \\( s_{t_i},\; i \in \{1,2,\ldots,n-1,n\} \\) of some dynamical process \\(s(t)\\)
that is assumed to be differentiable with respect to time \\(t\\),
and that you want to approximate the governing equation of \\(s(t)\\) in the form of (1).

When using SINDy we assume that \\( f(s_{t_i}) \\) can be approximated as

$$ f(s_{t_i}) \approx [\Theta\,\xi]_i \qquad (i'\text{th entry in vector})$$

where the \\( \Theta \in \mathbb R^{n \times p} \\) matrix has nonlinear functions of the observations as columns,
and \\( \xi \in \mathbb R^p \\) is a column vector of coefficients to be determined.<br>
For example using polynomial basis function in \\( \Theta \\) up to 4th order (no intercept)

$$ 
\Theta = \begin{bmatrix} s_{t_1} & s_{t_1}^2 & s_{t_1}^3 & s_{t_1}^4 \\
                            s_{t_2} & s_{t_2}^2 & s_{t_2}^3 & s_{t_2}^4\\
                            \vdots & \vdots & \vdots \\
                            s_{t_n} & s_{t_n}^2 & s_{t_n}^3 & s_{t_n}^4 \\
\end{bmatrix} 
$$

where we have \\( p = 4 \\).

Let \\( \dot S \in \mathbb R^n \\) be the column vector of derivatives, i.e. \\( \dot S = (\dot s_{t_1}, \dot s_{t_2},\ldots, \dot s_{t_n})^\top \\),
then these can be written down (according to our assumptions) as

$$ \dot S = \Theta\,\xi + \epsilon$$

where \\( \epsilon \\) is the residual vector.

If one estimates \\( \xi \\) (minimizes the residual vector \\( \epsilon \\) ) using ordinary least squares (OLS), i.e.

$$ 
\xi^* = \text{arg}\min\limits_{\xi}\, ||\dot S - \Theta\,\xi||_2^2
$$

then the solution 
$$ 
\xi^*=(\Theta^\top \Theta)^{-1}\Theta^\top \dot S
$$ 
will have \\( p \\) non-zero entries. Nothing sparse about that...<br>
Using SINDy, the point is rather to use sparse regression, that results in \\( \xi^* \\) having  \\( q,\; 0 < q < p \\) non-zero entries,
corresponding to \\( q \\) terms in the governing equation. <br>
If \\( q \\) is small, and the fit is acceptable, then we have succeeded in finding
a parsimonious approximation of the governing equations of \\( s(t) \\), i.e.

$$ \dot s(t) \approx \theta(s(t))^\top\,\xi^* $$

where \\( \theta(s(t))^\top \\) is a row vector of nonlinear functions of \\( s(t) \\).<br>
I.e. using up to 4th order polynomial basis functions as above we have

$$ \dot s(t) \approx ( s(t), s(t)^2, s(t)^3, s(t)^4 )\,\xi^*  $$

which are the same functions as in the columns of \\( \Theta \\).

The <em>Least Absolute Shrinkage and Selection Operator</em> (LASSO) [Tibshirani, 1996] is an ideal candidate method to use in the sparse regression step,
which also is noted by the authors of [Brunton et al. 2016].<br>
The LASSO finds the estimator by minimizing

$$ \xi^* = \text{arg}\min\limits_{\xi}\, ||\dot S - \Theta\,\xi||_2^2 + \lambda\,||\xi||_1 \qquad (2) $$

where \\( \lambda > 0 \\) imposes a penalty on the \\( \ell_1 \\) norm of the estimator. Using the \\( \ell_1 \\) norm results
in (depending on the magnitude of \\( \lambda \\)) some entries of \\( \xi^* \\) being effectively zero.
Thereby, the LASSO acts as a model selection method parameterized by \\( \lambda \\).
Note that the LASSO estimator converges to the OLS estimator as \\( \lambda \rightarrow 0 \\).

The solution to (2) has to be found numerically, and one can for example use the <samp>lars</samp> or <samp>glmnet</samp> packages in <samp>R</samp> to achieve this.<br>
For estimating derivatives one can for example use <samp>smooth.spline</samp> from the <samp>stats</samp> package.

## Conclusion
SINDy is a powerful, yet simple, method to approximate the governing equations of systems, for which there is some data available.<br>
In addition, is is quite modular, in the sense that there are many options available for estimating derivatives, constructing the \\( \Theta \\) matrix,
and obtaining the linear estimator \\( \xi \\).<br>
A crucial question is which basis functions to use in \\( \Theta \\). Here physical knowledge of the system under investigation can provide some clues.<br>
Obtaining a system of ordinary differential equations (ODEs) can be achieved (in a sequential manner) by a straight-forward generalization of the above &mdash;
just include all predictors (and functions/interactions thereof) in the \\( \Theta \\) matrix, and fit the derivatives using for example LASSO.

In this series of blog articles, we will be using SINDy to obtain a parsimonious model of the <a href="{% post_url 2021-12-14-networked-sir %}">Networked SIR Model</a>, i.e. we'll take data from simulations of a stochastic model, and obtain a continuous and deterministic model by way of SINDy. In doing so, we hope via dynamical system techniques, to obtain some insights regarding the effects of network topology on the dynamics of the simulation.

### References:
[Brunton et al. 2016] Brunton, S. L., Proctor, J. L., Kutz, J. N. & Bialek, W.
Discovering governing equations from data by sparse identification of nonlinear dynamical systems. Proc. Natl. Acad. Sci. U. S. A. 113, 3932–3937 (2016)

[Tibshirani, 1996] Tibshirani, R. Regression Shrinkage and Selection Via the Lasso. J. R. Stat. Soc. Ser. B 58, 267–288 (1996).
