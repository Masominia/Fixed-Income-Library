# Methods

Times are in years, zero rates are continuously compounded, P(t) is the discount factor.

## Curve

The curve stores zero rates at pillar times and interpolates y(t) = -ln P(t) with y(0) = 0.

- **Linear zero rates.** z(t) is linear between pillars and constant outside.
- **Log-linear.** y(t) is linear between pillars, including the segment from the origin: forward rates are constant on each segment. Beyond the last pillar the last forward continues.
- **Monotone cubic.** y(t) is a cubic Hermite spline with Fritsch-Carlson slopes (weighted harmonic mean of neighbouring secants, zero at a local extremum, limited one-sided slopes at the ends). If y increases at the pillars it increases everywhere. Beyond the last pillar the last slope continues.

The instantaneous forward rate is the analytic derivative of the interpolant, not a finite difference.

## Bootstrap

Each quote gives one pillar at its maturity and one equation: the rate implied by the curve equals the quoted rate. Implied rates are `(1/P(T) - 1)/T` for a deposit, `(P(s)/P(e) - 1)/(e - s)` for a forward rate agreement, and `(1 - P(T)) / annuity` for a swap. The pillars are solved one at a time by Brent's method with the others held fixed, and the sweep over all pillars is repeated until no zero rate moves by more than 1e-14.

## Bonds

With yield y compounded f times a year, `P = sum c_i (1 + y/f)^(-f t_i)`. Macaulay duration is the present-value-weighted mean time, modified duration divides it by `1 + y/f`, convexity is `sum t_i (t_i + 1/f) pv_i / (P (1 + y/f)^2)`. Effective and key-rate durations move zero rates by one basis point up and down, all together or one pillar at a time; with linear zero rates the one-pillar moves add up to the parallel move exactly.

## Nelson-Siegel

`z(t) = level + slope (1 - e^-x)/x + curvature ((1 - e^-x)/x - e^-x)`, x = t / tau. For fixed tau the coefficients solve a 3 by 3 linear system; tau is scanned on a logarithmic grid and refined by golden section.

## Black and Bachelier

Standard formulas on the forward. A caplet over [s, e] is `P(e) (e - s) Black(F, K, vol, s)`. A swaption is `annuity x Black(par rate, K, vol, expiry)`. Implied volatility is found by Brent's method after checking the no-arbitrage bounds.

## Short-rate models

Vasicek: `P(T) = exp(A - B r)`, `B = (1 - e^-aT)/a`, `A = (B - T)(a^2 b - sigma^2/2)/a^2 - sigma^2 B^2/(4a)`. Cox-Ingersoll-Ross: the usual formula with `g = sqrt(a^2 + 2 sigma^2)`. Simulation uses antithetic pairs and the trapezoid rule for the integral of the rate; the Vasicek transition is exact, the Cox-Ingersoll-Ross one is Euler with the rate floored at zero inside drift and diffusion.

## Hull-White

`dr = (theta(t) - a r) dt + sigma dW`, with theta implied by the initial curve.

- Bond at time t given r: `P(t,T) = P(0,T)/P(0,t) exp(B f(0,t) - sigma^2 (1 - e^-2at) B^2/(4a) - B r)`, B = B(T - t).
- Option on a zero-coupon bond: Black-type formula with `sigma_p = sigma sqrt((1 - e^-2aT)/(2a)) B(S - T)`.
- Caplet: `(1 + K tau)` puts on the bond maturing at the payment time, struck at `1/(1 + K tau)`.
- Swaption: all bond prices fall as r rises, so there is one rate r* at which the fixed leg is worth par. The swaption is a portfolio of bond options, one per payment, struck at the bond prices at r*.

## Tree

Equal time steps, spacing `sigma sqrt(3 dt)`, standard three-way branching, with the two alternative branchings at the edges once the tree reaches `ceil(0.184/(a dt))` nodes from the centre. The shift of each time slice is chosen by forward induction on Arrow-Debreu prices so that the next discount factor of the curve is matched.

Swaptions roll two value arrays back together: the remaining fixed leg including the notional, and the option. At an exercise time the option becomes the larger of its continuation value and par minus the fixed leg (payer) or the reverse (receiver); the coupon paid at that time is added to the fixed leg after the comparison, since it belongs to the period before. Payment and exercise times must fall on the time grid.

## Calibration

The quotes are Black volatilities of at-the-money swaptions, converted to prices. The objective is the sum of squared relative price errors. A coarse scan over the logarithms of a and sigma is followed by a pattern search that halves its step when no move improves.
