# CurveKit report

Times in years. Quotes in `configs/quotes.txt`: two deposits, one forward rate agreement, nine swaps out to 30 years.

## 1. Bootstrap

|  | sweeps | largest repricing error | largest forward jump, bp | lowest forward, % | highest forward, % |
|---|---|---|---|---|---|
| linear zero rates | 2 | 3.05e-16 | 45.48 | 3.09 | 4.51 |
| log-linear discount factors | 2 | 3.05e-16 | 32.88 | 3.09 | 4.39 |
| monotone cubic | 6 | 3.05e-16 | 8.05 | 2.95 | 4.43 |

Forward jump: largest change of the instantaneous forward rate between two points 0.05 years apart.

Pillars (monotone cubic):

| maturity, years | zero rate | discount factor |
|---|---|---|
| 0.250000 | 0.030880 | 0.992310 |
| 0.500000 | 0.032239 | 0.984010 |
| 1.000000 | 0.033566 | 0.966991 |
| 2.000000 | 0.035226 | 0.931972 |
| 3.000000 | 0.036541 | 0.896171 |
| 4.000000 | 0.037669 | 0.860128 |
| 5.000000 | 0.038602 | 0.824472 |
| 7.000000 | 0.039972 | 0.755933 |
| 10.000000 | 0.041151 | 0.662649 |
| 15.000000 | 0.042005 | 0.532553 |
| 20.000000 | 0.041604 | 0.435143 |
| 30.000000 | 0.040083 | 0.300442 |

## 2. Bond: 10 years, 4.0% coupon

Price 99.0172, yield 4.1209%. Macaulay duration 8.3293, modified 8.1611, convexity 78.682, value of a basis point 0.08081. Effective duration on the curve 8.3136.

| pillar, years | key-rate duration |
|---|---|
| 0.5000 | 9.94e-03 |
| 1.0000 | 0.0339 |
| 2.0000 | 0.0751 |
| 3.0000 | 0.1084 |
| 4.0000 | 0.1388 |
| 5.0000 | 0.2620 |
| 7.0000 | 0.5491 |
| 10.0000 | 7.1363 |

Repricing after a yield shift:

| yield shift, bp | exact price | duration only | duration and convexity | error, duration only | error, with convexity |
|---|---|---|---|---|---|
| -200.0000 | 116.8519 | 115.1791 | 116.7373 | -1.6728 | -0.1146 |
| -100.0000 | 107.5016 | 107.0982 | 107.4877 | -0.4035 | -0.0139 |
| -50.0000 | 103.1568 | 103.0577 | 103.1551 | -0.0991 | -1.72e-03 |
| -10.0000 | 99.8292 | 99.8253 | 99.8292 | -3.91e-03 | -1.36e-05 |
| 10.0000 | 98.2130 | 98.2091 | 98.2130 | -3.88e-03 | 1.35e-05 |
| 50.0000 | 95.0725 | 94.9768 | 95.0742 | -0.0957 | 1.67e-03 |
| 100.0000 | 91.3127 | 90.9363 | 91.3259 | -0.3764 | 0.0132 |
| 200.0000 | 84.3109 | 82.8554 | 84.4136 | -1.4556 | 0.1026 |

## 3. Nelson-Siegel fit

Level 3.4928%, slope -0.4003%, curvature 2.9002%, decay time 6.83 years. Root mean square error 2.66 bp, largest 6.35 bp.

## 4. Hull-White, a = 0.05, sigma = 0.01

Option expiring in 2 years on a 5-year swap at the forward rate 4.2214%. Closed forms: payer swaption 0.020274, zero-coupon bond option 0.017960, caplet 0.001813 (Black volatility 23.74%).

Relative error of the tree against the closed forms:

| steps per year | nodes across | bond option | caplet | swaption | tree payer - receiver | Bermudan / European |
|---|---|---|---|---|---|---|
| 4 | 31 | -0.01129 | 0.03110 | -0.01156 | 0.00000 | 1.38487 |
| 8 | 61 | -2.90e-03 | 0.01954 | -3.08e-03 | 0.00000 | 1.37321 |
| 16 | 119 | 1.92e-04 | 7.69e-03 | 7.67e-05 | 0.00000 | 1.36373 |
| 32 | 237 | 1.09e-03 | -2.40e-03 | 1.02e-03 | 0.00000 | 1.36079 |
| 64 | 473 | 1.10e-03 | 2.08e-03 | 1.06e-03 | 0.00000 | 1.35970 |
| 128 | 945 | 8.03e-04 | -2.96e-04 | 7.90e-04 | 0.00000 | 1.35951 |
| 256 | 1,887 | 4.44e-04 | 4.11e-04 | 4.47e-04 | 0.00000 | 1.35985 |
| 512 | 3,771 | 1.15e-04 | 5.58e-05 | 1.29e-04 | 0.00000 | 1.36016 |

## 5. Swaption volatilities implied by the model, and calibration back

| expiry \ tenor | 1 | 2 | 5 | 10 |
|---|---|---|---|---|
| 1.0000 | 0.2611 | 0.2472 | 0.2164 | 0.1872 |
| 2.0000 | 0.2406 | 0.2293 | 0.2048 | 0.1799 |
| 5.0000 | 0.2043 | 0.1984 | 0.1836 | 0.1652 |
| 10.0000 | 0.1804 | 0.1756 | 0.1652 | 0.1542 |

|  | exact quotes | quotes perturbed by 2% |
|---|---|---|
| mean reversion | 0.050000 | 0.052276 |
| sigma | 1.00e-02 | 0.010094 |
| volatility error, bp | 2.08e-04 | 25.047800 |

## 6. Short-rate models against Monte Carlo

| model | maturity | closed form | Monte Carlo | standard error | difference, bp of price | z |
|---|---|---|---|---|---|---|
| Vasicek | 1.000000 | 0.977554 | 0.977554 | 5.92e-08 | 2.95e-04 | 0.498010 |
| Vasicek | 2.000000 | 0.951424 | 0.951424 | 3.76e-07 | -1.24e-03 | -0.329566 |
| Vasicek | 5.000000 | 0.863221 | 0.863217 | 3.06e-06 | -0.035668 | -1.166453 |
| Vasicek | 10.000000 | 0.717207 | 0.717196 | 9.64e-06 | -0.114599 | -1.188761 |
| Cox-Ingersoll-Ross | 1.000000 | 0.977553 | 0.977545 | 1.64e-06 | -0.076810 | -4.690053 |
| Cox-Ingersoll-Ross | 2.000000 | 0.951426 | 0.951398 | 5.31e-06 | -0.275139 | -5.184534 |
| Cox-Ingersoll-Ross | 5.000000 | 0.863384 | 0.863290 | 1.80e-05 | -0.940406 | -5.235179 |
| Cox-Ingersoll-Ross | 10.000000 | 0.718265 | 0.718133 | 3.04e-05 | -1.318259 | -4.337090 |

## Validation

17 PASS, 0 WARN, 0 FAIL. See `validation.md`.
