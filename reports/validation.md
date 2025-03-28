# Validation checks

| area | check | status | detail |
|---|---|---|---|
| curve | every quote repriced under every interpolation | PASS | largest error 3.1e-16 |
| curve | the cubic removes the forward jumps | PASS | 8.0 bp between points 0.05 years apart, against 32.9 (log-linear) and 45.5 (linear zero) |
| curve | forward rates stay positive | PASS | lowest 2.95% |
| bond | key-rate durations add up to the effective duration | PASS | sum 8.313628 against 8.313629 |
| bond | convexity removes most of the error at 100 bp | PASS | 0.403 down to 0.0139 per 100 of face |
| bond | duration alone always understates the price | PASS | negative error for every shift |
| curve | Nelson-Siegel fits the zero rates within 10 bp | PASS | largest error 6.3 bp, root mean square 2.7 bp |
| Hull-White | tree reproduces the curve | PASS | largest relative error of a discount factor 3.1e-15 |
| Hull-White | closed-form payer and receiver satisfy parity | PASS | gap -1.1e-16 |
| Hull-White | at-the-money parity holds on the tree | PASS | largest gap 0.0e+00 |
| Hull-White | finest tree within 0.1% of the closed forms | PASS | largest relative error 1.3e-04 at 512 steps per year |
| Hull-White | Bermudan worth more than European | PASS | ratio 1.360 |
| calibration | parameters recovered from their own volatilities | PASS | a 0.050000, sigma 0.0100000 |
| calibration | stable under perturbed quotes | PASS | a 0.0523, sigma 0.01009 |
| calibration | volatility falls with expiry and tenor | PASS | 26.1% to 15.4% |
| short rate | Vasicek formula within Monte Carlo error | PASS | largest |z| 1.19 |
| short rate | Cox-Ingersoll-Ross formula within 2 bp of an Euler simulation | PASS | largest gap 1.32 bp of price |
