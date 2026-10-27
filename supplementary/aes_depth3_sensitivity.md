# S2. Parameter sensitivity for AES MixColumns at depth three

For AES MixColumns at depth three, each pair (R, W) is evaluated with L = 3000 and a nominal 1,800-second time limit per MILP solve. The table reports the best verified feasible circuit returned for each configuration.

- R: number of randomized seed traces.
- W: maximum number of expansion decompositions retained per target. Existing seed decompositions are not subject to this limit.
- OPT: the MILP proved optimality over the generated candidate space. This does not establish global circuit optimality.
- TL: a verified incumbent at the time limit.
- n2 and n3: numbers of XOR2 and XOR3 gates.
- Time (s) reproduces the reported elapsed time, including entries slightly above the nominal time limit.

| R | W | Status | GE | (n2, n3) | Time (s) |
| --- | --- | --- | --- | --- | --- |
| 1 | 64 | OPT | 258.46 | (27, 38) | 4.6 |
| 1 | 128 | OPT | 254.13 | (24, 39) | 23.4 |
| 1 | 256 | OPT | 252.46 | (25, 38) | 346.4 |
| 1 | 512 | OPT | 250.45 | (29, 35) | 757.1 |
| 1 | 1024 | TL | 251.44 | (34, 32) | 1800.0 |
| 3 | 64 | OPT | 249.12 | (27, 36) | 169.5 |
| 3 | 128 | TL | 249.12 | (27, 36) | 1800.0 |
| 3 | 256 | TL | 248.13 | (22, 39) | 1802.1 |
| 3 | 512 | TL | 248.13 | (22, 39) | 1800.0 |
| 3 | 1024 | TL | 255.12 | (29, 36) | 1800.0 |
| 5 | 64 | OPT | 249.12 | (27, 36) | 1607.9 |
| 5 | 128 | TL | 247.45 | (28, 35) | 1800.0 |
| 5 | 256 | TL | 247.45 | (28, 35) | 1800.0 |
| 5 | 512 | TL | 247.45 | (28, 35) | 1800.4 |
| 5 | 1024 | TL | 255.12 | (29, 36) | 1800.0 |
| 10 | 64 | TL | 247.11 | (31, 33) | 1800.0 |
| 10 | 128 | TL | 247.45 | (28, 35) | 1800.0 |
| 10 | 256 | TL | 247.45 | (28, 35) | 1800.0 |
| 10 | 512 | TL | 244.45 | (27, 35) | 1800.7 |
| 10 | 1024 | TL | 253.45 | (30, 35) | 1800.0 |
| 20 | 64 | TL | 247.79 | (25, 37) | 1804.6 |
| 20 | 128 | TL | 245.78 | (29, 34) | 1800.0 |
| 20 | 256 | TL | 247.45 | (28, 35) | 1800.0 |
| 20 | 512 | TL | 251.78 | (31, 34) | 1800.0 |
| 20 | 1024 | TL | 256.45 | (31, 35) | 1800.0 |

[Repository overview](../README.md)
