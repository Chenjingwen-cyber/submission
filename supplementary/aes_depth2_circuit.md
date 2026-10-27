# S3. Depth-two implementation of AES MixColumns

The following table gives the complete computation sequence reported for the depth-two AES MixColumns circuit. It uses 54 XOR3 gates, with a reported h-XOR cost of 252.18 GE under the UMC 180 nm cost model.

Inputs are x0–x31, intermediate signals are t0–t21, and outputs are y0–y31. The symbol ⊕ denotes XOR. Input signals have depth zero; the reported gate depths are shown in the final column. Bit-indexing conventions follow the accompanying manuscript.

| No. | Operation | Depth |
| --- | --- | --- |
| 0 | t0 = x0 ⊕ x7 ⊕ x9 | 1 |
| 1 | t1 = x2 ⊕ x7 ⊕ x11 | 1 |
| 2 | t2 = x8 ⊕ x15 ⊕ x17 | 1 |
| 3 | t3 = x10 ⊕ x15 ⊕ x19 | 1 |
| 4 | t4 = x3 ⊕ x7 ⊕ x20 | 1 |
| 5 | t5 = x4 ⊕ x13 ⊕ x21 | 1 |
| 6 | t6 = x5 ⊕ x14 ⊕ x22 | 1 |
| 7 | t7 = x4 ⊕ x19 ⊕ x23 | 1 |
| 8 | t8 = x15 ⊕ x16 ⊕ x24 | 1 |
| 9 | t9 = x2 ⊕ x10 ⊕ x25 | 1 |
| 10 | t10 = x16 ⊕ x23 ⊕ x25 | 1 |
| 11 | t11 = x9 ⊕ x18 ⊕ x26 | 1 |
| 12 | t12 = x18 ⊕ x23 ⊕ x27 | 1 |
| 13 | t13 = x11 ⊕ x15 ⊕ x28 | 1 |
| 14 | t14 = x5 ⊕ x20 ⊕ x29 | 1 |
| 15 | t15 = x7 ⊕ x15 ⊕ x30 | 1 |
| 16 | t16 = x6 ⊕ x21 ⊕ x30 | 1 |
| 17 | t17 = x0 ⊕ x8 ⊕ x31 | 1 |
| 18 | t18 = x14 ⊕ x23 ⊕ x31 | 1 |
| 19 | t19 = x1 ⊕ x24 ⊕ x31 | 1 |
| 20 | t20 = x3 ⊕ x26 ⊕ x31 | 1 |
| 21 | t21 = x12 ⊕ x27 ⊕ x31 | 1 |
| 22 | y0 = x7 ⊕ x8 ⊕ t8 | 2 |
| 23 | y8 = x0 ⊕ x23 ⊕ t8 | 2 |
| 24 | y1 = t0 ⊕ t2 ⊕ x25 | 2 |
| 25 | y26 = x1 ⊕ x18 ⊕ t9 | 2 |
| 26 | y9 = x1 ⊕ t2 ⊕ t10 | 2 |
| 27 | y2 = x1 ⊕ x10 ⊕ t11 | 2 |
| 28 | y10 = x2 ⊕ x17 ⊕ t11 | 2 |
| 29 | y18 = x17 ⊕ t9 ⊕ x26 | 2 |
| 30 | y3 = t1 ⊕ t3 ⊕ x27 | 2 |
| 31 | y11 = x3 ⊕ t3 ⊕ t12 | 2 |
| 32 | y4 = x12 ⊕ t4 ⊕ t13 | 2 |
| 33 | y29 = x5 ⊕ t5 ⊕ x28 | 2 |
| 34 | y12 = x20 ⊕ t7 ⊕ t13 | 2 |
| 35 | y5 = x12 ⊕ t5 ⊕ x29 | 2 |
| 36 | y13 = x12 ⊕ x21 ⊕ t14 | 2 |
| 37 | y30 = x6 ⊕ t6 ⊕ x29 | 2 |
| 38 | y21 = x13 ⊕ x28 ⊕ t14 | 2 |
| 39 | y6 = x13 ⊕ t6 ⊕ x30 | 2 |
| 40 | y14 = x13 ⊕ x22 ⊕ t16 | 2 |
| 41 | y31 = x6 ⊕ x23 ⊕ t15 | 2 |
| 42 | y22 = x14 ⊕ x29 ⊕ t16 | 2 |
| 43 | y24 = x7 ⊕ x16 ⊕ t17 | 2 |
| 44 | y7 = x6 ⊕ x15 ⊕ t18 | 2 |
| 45 | y15 = x7 ⊕ x22 ⊕ t18 | 2 |
| 46 | y25 = t0 ⊕ x17 ⊕ t19 | 2 |
| 47 | y16 = x23 ⊕ x24 ⊕ t17 | 2 |
| 48 | y17 = x9 ⊕ t10 ⊕ t19 | 2 |
| 49 | y27 = t1 ⊕ x19 ⊕ t20 | 2 |
| 50 | y28 = x4 ⊕ t4 ⊕ t21 | 2 |
| 51 | y19 = x11 ⊕ t12 ⊕ t20 | 2 |
| 52 | y20 = t7 ⊕ x28 ⊕ t21 | 2 |
| 53 | y23 = x22 ⊕ t15 ⊕ x31 | 2 |

[Repository overview](../README.md)
