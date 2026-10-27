# Source Code and Supplementary Experimental Material

This repository contains source code and supplementary experimental material accompanying the submitted manuscript.

## Source code

The source files are located in the repository root:

- `candidate_nodes.py`
- `milp_model.py`
- `run_compare.py`

## Supplementary experimental material

The following material is provided separately from the manuscript because of space limitations. The Markdown files can be read directly in the repository.

| Item | Material | Contents |
| --- | --- | --- |
| S1 | [Depth-three MDS results](supplementary/mds_depth3.md) | Results for all 5,550 MDS matrices from Li et al., summarized in 24 Hamming-weight groups, including HILL comparisons and best verified circuit costs. |
| S2 | [AES MixColumns parameter sensitivity](supplementary/aes_depth3_sensitivity.md) | Depth-three results for 25 combinations of randomized seed-trace count R and per-target expansion limit W, including cost, solver status, gate counts, and elapsed time. |
| S3 | [Depth-two AES MixColumns circuit](supplementary/aes_depth2_circuit.md) | Complete sequence of 54 XOR3 operations, with gate depths and a reported cost of 252.18 GE under the UMC 180 nm cost model. |

## Interpretation

- GE denotes gate-equivalent area; n2 and n3 denote XOR2 and XOR3 gate counts.
- OPT means that the MILP proved optimality within the generated candidate space. It does not imply global circuit optimality.
- TL denotes a verified incumbent at the time limit.
- References to Li et al., HILL, and the cost model follow the accompanying manuscript.
