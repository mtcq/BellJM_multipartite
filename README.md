# BellJM_multipartite

## Code to accompany: Code to accompany: [All incompatible measurements on qubits lead to multiparticle Bell nonlocality](https://arxiv.org/abs/2403.10564)

#### Martin Plávala, Otfried Gühne, Marco Túlio Quintino


This is a repository for the article [All incompatible measurements on qubits lead to multiparticle Bell nonlocality](https://arxiv.org/abs/2403.10564).

This repository consists of:

- [findJM3XZ.m](https://github.com/mtcq/BellJM_multipartite/blob/main/findJM3XZ.m):
MATLAB script which shows that the critical visibility for a pair of Pauli measurements to become JM3 is between eta = 0.7937 and eta = 0.7938.

- [findJM2T.m](https://github.com/mtcq/BellJM_multipartite/blob/main/findJM2T.m):
MATLAB script which shows that the critical visibility for the T measurements with 2 parties is \eta_2^{Bell} = sqrt(2/sqrt(7)) \approx 0.8694

- [findJM3T_GHZ_Y.m](https://github.com/mtcq/BellJM_multipartite/blob/main/findJM3T_GHZ_Y.m):
MATLAB script which shows that the critical visibility for the T measurements with 3 parties respects \eta_3^{Bell} < 0.8007

- [L22_22_22](https://github.com/mtcq/BellJM_multipartite/tree/main/L22_22_22):
A folder which contains all the tight Bell inequalities from the tripartite scenario where each party has access to two inputs and each input has two outcomes. These inequalities are known as "Sliwa inequalities" and were first published at [https://arxiv.org/abs/quant-ph/0305190](https://arxiv.org/abs/quant-ph/0305190).
Our database was obtained from [faacets](https://github.com/denisrosset/faacets-data/tree/master/solved/L22_22_22), an online database for Bell inequalities.

- [faacets_converter.sh](https://github.com/mtcq/BellJM_multipartite/blob/main/L22_22_22/faacets_converter.sh):
  A bash script which converts the inequalities from faacets from yaml format to MATLAB format. 
