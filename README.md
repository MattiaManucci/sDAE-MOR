[![arXiv][arxiv-shield]][arxiv-url]
[![DOI][doi-shield]][doi-url]

# [sDAE-MOR][arxiv-url]
This repository contains all scripts required to reproduce the figures and tables from the paper [Balancing-based model reduction for switched descriptor systems][arxiv-url]. The code implements the proposed methodology, runs the numerical experiments, and ensures full reproducibility of the results reported in the manuscript. 


##Code info:

* **FigX**: the script generates Figure X ($X\in{1,2,3,4}$) of the work [Balancing-based model reduction for switched descriptor systems][arxiv-url].

Functions:

* **Solve\_LS_GLE_2**: solve Generalized Lyapunov Equations with stationary iteration algorithm, see [preprint][arxiv-url]----> $AX+X A^T +\Sigma_{j=1}^{M} (D_j X D_{j}^{T}+B_j B_{j}^{T})=0$ 
* **solve_KS**: solve large-scale Lyapunov equation of type $AX+XA^T+BB^T=0$ using a standard polynomial Krylov method; 
* **solve\_KS_t**: solve large-scale Lyapunov equation of type $A^T X+XA+C^T C=0$ using a standard polynomial Krylov method.


## Citing
If you use this project for academic work, please consider citing our
[publication][arxiv-url]:

    M. Manucci and B. Unger
    Balancing-based model reduction for switched descriptor systems
    ArXiv e-print **number**, 2026.
    
## License
Distributed under the MIT License. See `LICENSE` for more information.


## Contacts

* Mattia Manucci - [mattia.manucci@simtech.uni-stuttgart.de](mattia.manucci@simtech.uni-stuttgart.de)
* Benjamin Unger - [benjamin.unger@kit.edu](benjamin.unger@kit.edu)



[doi-shield]: https://img.shields.io/badge/DOI-10.5281%20%2F%20zenodo.8335231-blue.svg?style=for-the-badge
[doi-url]: ???
[arxiv-shield]: https://img.shields.io/badge/arXiv-2204.13474-b31b1b.svg?style=for-the-badge
[arxiv-url]:???







