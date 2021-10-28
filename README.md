# Microvascular Blood & Interstitial Flow Simulator

***A 2nd generation, user-friendly version of this code is currently in development.***

This C++ library has been used to simulate fluid transport in microvascular tissue [here](http://www.nature.com/articles/s41551-018-0306-y). Vascular networks are imported as weighted, undirected graphs which represent a network(s) of blood vessels generated either synthetically or by segmenting and skeletonising biomedical images. 

This microvascular flow solver is the 1st generation version of the [code](https://zenodo.org/record/1414160#.YXbN7y1Q1bV) which forms the basis of the REANIMATE (**Rea**listic **N**umerical **I**mage-based **M**odelling of Biologic**a**l **T**issue Substrat**e**s) framework published [here](http://www.nature.com/articles/s41551-018-0306-y) and corresponding mathematical methods [here](https://journals.plos.org/ploscompbiol/article/comments?id=10.1371/journal.pcbi.1006751).

The flow solver implements several mathematical models in research literature including:
* [Poiseuille's Law](https://www.annualreviews.org/doi/10.1146/annurev.fl.25.010193.000245) for steady-state,axisymmetric, laminar pipe flow in 1D networks with known boundary conditions.
* [Empirical blood viscosity laws](https://journals.physiology.org/doi/full/10.1152/ajpheart.00297.2005) to compute bulk blood viscosity as a function of vessel diameter and haematocrit, thereby capturing the Fåhraeus-Lindqvist effect.
* [Empirical phase separation law](https://www.ahajournals.org/doi/10.1161/01.res.67.4.826?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed) to calculate the disproportion distribution of haematocrit at microvascular bifurcations.
* [Blood flow estimation model](https://onlinelibrary.wiley.com/doi/10.1111/j.1549-8719.2012.00184.x) to simulate blood flow in a microvascular network with incomplete boundary conditions.
* [Interstitial flow model](https://journals.plos.org/ploscompbiol/article/comments?id=10.1371/journal.pcbi.1006751) which uses a Green's function method to simulate transvascular fluid transport.

If you wish to cite the published articles where this software has been used, please use the following references:
> [Computational fluid dynamics with imaging of cleared tissue and of in vivo perfusion predicts drug uptake and treatment responses in tumours](http://www.nature.com/articles/s41551-018-0306-y)<br>
> Angela d'Esposito & Paul W. Sweeney et al.

> [Modelling the transport of fluid through heterogeneous, whole tumours in silico](https://journals.plos.org/ploscompbiol/article/comments?id=10.1371/journal.pcbi.1006751)<br>
> Paul W. Sweeney et al.

If you wish to cite the software directly in your work, please use the following reference:
> [Vascular and interstitial flow solver for discrete microvascular networks](http://doi.org/10.5281/zenodo.1414160)<br>
> Paul W. Sweeney, Simon Walker-Samuel & Rebecca J. Shipley. 

## Installation
This software is compatible with C++11, and has been tested on Ubuntu 18.04 LTS and macOS Big Sur. 
Other distributions of Linux and Windows should work as well.

To install the scripts from source, download zip file on GitHub page or run the following in a terminal:

```bash
git clone https://github.com/psweens/VF_NatureBioEng.git
```

## Contributing
This C++ library is an open-source project started by [Dr Paul Sweeney](www.psweeney.co.uk) during his PhD at University College London under the supervision of [Prof. Rebecca Shipley](https://mecheng.ucl.ac.uk/people/profile/dr-rebecca-shipley/) and [Prof. Simon Walker-Samuel](http://simonwalkersamuel.com). 

This Github repository **is not maintained**. However, a 2nd generation, user-friendly version of this code is currently in development. Please contact me directly for further details.

## Acknowledgements
I would like to acknowledge that this C++ library utilises several portions of code originally written by [Prof. Timothy Secomb](https://github.com/secomb) but amended to allow integration of the [Armadillo C++ library](http://arma.sourceforge.net/).
