# Dilawar's PySCF fork: geodesic molecular paths

**Start with the geometry-only butane explorer:** 101 interpolated geometries,
interactive molecular visualization, and bond-length, angle and dihedral tracking.
No electronic energy or gradient calculations are needed for this example.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Dilawars9/pyscf/blob/dilawar7/geodesic-interpolation/dilawar7/geodesic/notebooks/butane_path_explorer.ipynb)

[Worked examples and code guide](dilawar7/geodesic/README.md) ·
[Butane walkthrough](dilawar7/geodesic/BUTANE.md) ·
[Colab/Jupyter notebook](dilawar7/geodesic/notebooks/butane_path_explorer.ipynb) ·
[Validation](dilawar7/geodesic/VALIDATION.md)

**Looking for the measurement calculations?** They are in
[`pyscf/geomopt/path_analysis.py`](pyscf/geomopt/path_analysis.py):
`distances()`, `angles()`, `dihedrals()`, `unwrap_dihedrals()` and `aligned_rmsd()`.
The interpolation adapter is [`pyscf/geomopt/geodesic.py`](pyscf/geomopt/geodesic.py).

This is a development fork. The geodesic algorithm comes from
[geodesic-interpolate](https://github.com/virtualzx-nad/geodesic-interpolate);
our work connects it to PySCF and adds geometry analysis and visualization.
The upstream PySCF project information follows below.

---

<div align="left">
  <img src="https://github.com/pyscf/pyscf-doc/blob/master/logo/pyscf-logo.png" height="80px"/>
</div>

Python-based Simulations of Chemistry Framework
-----------------------------------------------
[![Build Status](https://github.com/pyscf/pyscf/workflows/CI/badge.svg)](https://github.com/pyscf/pyscf/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/pyscf/pyscf/branch/master/graph/badge.svg)](https://codecov.io/gh/pyscf/pyscf)

2026-11-05

* [Stable release 2.15.0](https://github.com/pyscf/pyscf/releases/tag/v2.15.0)
* [Changelog](../master/CHANGELOG)
* [Documentation](http://www.pyscf.org)
* [Installation](#installation)
* [Features](../master/FEATURES)
* [News](https://pyscf.org/news.html): **3rd PySCF Developers Meeting!**


# Installation

* Install stable release:

        pip install pyscf

* New features developed in recent years are available in the pyscf-forge package:

        pip install pyscf-forge

* Certain modules are maintained as extensions of PySCF, such as dispersion,
  dmrgscf, fciqmc, icmpspt, properties, semiempirical, shciscf ... (more on
  https://github.com/pyscf) can be installed using pip:

        pip install pyscf[all]

  An individual extension can be installed:

        pip install pyscf[dispersion]

* More details of custom installation can be found in
  [installation manual](http://pyscf.org/user/install.html#build-from-source)


# Citing PySCF

## Base PySCF
The following paper should be cited in publications utilizing the PySCF program package:

[The Python simulations of chemistry framework: 10 years of an open-source quantum chemistry project](https://doi.org/10.1063/5.0337441),
Qiming Sun, Matthew R Hermes, Xiaojie Wu, Huanchen Zhai, Xing Zhang, Abdelrahman M. Ahmed, Juan José Aucar, Oliver J. Backhouse, Samragni Banerjee, Peng Bao, Nikolay A. Bogdanov, Kyle Bystrom, Frédéric Chapoton, Ning-Yuan Chen, Ivan Yu. Chernyshov, Helen S. Clifford, Sander Cohen-Janes, Zhi-Hao Cui, Yann D. Damour, Nike Dattani, Linus Bjarne Dittmer, Sebastian Ehlert, Janus Juul Eriksen, Francesco A. Evangelista, Simon A. Ewing, Ardavan Farahvash, Kevin Focke, Yang Gao, Kevin E. Gasperich, Nathan Gillispie, Jonas Greiner, Matthew R. Hennefarth, Jan Hermann, Christopher Hillenbrand, Joonatan Huhtasalo, Basil Ibrahim, Bhavnesh Jangid, Alireza Nejati Javaremi, Andrew J. Jenkins, Yu Jin, Daniel S. King, Derk Pieter Kooi, Jo S. Kurian, Henrik R. Larsson, Bryan Tak Gwong Lau, Seunghoon Lee, Susi Lehtola, Chenghan Li, Hao Li, Jiachen Li, Rui Li, Shuhang Li, Aleksandr O. Lykhin, Ankit Mahajan, Nastasia Mauger, Pablo del Mazo-Sevillano, Jonathan Moussa, Kousuke Nakano, Verena A. Neufeld, Linqing Peng, Hung Q. Pham, Peter Pinski, Pavel Pokhilko, Zhichen Pu, Yubing Qian, Stephen Jon Quiton, Wanja T. Schulze, Thais R. Scott, Aniruddha Seal, James D. Serna, James E. T. Smith, Kori E. Smyser, Terrence Stahl, Chong Sun, Kevin J. Sung, Egor Trushin, Shiv Upadhyay, Ethan A. Vo, Thijs Vogels, Shirong Wang, Tai Wang, Xiao Wang, Xubo Wang, Yuanheng Wang, Mark Williamson, Junjie Yang, Hong-Zhou Ye, Chia-Nan Yeh, Haiyang Yu, Jincheng Yu, Victor Wen-zhe Yu, Chaoqun Zhang, Dayou Zhang, Yichi Zhang, Zijun Zhao, Zehao Zhou, Andrew J. Zhu, Tianyu Zhu, Timothy C. Berkelbach, Laura Gagliardi , Sandeep Sharma, Alexander Sokolov, Garnet Kin-Lic Chan,
*J. Chem. Phys.*, **165**, 102502 (2026). doi:[10.1063/5.0337441](https://doi.org/10.1063/5.0337441)

## Density functional calculations

As PySCF does not implement density functionals, instead employing external libraries to handle their evaluation, these libraries should also be cited in publications employing PySCF for density functional calculations.

If your calculation employed Libxc, cite

[Recent developments in libxc — A comprehensive library of functionals for density functional theory](https://doi.org/10.1016/j.softx.2017.11.002),
Susi Lehtola, Conrad Steigemann, Micael J.T. Oliveira, and Miguel A.L. Marques,
*SoftwareX* **7**, 1 (2018). doi:[10.1016/j.softx.2017.11.002](https://doi.org/10.1016/j.softx.2017.11.002)

If your calculation employed XCFun, cite

[Arbitrary-order density functional response theory from automatic differentiation](https://doi.org/10.1021/ct100117s),
Ulf Ekström, Lucas Visscher, Radovan Bast, Andreas J. Thorvaldsen, and Kenneth Ruud,
*J. Chem. Theory Comput.* **6**, 1971 (2010). doi:[10.1021/ct100117s](https://doi.org/10.1021/ct100117s)

# Bug reports and feature requests

Please submit tickets on the [issues](https://github.com/pyscf/pyscf/issues) page.

