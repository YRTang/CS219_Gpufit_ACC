# CS219_Gpufit_ACC: GPU acceleration for 5G Cloud RAN 

Open RAN (O-RAN) provides an open RAN system for the 5G/4G infrastructures that are traditionally constructed with commodity and vendor products. It further splits the functions into non real-time and near real-time components. In 5G Cloud Radio Access Network (RAN), a major portion of computation workload has been shifted to the edge/cloud servers, where acceleration is needed for the near real-time components. In this project, you will be working on applying two acceleration approaches to one such example component, i.e., 5G channel estimator, as the case study.

## Goal
Implement two acceleration approaches for the given real-time module of channel estimator and compare their performance gains

## Tasks
- Preparation: wrap the estimator in a bare minimal server program which handles estimation requests (function input) and return the results (function output) sequentially. This works as the baseline.
    - Equation to solve (in latex): $y_k = \sum_{p=1}^Pe^{-j2\pi(m_k\tau_p-n_k\nu_p)}h_p$
    - Input: $m_k, n_k, y_k$ for $k=1,\cdots,K$, path number $P$ and pilot number $K$
    - Output: $\tau_p, \nu_p, h_p$ for $p=1,\cdots,P$
- **[Focus 1]** CPU acceleration: use thread pooling to parallelize concurrent estimator tasks
- **[Focus 2]** GPU acceleration: Integrate GPU acceleration with gpufit() library
    - Tutorials given at https://gpufit.readthedocs.io/en/latest/index.html
- Evaluation Metrics: processing time under a given estimator setting.

## Installation
### [Documents](https://gpufit.readthedocs.io/en/latest/installation.html#compiling-gpufit-on-linux)
```
git clone https://github.com/gpufit/Gpufit.git Gpufit
mkdir Gpufit-build
cd Gpufit-build
```
### [Bug-fix](https://github.com/gpufit/Gpufit/issues/129)
```
cd ../
cd Gpufit
git reset --hard 12b3cf4
```
### [Rerun](https://gpufit.readthedocs.io/en/latest/installation.html#compiling-gpufit-on-linux)
```
cd ../
cd Gpufit-build/
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_C_COMPILER=gcc-12 ../Gpufit
make
```
### Test (the original example) 
```
./Gpufit_Cpufit_Performance_Comparison
```

## Simple Usage Example
* Create example cpp file helloworld under Gpufit/examples/c++/
* Edit CMakeLists.txt
    - At the bottom of the file, add “add_example(Gpufit hello_world) “
* Cmake and make like above procedure
