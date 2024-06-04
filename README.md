# CS219_Gpufit_ACC

## Installation ##
[Documents](https://gpufit.readthedocs.io/en/latest/installation.html#compiling-gpufit-on-linux)
```
git clone https://github.com/gpufit/Gpufit.git Gpufit
mkdir Gpufit-build
cd Gpufit-build
```
* [Bug-fix](https://github.com/gpufit/Gpufit/issues/129)
```
cd ../
cd Gpufit
git reset --hard 12b3cf4
```
* [Rerun](https://gpufit.readthedocs.io/en/latest/installation.html#compiling-gpufit-on-linux)
```
cd ../
cd Gpufit-build/
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_C_COMPILER=gcc-12 ../Gpufit
make
```
* Test (if runs then success) 
```
./Gpufit_Cpufit_Performance_Comparison
```

## Simple Usage Example ##
* Create example cpp file helloworld under Gpufit/examples/c++/
* Edit CMakeLists.txt
- At the bottom of the file, add “add_example(Gpufit hello_world) “
* Cmake and make like above procedure
