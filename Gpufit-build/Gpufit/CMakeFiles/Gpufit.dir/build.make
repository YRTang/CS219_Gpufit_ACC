# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/test/sarahTang/Gpufit

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/test/sarahTang/Gpufit-build

# Include any dependencies generated for this target.
include Gpufit/CMakeFiles/Gpufit.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include Gpufit/CMakeFiles/Gpufit.dir/compiler_depend.make

# Include the progress variables for this target.
include Gpufit/CMakeFiles/Gpufit.dir/progress.make

# Include the compile flags for this target's objects.
include Gpufit/CMakeFiles/Gpufit.dir/flags.make

Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_lm_fit_cuda.cu.o: /home/test/sarahTang/Gpufit/Gpufit/lm_fit_cuda.cu
Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_lm_fit_cuda.cu.o: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_lm_fit_cuda.cu.o.depend
Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_lm_fit_cuda.cu.o: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_lm_fit_cuda.cu.o.RELEASE.cmake
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold --progress-dir=/home/test/sarahTang/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building NVCC (Device) object Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_lm_fit_cuda.cu.o"
	cd /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir && /usr/bin/cmake -E make_directory /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//.
	cd /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir && /usr/bin/cmake -D verbose:BOOL=$(VERBOSE) -D build_configuration:STRING=RELEASE -D generated_file:STRING=/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//./Gpufit_generated_lm_fit_cuda.cu.o -D generated_cubin_file:STRING=/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//./Gpufit_generated_lm_fit_cuda.cu.o.cubin.txt -P /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//Gpufit_generated_lm_fit_cuda.cu.o.RELEASE.cmake

Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_kernels.cu.o: /home/test/sarahTang/Gpufit/Gpufit/cuda_kernels.cu
Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_kernels.cu.o: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_kernels.cu.o.depend
Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_kernels.cu.o: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_kernels.cu.o.RELEASE.cmake
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold --progress-dir=/home/test/sarahTang/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building NVCC (Device) object Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_kernels.cu.o"
	cd /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir && /usr/bin/cmake -E make_directory /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//.
	cd /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir && /usr/bin/cmake -D verbose:BOOL=$(VERBOSE) -D build_configuration:STRING=RELEASE -D generated_file:STRING=/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//./Gpufit_generated_cuda_kernels.cu.o -D generated_cubin_file:STRING=/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//./Gpufit_generated_cuda_kernels.cu.o.cubin.txt -P /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//Gpufit_generated_cuda_kernels.cu.o.RELEASE.cmake

Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_info.cu.o: /home/test/sarahTang/Gpufit/Gpufit/info.cu
Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_info.cu.o: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_info.cu.o.depend
Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_info.cu.o: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_info.cu.o.RELEASE.cmake
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold --progress-dir=/home/test/sarahTang/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building NVCC (Device) object Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_info.cu.o"
	cd /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir && /usr/bin/cmake -E make_directory /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//.
	cd /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir && /usr/bin/cmake -D verbose:BOOL=$(VERBOSE) -D build_configuration:STRING=RELEASE -D generated_file:STRING=/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//./Gpufit_generated_info.cu.o -D generated_cubin_file:STRING=/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//./Gpufit_generated_info.cu.o.cubin.txt -P /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//Gpufit_generated_info.cu.o.RELEASE.cmake

Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_gpu_data.cu.o: /home/test/sarahTang/Gpufit/Gpufit/gpu_data.cu
Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_gpu_data.cu.o: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_gpu_data.cu.o.depend
Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_gpu_data.cu.o: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_gpu_data.cu.o.RELEASE.cmake
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold --progress-dir=/home/test/sarahTang/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building NVCC (Device) object Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_gpu_data.cu.o"
	cd /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir && /usr/bin/cmake -E make_directory /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//.
	cd /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir && /usr/bin/cmake -D verbose:BOOL=$(VERBOSE) -D build_configuration:STRING=RELEASE -D generated_file:STRING=/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//./Gpufit_generated_gpu_data.cu.o -D generated_cubin_file:STRING=/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//./Gpufit_generated_gpu_data.cu.o.cubin.txt -P /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//Gpufit_generated_gpu_data.cu.o.RELEASE.cmake

Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_gaussjordan.cu.o: /home/test/sarahTang/Gpufit/Gpufit/cuda_gaussjordan.cu
Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_gaussjordan.cu.o: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_gaussjordan.cu.o.depend
Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_gaussjordan.cu.o: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_gaussjordan.cu.o.RELEASE.cmake
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold --progress-dir=/home/test/sarahTang/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building NVCC (Device) object Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_gaussjordan.cu.o"
	cd /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir && /usr/bin/cmake -E make_directory /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//.
	cd /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir && /usr/bin/cmake -D verbose:BOOL=$(VERBOSE) -D build_configuration:STRING=RELEASE -D generated_file:STRING=/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//./Gpufit_generated_cuda_gaussjordan.cu.o -D generated_cubin_file:STRING=/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//./Gpufit_generated_cuda_gaussjordan.cu.o.cubin.txt -P /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir//Gpufit_generated_cuda_gaussjordan.cu.o.RELEASE.cmake

Gpufit/CMakeFiles/Gpufit.dir/gpufit.cpp.o: Gpufit/CMakeFiles/Gpufit.dir/flags.make
Gpufit/CMakeFiles/Gpufit.dir/gpufit.cpp.o: /home/test/sarahTang/Gpufit/Gpufit/gpufit.cpp
Gpufit/CMakeFiles/Gpufit.dir/gpufit.cpp.o: Gpufit/CMakeFiles/Gpufit.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/test/sarahTang/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object Gpufit/CMakeFiles/Gpufit.dir/gpufit.cpp.o"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Gpufit/CMakeFiles/Gpufit.dir/gpufit.cpp.o -MF CMakeFiles/Gpufit.dir/gpufit.cpp.o.d -o CMakeFiles/Gpufit.dir/gpufit.cpp.o -c /home/test/sarahTang/Gpufit/Gpufit/gpufit.cpp

Gpufit/CMakeFiles/Gpufit.dir/gpufit.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Gpufit.dir/gpufit.cpp.i"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/test/sarahTang/Gpufit/Gpufit/gpufit.cpp > CMakeFiles/Gpufit.dir/gpufit.cpp.i

Gpufit/CMakeFiles/Gpufit.dir/gpufit.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Gpufit.dir/gpufit.cpp.s"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/test/sarahTang/Gpufit/Gpufit/gpufit.cpp -o CMakeFiles/Gpufit.dir/gpufit.cpp.s

Gpufit/CMakeFiles/Gpufit.dir/info.cpp.o: Gpufit/CMakeFiles/Gpufit.dir/flags.make
Gpufit/CMakeFiles/Gpufit.dir/info.cpp.o: /home/test/sarahTang/Gpufit/Gpufit/info.cpp
Gpufit/CMakeFiles/Gpufit.dir/info.cpp.o: Gpufit/CMakeFiles/Gpufit.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/test/sarahTang/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object Gpufit/CMakeFiles/Gpufit.dir/info.cpp.o"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Gpufit/CMakeFiles/Gpufit.dir/info.cpp.o -MF CMakeFiles/Gpufit.dir/info.cpp.o.d -o CMakeFiles/Gpufit.dir/info.cpp.o -c /home/test/sarahTang/Gpufit/Gpufit/info.cpp

Gpufit/CMakeFiles/Gpufit.dir/info.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Gpufit.dir/info.cpp.i"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/test/sarahTang/Gpufit/Gpufit/info.cpp > CMakeFiles/Gpufit.dir/info.cpp.i

Gpufit/CMakeFiles/Gpufit.dir/info.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Gpufit.dir/info.cpp.s"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/test/sarahTang/Gpufit/Gpufit/info.cpp -o CMakeFiles/Gpufit.dir/info.cpp.s

Gpufit/CMakeFiles/Gpufit.dir/lm_fit.cpp.o: Gpufit/CMakeFiles/Gpufit.dir/flags.make
Gpufit/CMakeFiles/Gpufit.dir/lm_fit.cpp.o: /home/test/sarahTang/Gpufit/Gpufit/lm_fit.cpp
Gpufit/CMakeFiles/Gpufit.dir/lm_fit.cpp.o: Gpufit/CMakeFiles/Gpufit.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/test/sarahTang/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object Gpufit/CMakeFiles/Gpufit.dir/lm_fit.cpp.o"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Gpufit/CMakeFiles/Gpufit.dir/lm_fit.cpp.o -MF CMakeFiles/Gpufit.dir/lm_fit.cpp.o.d -o CMakeFiles/Gpufit.dir/lm_fit.cpp.o -c /home/test/sarahTang/Gpufit/Gpufit/lm_fit.cpp

Gpufit/CMakeFiles/Gpufit.dir/lm_fit.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Gpufit.dir/lm_fit.cpp.i"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/test/sarahTang/Gpufit/Gpufit/lm_fit.cpp > CMakeFiles/Gpufit.dir/lm_fit.cpp.i

Gpufit/CMakeFiles/Gpufit.dir/lm_fit.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Gpufit.dir/lm_fit.cpp.s"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/test/sarahTang/Gpufit/Gpufit/lm_fit.cpp -o CMakeFiles/Gpufit.dir/lm_fit.cpp.s

Gpufit/CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.o: Gpufit/CMakeFiles/Gpufit.dir/flags.make
Gpufit/CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.o: /home/test/sarahTang/Gpufit/Gpufit/lm_fit_cuda.cpp
Gpufit/CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.o: Gpufit/CMakeFiles/Gpufit.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/test/sarahTang/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object Gpufit/CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.o"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Gpufit/CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.o -MF CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.o.d -o CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.o -c /home/test/sarahTang/Gpufit/Gpufit/lm_fit_cuda.cpp

Gpufit/CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.i"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/test/sarahTang/Gpufit/Gpufit/lm_fit_cuda.cpp > CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.i

Gpufit/CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.s"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/test/sarahTang/Gpufit/Gpufit/lm_fit_cuda.cpp -o CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.s

Gpufit/CMakeFiles/Gpufit.dir/interface.cpp.o: Gpufit/CMakeFiles/Gpufit.dir/flags.make
Gpufit/CMakeFiles/Gpufit.dir/interface.cpp.o: /home/test/sarahTang/Gpufit/Gpufit/interface.cpp
Gpufit/CMakeFiles/Gpufit.dir/interface.cpp.o: Gpufit/CMakeFiles/Gpufit.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/test/sarahTang/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object Gpufit/CMakeFiles/Gpufit.dir/interface.cpp.o"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Gpufit/CMakeFiles/Gpufit.dir/interface.cpp.o -MF CMakeFiles/Gpufit.dir/interface.cpp.o.d -o CMakeFiles/Gpufit.dir/interface.cpp.o -c /home/test/sarahTang/Gpufit/Gpufit/interface.cpp

Gpufit/CMakeFiles/Gpufit.dir/interface.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Gpufit.dir/interface.cpp.i"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/test/sarahTang/Gpufit/Gpufit/interface.cpp > CMakeFiles/Gpufit.dir/interface.cpp.i

Gpufit/CMakeFiles/Gpufit.dir/interface.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Gpufit.dir/interface.cpp.s"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/test/sarahTang/Gpufit/Gpufit/interface.cpp -o CMakeFiles/Gpufit.dir/interface.cpp.s

# Object files for target Gpufit
Gpufit_OBJECTS = \
"CMakeFiles/Gpufit.dir/gpufit.cpp.o" \
"CMakeFiles/Gpufit.dir/info.cpp.o" \
"CMakeFiles/Gpufit.dir/lm_fit.cpp.o" \
"CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.o" \
"CMakeFiles/Gpufit.dir/interface.cpp.o"

# External object files for target Gpufit
Gpufit_EXTERNAL_OBJECTS = \
"/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_lm_fit_cuda.cu.o" \
"/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_kernels.cu.o" \
"/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_info.cu.o" \
"/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_gpu_data.cu.o" \
"/home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_gaussjordan.cu.o"

Gpufit/libGpufit.so: Gpufit/CMakeFiles/Gpufit.dir/gpufit.cpp.o
Gpufit/libGpufit.so: Gpufit/CMakeFiles/Gpufit.dir/info.cpp.o
Gpufit/libGpufit.so: Gpufit/CMakeFiles/Gpufit.dir/lm_fit.cpp.o
Gpufit/libGpufit.so: Gpufit/CMakeFiles/Gpufit.dir/lm_fit_cuda.cpp.o
Gpufit/libGpufit.so: Gpufit/CMakeFiles/Gpufit.dir/interface.cpp.o
Gpufit/libGpufit.so: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_lm_fit_cuda.cu.o
Gpufit/libGpufit.so: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_kernels.cu.o
Gpufit/libGpufit.so: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_info.cu.o
Gpufit/libGpufit.so: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_gpu_data.cu.o
Gpufit/libGpufit.so: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_gaussjordan.cu.o
Gpufit/libGpufit.so: Gpufit/CMakeFiles/Gpufit.dir/build.make
Gpufit/libGpufit.so: /usr/lib/x86_64-linux-gnu/libcudart_static.a
Gpufit/libGpufit.so: /usr/lib/x86_64-linux-gnu/librt.a
Gpufit/libGpufit.so: /home/test/sarahTang/Gpufit/Gpufit/Gpufit.def
Gpufit/libGpufit.so: Gpufit/CMakeFiles/Gpufit.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/test/sarahTang/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Linking CXX shared library libGpufit.so"
	cd /home/test/sarahTang/Gpufit-build/Gpufit && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Gpufit.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Gpufit/CMakeFiles/Gpufit.dir/build: Gpufit/libGpufit.so
.PHONY : Gpufit/CMakeFiles/Gpufit.dir/build

Gpufit/CMakeFiles/Gpufit.dir/clean:
	cd /home/test/sarahTang/Gpufit-build/Gpufit && $(CMAKE_COMMAND) -P CMakeFiles/Gpufit.dir/cmake_clean.cmake
.PHONY : Gpufit/CMakeFiles/Gpufit.dir/clean

Gpufit/CMakeFiles/Gpufit.dir/depend: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_gaussjordan.cu.o
Gpufit/CMakeFiles/Gpufit.dir/depend: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_cuda_kernels.cu.o
Gpufit/CMakeFiles/Gpufit.dir/depend: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_gpu_data.cu.o
Gpufit/CMakeFiles/Gpufit.dir/depend: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_info.cu.o
Gpufit/CMakeFiles/Gpufit.dir/depend: Gpufit/CMakeFiles/Gpufit.dir/Gpufit_generated_lm_fit_cuda.cu.o
	cd /home/test/sarahTang/Gpufit-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/test/sarahTang/Gpufit /home/test/sarahTang/Gpufit/Gpufit /home/test/sarahTang/Gpufit-build /home/test/sarahTang/Gpufit-build/Gpufit /home/test/sarahTang/Gpufit-build/Gpufit/CMakeFiles/Gpufit.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : Gpufit/CMakeFiles/Gpufit.dir/depend

