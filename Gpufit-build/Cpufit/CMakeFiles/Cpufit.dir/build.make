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
CMAKE_SOURCE_DIR = /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build

# Include any dependencies generated for this target.
include Cpufit/CMakeFiles/Cpufit.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include Cpufit/CMakeFiles/Cpufit.dir/compiler_depend.make

# Include the progress variables for this target.
include Cpufit/CMakeFiles/Cpufit.dir/progress.make

# Include the compile flags for this target's objects.
include Cpufit/CMakeFiles/Cpufit.dir/flags.make

Cpufit/CMakeFiles/Cpufit.dir/cpufit.cpp.o: Cpufit/CMakeFiles/Cpufit.dir/flags.make
Cpufit/CMakeFiles/Cpufit.dir/cpufit.cpp.o: /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/cpufit.cpp
Cpufit/CMakeFiles/Cpufit.dir/cpufit.cpp.o: Cpufit/CMakeFiles/Cpufit.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Cpufit/CMakeFiles/Cpufit.dir/cpufit.cpp.o"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Cpufit/CMakeFiles/Cpufit.dir/cpufit.cpp.o -MF CMakeFiles/Cpufit.dir/cpufit.cpp.o.d -o CMakeFiles/Cpufit.dir/cpufit.cpp.o -c /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/cpufit.cpp

Cpufit/CMakeFiles/Cpufit.dir/cpufit.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Cpufit.dir/cpufit.cpp.i"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/cpufit.cpp > CMakeFiles/Cpufit.dir/cpufit.cpp.i

Cpufit/CMakeFiles/Cpufit.dir/cpufit.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Cpufit.dir/cpufit.cpp.s"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/cpufit.cpp -o CMakeFiles/Cpufit.dir/cpufit.cpp.s

Cpufit/CMakeFiles/Cpufit.dir/info.cpp.o: Cpufit/CMakeFiles/Cpufit.dir/flags.make
Cpufit/CMakeFiles/Cpufit.dir/info.cpp.o: /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/info.cpp
Cpufit/CMakeFiles/Cpufit.dir/info.cpp.o: Cpufit/CMakeFiles/Cpufit.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object Cpufit/CMakeFiles/Cpufit.dir/info.cpp.o"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Cpufit/CMakeFiles/Cpufit.dir/info.cpp.o -MF CMakeFiles/Cpufit.dir/info.cpp.o.d -o CMakeFiles/Cpufit.dir/info.cpp.o -c /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/info.cpp

Cpufit/CMakeFiles/Cpufit.dir/info.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Cpufit.dir/info.cpp.i"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/info.cpp > CMakeFiles/Cpufit.dir/info.cpp.i

Cpufit/CMakeFiles/Cpufit.dir/info.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Cpufit.dir/info.cpp.s"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/info.cpp -o CMakeFiles/Cpufit.dir/info.cpp.s

Cpufit/CMakeFiles/Cpufit.dir/lm_fit.cpp.o: Cpufit/CMakeFiles/Cpufit.dir/flags.make
Cpufit/CMakeFiles/Cpufit.dir/lm_fit.cpp.o: /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/lm_fit.cpp
Cpufit/CMakeFiles/Cpufit.dir/lm_fit.cpp.o: Cpufit/CMakeFiles/Cpufit.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object Cpufit/CMakeFiles/Cpufit.dir/lm_fit.cpp.o"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Cpufit/CMakeFiles/Cpufit.dir/lm_fit.cpp.o -MF CMakeFiles/Cpufit.dir/lm_fit.cpp.o.d -o CMakeFiles/Cpufit.dir/lm_fit.cpp.o -c /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/lm_fit.cpp

Cpufit/CMakeFiles/Cpufit.dir/lm_fit.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Cpufit.dir/lm_fit.cpp.i"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/lm_fit.cpp > CMakeFiles/Cpufit.dir/lm_fit.cpp.i

Cpufit/CMakeFiles/Cpufit.dir/lm_fit.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Cpufit.dir/lm_fit.cpp.s"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/lm_fit.cpp -o CMakeFiles/Cpufit.dir/lm_fit.cpp.s

Cpufit/CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.o: Cpufit/CMakeFiles/Cpufit.dir/flags.make
Cpufit/CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.o: /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/lm_fit_cpp.cpp
Cpufit/CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.o: Cpufit/CMakeFiles/Cpufit.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object Cpufit/CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.o"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Cpufit/CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.o -MF CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.o.d -o CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.o -c /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/lm_fit_cpp.cpp

Cpufit/CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.i"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/lm_fit_cpp.cpp > CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.i

Cpufit/CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.s"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/lm_fit_cpp.cpp -o CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.s

Cpufit/CMakeFiles/Cpufit.dir/interface.cpp.o: Cpufit/CMakeFiles/Cpufit.dir/flags.make
Cpufit/CMakeFiles/Cpufit.dir/interface.cpp.o: /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/interface.cpp
Cpufit/CMakeFiles/Cpufit.dir/interface.cpp.o: Cpufit/CMakeFiles/Cpufit.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object Cpufit/CMakeFiles/Cpufit.dir/interface.cpp.o"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Cpufit/CMakeFiles/Cpufit.dir/interface.cpp.o -MF CMakeFiles/Cpufit.dir/interface.cpp.o.d -o CMakeFiles/Cpufit.dir/interface.cpp.o -c /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/interface.cpp

Cpufit/CMakeFiles/Cpufit.dir/interface.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Cpufit.dir/interface.cpp.i"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/interface.cpp > CMakeFiles/Cpufit.dir/interface.cpp.i

Cpufit/CMakeFiles/Cpufit.dir/interface.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Cpufit.dir/interface.cpp.s"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/interface.cpp -o CMakeFiles/Cpufit.dir/interface.cpp.s

# Object files for target Cpufit
Cpufit_OBJECTS = \
"CMakeFiles/Cpufit.dir/cpufit.cpp.o" \
"CMakeFiles/Cpufit.dir/info.cpp.o" \
"CMakeFiles/Cpufit.dir/lm_fit.cpp.o" \
"CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.o" \
"CMakeFiles/Cpufit.dir/interface.cpp.o"

# External object files for target Cpufit
Cpufit_EXTERNAL_OBJECTS =

Cpufit/libCpufit.so: Cpufit/CMakeFiles/Cpufit.dir/cpufit.cpp.o
Cpufit/libCpufit.so: Cpufit/CMakeFiles/Cpufit.dir/info.cpp.o
Cpufit/libCpufit.so: Cpufit/CMakeFiles/Cpufit.dir/lm_fit.cpp.o
Cpufit/libCpufit.so: Cpufit/CMakeFiles/Cpufit.dir/lm_fit_cpp.cpp.o
Cpufit/libCpufit.so: Cpufit/CMakeFiles/Cpufit.dir/interface.cpp.o
Cpufit/libCpufit.so: Cpufit/CMakeFiles/Cpufit.dir/build.make
Cpufit/libCpufit.so: /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit/Cpufit.def
Cpufit/libCpufit.so: Cpufit/CMakeFiles/Cpufit.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX shared library libCpufit.so"
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Cpufit.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Cpufit/CMakeFiles/Cpufit.dir/build: Cpufit/libCpufit.so
.PHONY : Cpufit/CMakeFiles/Cpufit.dir/build

Cpufit/CMakeFiles/Cpufit.dir/clean:
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit && $(CMAKE_COMMAND) -P CMakeFiles/Cpufit.dir/cmake_clean.cmake
.PHONY : Cpufit/CMakeFiles/Cpufit.dir/clean

Cpufit/CMakeFiles/Cpufit.dir/depend:
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/Cpufit /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/Cpufit/CMakeFiles/Cpufit.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : Cpufit/CMakeFiles/Cpufit.dir/depend

