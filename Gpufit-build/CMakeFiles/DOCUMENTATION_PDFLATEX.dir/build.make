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

# Utility rule file for DOCUMENTATION_PDFLATEX.

# Include any custom commands dependencies for this target.
include CMakeFiles/DOCUMENTATION_PDFLATEX.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/DOCUMENTATION_PDFLATEX.dir/progress.make

CMakeFiles/DOCUMENTATION_PDFLATEX:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold --progress-dir=/home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Converting documentation to PDF"
	/usr/bin/cmake -E chdir /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/docs/_build/latex /usr/bin/pdflatex -interaction=nonstopmode Gpufit.tex
	/usr/bin/cmake -E chdir /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit/docs/_build/latex /usr/bin/pdflatex -interaction=nonstopmode Gpufit.tex

DOCUMENTATION_PDFLATEX: CMakeFiles/DOCUMENTATION_PDFLATEX
DOCUMENTATION_PDFLATEX: CMakeFiles/DOCUMENTATION_PDFLATEX.dir/build.make
.PHONY : DOCUMENTATION_PDFLATEX

# Rule to build all files generated by this target.
CMakeFiles/DOCUMENTATION_PDFLATEX.dir/build: DOCUMENTATION_PDFLATEX
.PHONY : CMakeFiles/DOCUMENTATION_PDFLATEX.dir/build

CMakeFiles/DOCUMENTATION_PDFLATEX.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/DOCUMENTATION_PDFLATEX.dir/cmake_clean.cmake
.PHONY : CMakeFiles/DOCUMENTATION_PDFLATEX.dir/clean

CMakeFiles/DOCUMENTATION_PDFLATEX.dir/depend:
	cd /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build /home/test/kexinZheng/CS219_Gpufit_ACC/Gpufit-build/CMakeFiles/DOCUMENTATION_PDFLATEX.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/DOCUMENTATION_PDFLATEX.dir/depend

