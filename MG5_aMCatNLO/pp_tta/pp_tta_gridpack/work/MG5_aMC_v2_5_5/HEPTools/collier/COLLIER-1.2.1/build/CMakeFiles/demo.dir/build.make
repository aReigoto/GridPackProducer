# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /afs/cern.ch/user/a/areigoto/developing/FCC/Gridpack_Producer_Folders/GridPackProducer/MG5_aMCatNLO/pp_tta/pp_tta_gridpack/work/MG5_aMC_v2_5_5/HEPTools/collier/COLLIER-1.2.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /afs/cern.ch/user/a/areigoto/developing/FCC/Gridpack_Producer_Folders/GridPackProducer/MG5_aMCatNLO/pp_tta/pp_tta_gridpack/work/MG5_aMC_v2_5_5/HEPTools/collier/COLLIER-1.2.1/build

# Include any dependencies generated for this target.
include CMakeFiles/demo.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/demo.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/demo.dir/flags.make

CMakeFiles/demo.dir/demos/demo.f90.o: CMakeFiles/demo.dir/flags.make
CMakeFiles/demo.dir/demos/demo.f90.o: ../demos/demo.f90
	$(CMAKE_COMMAND) -E cmake_progress_report /afs/cern.ch/user/a/areigoto/developing/FCC/Gridpack_Producer_Folders/GridPackProducer/MG5_aMCatNLO/pp_tta/pp_tta_gridpack/work/MG5_aMC_v2_5_5/HEPTools/collier/COLLIER-1.2.1/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/demo.dir/demos/demo.f90.o"
	/cvmfs/sft.cern.ch/lcg/releases/LCG_88/gcc/4.9.1/x86_64-slc6/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /afs/cern.ch/user/a/areigoto/developing/FCC/Gridpack_Producer_Folders/GridPackProducer/MG5_aMCatNLO/pp_tta/pp_tta_gridpack/work/MG5_aMC_v2_5_5/HEPTools/collier/COLLIER-1.2.1/demos/demo.f90 -o CMakeFiles/demo.dir/demos/demo.f90.o

CMakeFiles/demo.dir/demos/demo.f90.o.requires:
.PHONY : CMakeFiles/demo.dir/demos/demo.f90.o.requires

CMakeFiles/demo.dir/demos/demo.f90.o.provides: CMakeFiles/demo.dir/demos/demo.f90.o.requires
	$(MAKE) -f CMakeFiles/demo.dir/build.make CMakeFiles/demo.dir/demos/demo.f90.o.provides.build
.PHONY : CMakeFiles/demo.dir/demos/demo.f90.o.provides

CMakeFiles/demo.dir/demos/demo.f90.o.provides.build: CMakeFiles/demo.dir/demos/demo.f90.o

# Object files for target demo
demo_OBJECTS = \
"CMakeFiles/demo.dir/demos/demo.f90.o"

# External object files for target demo
demo_EXTERNAL_OBJECTS =

../demos/demo: CMakeFiles/demo.dir/demos/demo.f90.o
../demos/demo: CMakeFiles/demo.dir/build.make
../demos/demo: ../libcollier.a
../demos/demo: CMakeFiles/demo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking Fortran executable ../demos/demo"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/demo.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/demo.dir/build: ../demos/demo
.PHONY : CMakeFiles/demo.dir/build

CMakeFiles/demo.dir/requires: CMakeFiles/demo.dir/demos/demo.f90.o.requires
.PHONY : CMakeFiles/demo.dir/requires

CMakeFiles/demo.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/demo.dir/cmake_clean.cmake
.PHONY : CMakeFiles/demo.dir/clean

CMakeFiles/demo.dir/depend:
	cd /afs/cern.ch/user/a/areigoto/developing/FCC/Gridpack_Producer_Folders/GridPackProducer/MG5_aMCatNLO/pp_tta/pp_tta_gridpack/work/MG5_aMC_v2_5_5/HEPTools/collier/COLLIER-1.2.1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /afs/cern.ch/user/a/areigoto/developing/FCC/Gridpack_Producer_Folders/GridPackProducer/MG5_aMCatNLO/pp_tta/pp_tta_gridpack/work/MG5_aMC_v2_5_5/HEPTools/collier/COLLIER-1.2.1 /afs/cern.ch/user/a/areigoto/developing/FCC/Gridpack_Producer_Folders/GridPackProducer/MG5_aMCatNLO/pp_tta/pp_tta_gridpack/work/MG5_aMC_v2_5_5/HEPTools/collier/COLLIER-1.2.1 /afs/cern.ch/user/a/areigoto/developing/FCC/Gridpack_Producer_Folders/GridPackProducer/MG5_aMCatNLO/pp_tta/pp_tta_gridpack/work/MG5_aMC_v2_5_5/HEPTools/collier/COLLIER-1.2.1/build /afs/cern.ch/user/a/areigoto/developing/FCC/Gridpack_Producer_Folders/GridPackProducer/MG5_aMCatNLO/pp_tta/pp_tta_gridpack/work/MG5_aMC_v2_5_5/HEPTools/collier/COLLIER-1.2.1/build /afs/cern.ch/user/a/areigoto/developing/FCC/Gridpack_Producer_Folders/GridPackProducer/MG5_aMCatNLO/pp_tta/pp_tta_gridpack/work/MG5_aMC_v2_5_5/HEPTools/collier/COLLIER-1.2.1/build/CMakeFiles/demo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/demo.dir/depend

