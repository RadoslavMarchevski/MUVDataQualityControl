# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

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
CMAKE_COMMAND = /afs/cern.ch/sw/lcg/contrib/CMake/3.2.3/Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /afs/cern.ch/sw/lcg/contrib/CMake/3.2.3/Linux-x86_64/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/analysis_new

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/analysis_new/build

# Include any dependencies generated for this target.
include CMakeFiles/Kmu2.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Kmu2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Kmu2.dir/flags.make

CMakeFiles/Kmu2.dir/main.cc.o: CMakeFiles/Kmu2.dir/flags.make
CMakeFiles/Kmu2.dir/main.cc.o: ../main.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/analysis_new/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/Kmu2.dir/main.cc.o"
	/afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Kmu2.dir/main.cc.o -c /afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/analysis_new/main.cc

CMakeFiles/Kmu2.dir/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Kmu2.dir/main.cc.i"
	/afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/analysis_new/main.cc > CMakeFiles/Kmu2.dir/main.cc.i

CMakeFiles/Kmu2.dir/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Kmu2.dir/main.cc.s"
	/afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/analysis_new/main.cc -o CMakeFiles/Kmu2.dir/main.cc.s

CMakeFiles/Kmu2.dir/main.cc.o.requires:
.PHONY : CMakeFiles/Kmu2.dir/main.cc.o.requires

CMakeFiles/Kmu2.dir/main.cc.o.provides: CMakeFiles/Kmu2.dir/main.cc.o.requires
	$(MAKE) -f CMakeFiles/Kmu2.dir/build.make CMakeFiles/Kmu2.dir/main.cc.o.provides.build
.PHONY : CMakeFiles/Kmu2.dir/main.cc.o.provides

CMakeFiles/Kmu2.dir/main.cc.o.provides.build: CMakeFiles/Kmu2.dir/main.cc.o

# Object files for target Kmu2
Kmu2_OBJECTS = \
"CMakeFiles/Kmu2.dir/main.cc.o"

# External object files for target Kmu2
Kmu2_EXTERNAL_OBJECTS =

../Kmu2: CMakeFiles/Kmu2.dir/main.cc.o
../Kmu2: CMakeFiles/Kmu2.dir/build.make
../Kmu2: ../lib/liblKmu2-static.a
../Kmu2: CMakeFiles/Kmu2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../Kmu2"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Kmu2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Kmu2.dir/build: ../Kmu2
.PHONY : CMakeFiles/Kmu2.dir/build

CMakeFiles/Kmu2.dir/requires: CMakeFiles/Kmu2.dir/main.cc.o.requires
.PHONY : CMakeFiles/Kmu2.dir/requires

CMakeFiles/Kmu2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Kmu2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Kmu2.dir/clean

CMakeFiles/Kmu2.dir/depend:
	cd /afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/analysis_new/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/analysis_new /afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/analysis_new /afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/analysis_new/build /afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/analysis_new/build /afs/cern.ch/user/r/rmarchev/work/New_NA62Fw/analysis_new/build/CMakeFiles/Kmu2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Kmu2.dir/depend
