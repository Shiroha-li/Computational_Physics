# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_SOURCE_DIR = /home/shiroha/Workspace/cp/Ising_Spin

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/shiroha/Workspace/cp/Ising_Spin/build

# Include any dependencies generated for this target.
include CMakeFiles/UnitTests2_catch2.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/UnitTests2_catch2.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/UnitTests2_catch2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/UnitTests2_catch2.dir/flags.make

CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.o: CMakeFiles/UnitTests2_catch2.dir/flags.make
CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.o: ../UnitTests_IsingSystem.cpp
CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.o: CMakeFiles/UnitTests2_catch2.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/shiroha/Workspace/cp/Ising_Spin/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.o -MF CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.o.d -o CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.o -c /home/shiroha/Workspace/cp/Ising_Spin/UnitTests_IsingSystem.cpp

CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/shiroha/Workspace/cp/Ising_Spin/UnitTests_IsingSystem.cpp > CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.i

CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/shiroha/Workspace/cp/Ising_Spin/UnitTests_IsingSystem.cpp -o CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.s

# Object files for target UnitTests2_catch2
UnitTests2_catch2_OBJECTS = \
"CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.o"

# External object files for target UnitTests2_catch2
UnitTests2_catch2_EXTERNAL_OBJECTS =

UnitTests2_catch2: CMakeFiles/UnitTests2_catch2.dir/UnitTests_IsingSystem.cpp.o
UnitTests2_catch2: CMakeFiles/UnitTests2_catch2.dir/build.make
UnitTests2_catch2: _deps/catch2-build/src/libCatch2Main.a
UnitTests2_catch2: _deps/catch2-build/src/libCatch2.a
UnitTests2_catch2: CMakeFiles/UnitTests2_catch2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/shiroha/Workspace/cp/Ising_Spin/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable UnitTests2_catch2"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/UnitTests2_catch2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/UnitTests2_catch2.dir/build: UnitTests2_catch2
.PHONY : CMakeFiles/UnitTests2_catch2.dir/build

CMakeFiles/UnitTests2_catch2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/UnitTests2_catch2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/UnitTests2_catch2.dir/clean

CMakeFiles/UnitTests2_catch2.dir/depend:
	cd /home/shiroha/Workspace/cp/Ising_Spin/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shiroha/Workspace/cp/Ising_Spin /home/shiroha/Workspace/cp/Ising_Spin /home/shiroha/Workspace/cp/Ising_Spin/build /home/shiroha/Workspace/cp/Ising_Spin/build /home/shiroha/Workspace/cp/Ising_Spin/build/CMakeFiles/UnitTests2_catch2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/UnitTests2_catch2.dir/depend

