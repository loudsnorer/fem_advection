# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/fariz/CLionProjects/fem_advection

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/fariz/CLionProjects/fem_advection/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/fem_advection.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/fem_advection.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/fem_advection.dir/flags.make

CMakeFiles/fem_advection.dir/main.cpp.o: CMakeFiles/fem_advection.dir/flags.make
CMakeFiles/fem_advection.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/fariz/CLionProjects/fem_advection/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/fem_advection.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fem_advection.dir/main.cpp.o -c /Users/fariz/CLionProjects/fem_advection/main.cpp

CMakeFiles/fem_advection.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fem_advection.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/fariz/CLionProjects/fem_advection/main.cpp > CMakeFiles/fem_advection.dir/main.cpp.i

CMakeFiles/fem_advection.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fem_advection.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/fariz/CLionProjects/fem_advection/main.cpp -o CMakeFiles/fem_advection.dir/main.cpp.s

CMakeFiles/fem_advection.dir/solve.cpp.o: CMakeFiles/fem_advection.dir/flags.make
CMakeFiles/fem_advection.dir/solve.cpp.o: ../solve.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/fariz/CLionProjects/fem_advection/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/fem_advection.dir/solve.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fem_advection.dir/solve.cpp.o -c /Users/fariz/CLionProjects/fem_advection/solve.cpp

CMakeFiles/fem_advection.dir/solve.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fem_advection.dir/solve.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/fariz/CLionProjects/fem_advection/solve.cpp > CMakeFiles/fem_advection.dir/solve.cpp.i

CMakeFiles/fem_advection.dir/solve.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fem_advection.dir/solve.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/fariz/CLionProjects/fem_advection/solve.cpp -o CMakeFiles/fem_advection.dir/solve.cpp.s

# Object files for target fem_advection
fem_advection_OBJECTS = \
"CMakeFiles/fem_advection.dir/main.cpp.o" \
"CMakeFiles/fem_advection.dir/solve.cpp.o"

# External object files for target fem_advection
fem_advection_EXTERNAL_OBJECTS =

fem_advection: CMakeFiles/fem_advection.dir/main.cpp.o
fem_advection: CMakeFiles/fem_advection.dir/solve.cpp.o
fem_advection: CMakeFiles/fem_advection.dir/build.make
fem_advection: CMakeFiles/fem_advection.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/fariz/CLionProjects/fem_advection/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable fem_advection"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fem_advection.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/fem_advection.dir/build: fem_advection

.PHONY : CMakeFiles/fem_advection.dir/build

CMakeFiles/fem_advection.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fem_advection.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fem_advection.dir/clean

CMakeFiles/fem_advection.dir/depend:
	cd /Users/fariz/CLionProjects/fem_advection/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/fariz/CLionProjects/fem_advection /Users/fariz/CLionProjects/fem_advection /Users/fariz/CLionProjects/fem_advection/cmake-build-debug /Users/fariz/CLionProjects/fem_advection/cmake-build-debug /Users/fariz/CLionProjects/fem_advection/cmake-build-debug/CMakeFiles/fem_advection.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fem_advection.dir/depend
