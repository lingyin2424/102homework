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
CMAKE_SOURCE_DIR = /home/lingyin/myCode/hw_102/all/code2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lingyin/myCode/hw_102/all/code2/build

# Include any dependencies generated for this target.
include CMakeFiles/code2.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/code2.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/code2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/code2.dir/flags.make

CMakeFiles/code2.dir/code2.o: CMakeFiles/code2.dir/flags.make
CMakeFiles/code2.dir/code2.o: ../code2.cpp
CMakeFiles/code2.dir/code2.o: CMakeFiles/code2.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lingyin/myCode/hw_102/all/code2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/code2.dir/code2.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/code2.dir/code2.o -MF CMakeFiles/code2.dir/code2.o.d -o CMakeFiles/code2.dir/code2.o -c /home/lingyin/myCode/hw_102/all/code2/code2.cpp

CMakeFiles/code2.dir/code2.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/code2.dir/code2.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lingyin/myCode/hw_102/all/code2/code2.cpp > CMakeFiles/code2.dir/code2.i

CMakeFiles/code2.dir/code2.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/code2.dir/code2.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lingyin/myCode/hw_102/all/code2/code2.cpp -o CMakeFiles/code2.dir/code2.s

# Object files for target code2
code2_OBJECTS = \
"CMakeFiles/code2.dir/code2.o"

# External object files for target code2
code2_EXTERNAL_OBJECTS =

code2: CMakeFiles/code2.dir/code2.o
code2: CMakeFiles/code2.dir/build.make
code2: /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.5.15.3
code2: /usr/lib/x86_64-linux-gnu/libQt5Gui.so.5.15.3
code2: /usr/lib/x86_64-linux-gnu/libQt5Core.so.5.15.3
code2: CMakeFiles/code2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lingyin/myCode/hw_102/all/code2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable code2"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/code2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/code2.dir/build: code2
.PHONY : CMakeFiles/code2.dir/build

CMakeFiles/code2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/code2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/code2.dir/clean

CMakeFiles/code2.dir/depend:
	cd /home/lingyin/myCode/hw_102/all/code2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lingyin/myCode/hw_102/all/code2 /home/lingyin/myCode/hw_102/all/code2 /home/lingyin/myCode/hw_102/all/code2/build /home/lingyin/myCode/hw_102/all/code2/build /home/lingyin/myCode/hw_102/all/code2/build/CMakeFiles/code2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/code2.dir/depend

