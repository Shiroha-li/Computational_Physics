# CMake generated Testfile for 
# Source directory: /home/shiroha/Workspace/cp/Ising_Spin
# Build directory: /home/shiroha/Workspace/cp/Ising_Spin/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTests1 "UnitTests1_catch2")
set_tests_properties(UnitTests1 PROPERTIES  _BACKTRACE_TRIPLES "/home/shiroha/Workspace/cp/Ising_Spin/CMakeLists.txt;21;add_test;/home/shiroha/Workspace/cp/Ising_Spin/CMakeLists.txt;0;")
add_test(UnitTests2 "UnitTests2_catch2")
set_tests_properties(UnitTests2 PROPERTIES  _BACKTRACE_TRIPLES "/home/shiroha/Workspace/cp/Ising_Spin/CMakeLists.txt;26;add_test;/home/shiroha/Workspace/cp/Ising_Spin/CMakeLists.txt;0;")
subdirs("_deps/catch2-build")
