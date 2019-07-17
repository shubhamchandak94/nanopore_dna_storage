# CMake generated Testfile for 
# Source directory: /raid/nanopore/shubham/nanopore_dna_storage/flappie
# Build directory: /raid/nanopore/shubham/nanopore_dna_storage/flappie/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(unittest "/raid/nanopore/shubham/nanopore_dna_storage/flappie/build/flappie_unittest")
set_tests_properties(unittest PROPERTIES  WORKING_DIRECTORY "/raid/nanopore/shubham/nanopore_dna_storage/flappie/src/test/")
add_test(test_call "flappie" "/raid/nanopore/shubham/nanopore_dna_storage/flappie/reads")
add_test(test_licence "flappie" "--licence")
add_test(test_license "flappie" "--license")
add_test(test_help "flappie" "--help")
add_test(test_version "flappie" "--version")
