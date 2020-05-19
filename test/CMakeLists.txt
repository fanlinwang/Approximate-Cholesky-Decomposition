<<<<<<< HEAD

add_executable(test_LDLsolver test_LDLsolver.cpp)
target_link_libraries(test_LDLsolver PRIVATE approxChol)

add_executable(test_approxChol test_approxChol.cpp)  
target_link_libraries(test_approxChol PRIVATE approxChol)

add_executable(test_LLMatOrd test_LLMatOrd.cpp)
target_link_libraries(test_LLMatOrd PRIVATE approxChol)

add_executable(test_pcg test_pcg.cpp)
target_link_libraries(test_pcg PRIVATE approxChol)

add_executable(test_flops test_flops.cpp)
=======

add_executable(test_LDLsolver test_LDLsolver.cpp)
target_link_libraries(test_LDLsolver PRIVATE approxChol)

add_executable(test_approxChol test_approxChol.cpp)  
target_link_libraries(test_approxChol PRIVATE approxChol)

add_executable(test_LLMatOrd test_LLMatOrd.cpp)
target_link_libraries(test_LLMatOrd PRIVATE approxChol)

add_executable(test_pcg test_pcg.cpp)
target_link_libraries(test_pcg PRIVATE approxChol)

add_executable(test_validation test_validation.cpp)
target_link_libraries(test_validation PRIVATE approxChol)
>>>>>>> 9e4f2a498c3ed036f6b87b3b34b16da81fb84057
