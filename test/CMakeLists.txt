
add_executable(test_LDLsolver test_LDLsolver.cpp)
target_link_libraries(test_LDLsolver PRIVATE approxChol)

add_executable(test_approxChol test_approxChol.cpp)  
target_link_libraries(test_approxChol PRIVATE approxChol)

add_executable(test_LLMatOrd test_LLMatOrd.cpp)
target_link_libraries(test_LLMatOrd PRIVATE approxChol)

add_executable(test_pcg test_pcg.cpp)
target_link_libraries(test_pcg PRIVATE approxChol)

add_executable(test_perf test_solver.cpp)
target_link_libraries(test_perf PRIVATE approxChol)
