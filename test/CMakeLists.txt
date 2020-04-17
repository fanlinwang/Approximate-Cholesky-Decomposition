add_executable(test_LDLsolver test_LDLsolver.cpp ${PROJECT_SOURCES})
target_link_libraries(test_LDLsolver PRIVATE approxChol)
