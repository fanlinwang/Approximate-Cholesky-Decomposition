#[[
add_executable(test_LDLsolver test_LDLsolver.cpp ${PROJECT_SOURCES})
target_link_libraries(test_LDLsolver PRIVATE approxChol)
add_executable(test_approxChol test_approxChol.cpp ${PROJECT_SOURCES})  
target_link_libraries(test_approxChol PRIVATE approxChol)
]]
add_executable(test_LLMatOrd test_LLMatOrd.cpp ${PROJECT_SOURCES})
target_link_libraries(test_LLMatOrd PRIVATE approxChol)