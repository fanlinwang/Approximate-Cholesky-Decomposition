add_library(approxChol STATIC approxChol.cpp approxCholTypes.cpp conditionNumber.cpp)
target_include_directories(approxChol PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# ${Boost_INCLUDE_DIRS}

add_executable(main main.cpp)
target_link_libraries(main PRIVATE approxChol)

