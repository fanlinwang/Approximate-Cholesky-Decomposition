#pragma once

#include <stdexcept>
#include <string>
#include <chrono>

struct test_exception 
{
   test_exception(const std::string& s) : error_message(s) {}
   std::string error_message; 
};

inline void test(const bool pred, const std::string& error_message = "")
{
    std::string msg;
    if(!error_message.empty()) {
        msg = "test failed: " + error_message;
    } else {
        msg = "test failed.";
    }

    if (!pred)
        throw std::runtime_error(msg);
}

typedef LDLinv(*func1)(LLMatOrd a);
typedef LDLinv(*func2)(LLMatOrd_vector2 a);

void elapsed_time(func1 f1, LLMatOrd llmat, std::string func_name) {
    auto start = std::chrono::steady_clock::now();
    f1(llmat);
    auto end = std::chrono::steady_clock::now();        
	std::cout << func_name << " elapsed time in nanoseconds : " 
		<< std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;

}
void elapsed_time(func2 f2, LLMatOrd_vector2 llmat, std::string func_name) {
    auto start = std::chrono::steady_clock::now();
    f2(llmat);
    auto end = std::chrono::steady_clock::now();        
    std::cout << func_name << " elapsed time in nanoseconds : " 
		<< std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
}

