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

typedef LDLinv(*func1)(LLMatOrd &a);
typedef LDLinv(*func2)(LLMatOrd_vector2 &a);

void elapsed_time(func1 f1, LLMatOrd llmat, std::string func_name, int rep) {
    long count = 0;
    std::vector<LLMatOrd> llmats(rep, llmat);
    for (int i = 0; i < rep; ++i){
        auto start = std::chrono::steady_clock::now();
        f1(llmats[i]);
        auto end = std::chrono::steady_clock::now();   
        count += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    }
    std::cout << " " <<  count/(float)(rep)/1e9 << " ";
}

void elapsed_time(func2 f2, LLMatOrd_vector2 llmat, std::string func_name, int rep) {
    long count = 0;
    std::vector<LLMatOrd_vector2> llmats(rep, llmat);
    for (int i = 0; i < rep; ++i){
        auto start = std::chrono::steady_clock::now();
        f2(llmats[i]);
        auto end = std::chrono::steady_clock::now();   
        count += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    }
    std::cout << " " <<  count/(float)(rep)/1e9 << " ";
}

