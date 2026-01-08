//
// Created by xumiz on 2026/1/1.
//

#ifndef TAIYI_PROFILER_H
#define TAIYI_PROFILER_H

#include <chrono>
#include <cstdint>

class Profiler {
public:
    class ScopedTimer {
    public:
        ScopedTimer(const char* name, const char* file, int line);
        ~ScopedTimer();

        ScopedTimer(const ScopedTimer&) = delete;
        ScopedTimer& operator=(const ScopedTimer&) = delete;

    private:
        const char* name_;
        const char* file_;
        int line_;
        std::chrono::steady_clock::time_point start_;
    };

    static void BeginFrame();
    static void EndFrame();
    static void Shutdown();

private:
    static void Record(const char* name, const char* file, int line, double elapsed_ms);
};

#define PROFILE_SCOPE(name) ::Profiler::ScopedTimer profilerScope##__LINE__(name, __FILE__, __LINE__)
#define PROFILE_FUNCTION() PROFILE_SCOPE(__func__)

#endif // TAIYI_PROFILER_H
