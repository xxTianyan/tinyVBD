//
// Created by xumiz on 2026/1/1.
//

#include "Profiler.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
    struct Stat {
        double total_ms = 0.0;
        double max_ms = 0.0;
        std::uint64_t count = 0;
    };

    struct ProfilerState {
        std::unordered_map<std::string, Stat> stats;
        std::vector<double> frame_times_ms;
        std::ofstream frame_csv;
        std::uint64_t frame_index = 0;
        std::chrono::steady_clock::time_point frame_start;
        bool frame_active = false;
    };

    ProfilerState& State() {
        static ProfilerState state;
        return state;
    }

    std::string MakeKey(const char* name, const char* file, int line) {
        return std::string(name) + " @ " + file + ":" + std::to_string(line);
    }

    void EnsureCsvOpen() {
        auto& state = State();
        if (!state.frame_csv.is_open()) {
            state.frame_csv.open("profiling_frames.csv", std::ios::out | std::ios::trunc);
            if (state.frame_csv.is_open()) {
                state.frame_csv << "frame,frame_ms\n";
            }
        }
    }
} // namespace

Profiler::ScopedTimer::ScopedTimer(const char* name, const char* file, int line)
    : name_(name),
      file_(file),
      line_(line),
      start_(std::chrono::steady_clock::now()) {}

Profiler::ScopedTimer::~ScopedTimer() {
    const auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double, std::milli> elapsed = end - start_;
    Profiler::Record(name_, file_, line_, elapsed.count());
}

void Profiler::BeginFrame() {
    auto& state = State();
    state.frame_start = std::chrono::steady_clock::now();
    state.frame_active = true;
    state.frame_index += 1;
}

void Profiler::EndFrame() {
    auto& state = State();
    if (!state.frame_active) {
        return;
    }
    const auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double, std::milli> elapsed = end - state.frame_start;
    const double frame_ms = elapsed.count();
    state.frame_times_ms.push_back(frame_ms);
    EnsureCsvOpen();
    if (state.frame_csv.is_open()) {
        state.frame_csv << state.frame_index << "," << std::fixed << std::setprecision(3) << frame_ms << "\n";
    }
    state.frame_active = false;
}

void Profiler::Shutdown() {
    auto& state = State();
    if (state.frame_active) {
        EndFrame();
    }

    std::vector<std::pair<std::string, Stat>> entries;
    entries.reserve(state.stats.size());
    for (const auto& entry : state.stats) {
        entries.emplace_back(entry.first, entry.second);
    }
    std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
        return a.second.total_ms > b.second.total_ms;
    });

    std::ofstream report("profiling_report.txt", std::ios::out | std::ios::trunc);
    if (report.is_open()) {
        report << "==== Profiling Report ====\n";
        report << "Scopes sorted by total time (ms):\n";
        for (const auto& entry : entries) {
            const auto& stat = entry.second;
            const double avg_ms = (stat.count > 0) ? (stat.total_ms / static_cast<double>(stat.count)) : 0.0;
            report << std::setw(10) << std::fixed << std::setprecision(3) << stat.total_ms << " ms | "
                   << "avg " << std::setw(7) << avg_ms << " ms | "
                   << "max " << std::setw(7) << stat.max_ms << " ms | "
                   << "count " << stat.count << " | "
                   << entry.first << "\n";
        }
        report << "\nFrame times stored in profiling_frames.csv\n";
    }

    if (!state.frame_times_ms.empty()) {
        double sum = 0.0;
        double min_ms = state.frame_times_ms.front();
        double max_ms = state.frame_times_ms.front();
        for (double t : state.frame_times_ms) {
            sum += t;
            min_ms = std::min(min_ms, t);
            max_ms = std::max(max_ms, t);
        }
        const double avg_ms = sum / static_cast<double>(state.frame_times_ms.size());
        std::cout << "Profiling report written to profiling_report.txt\n";
        std::cout << "Frame times CSV written to profiling_frames.csv\n";
        std::cout << "Frame stats: avg " << std::fixed << std::setprecision(3) << avg_ms
                  << " ms, min " << min_ms << " ms, max " << max_ms << " ms\n";
    }

    if (state.frame_csv.is_open()) {
        state.frame_csv.flush();
        state.frame_csv.close();
    }
}

void Profiler::Record(const char* name, const char* file, int line, double elapsed_ms) {
    auto& state = State();
    const std::string key = MakeKey(name, file, line);
    auto& stat = state.stats[key];
    stat.total_ms += elapsed_ms;
    stat.max_ms = std::max(stat.max_ms, elapsed_ms);
    stat.count += 1;
}
