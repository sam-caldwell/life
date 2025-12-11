/**
 * @file Logger.cpp
 */
#include "Logger.h"

#include <mutex>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <thread>

namespace {
std::mutex g_logMtx;
std::ofstream g_log;
bool g_inited = false;

static std::string basenameFromPath(const std::string& p) {
    if (p.empty()) return std::string("app");
    size_t pos = p.find_last_of("/\\");
    if (pos == std::string::npos) return p;
    if (pos + 1 >= p.size()) return std::string("app");
    return p.substr(pos + 1);
}

static std::string nowTs() {
    using namespace std::chrono;
    auto tp = system_clock::now();
    auto t = system_clock::to_time_t(tp);
    auto ms = duration_cast<milliseconds>(tp.time_since_epoch()) % 1000;
    std::tm tmv;
#if defined(_WIN32)
    localtime_s(&tmv, &t);
#else
    localtime_r(&t, &tmv);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tmv, "%Y-%m-%d %H:%M:%S")
        << '.' << std::setw(3) << std::setfill('0') << ms.count();
    return oss.str();
}
}

void Logger::initFromArgv0(const char* argv0) {
    std::string base = basenameFromPath(argv0 ? std::string(argv0) : std::string("app"));
    if (base.empty()) base = "app";
    std::string file = std::string("./") + base + ".log";
    init(file);
}

void Logger::init(const std::string& filename) {
    std::lock_guard<std::mutex> lock(g_logMtx);
    if (g_inited) return;
    // Append to preserve prior runs; we also mark a session header.
    g_log.open(filename, std::ios::out | std::ios::app);
    if (g_log.is_open()) {
        g_inited = true;
        g_log << "===== session start " << nowTs() << " =====" << std::endl;
        g_log.flush();
    }
}

void Logger::shutdown() {
    std::lock_guard<std::mutex> lock(g_logMtx);
    if (!g_inited) return;
    g_log << "===== session end   " << nowTs() << " =====" << std::endl;
    g_log.flush();
    g_log.close();
    g_inited = false;
}

void Logger::logImpl(const char* level, const std::string& msg) {
    std::lock_guard<std::mutex> lock(g_logMtx);
    if (!g_inited || !g_log.is_open()) return;
    std::ostringstream tid;
    tid << std::this_thread::get_id();
    g_log << nowTs() << " [" << level << "] [t:" << tid.str() << "] " << msg << std::endl;
    g_log.flush();
}

void Logger::info(const std::string& msg) { logImpl("INFO", msg); }
void Logger::warn(const std::string& msg) { logImpl("WARN", msg); }
void Logger::error(const std::string& msg) { logImpl("ERROR", msg); }
void Logger::debug(const std::string& msg) { logImpl("DEBUG", msg); }

void Logger::logException(const std::string& where, const std::exception& e) {
    logImpl("ERROR", where + ": " + e.what());
}

void Logger::logUnknownException(const std::string& where) {
    logImpl("ERROR", where + ": unknown exception");
}
