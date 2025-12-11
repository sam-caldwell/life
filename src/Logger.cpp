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
#include <cstdlib>

namespace {
std::mutex g_logMtx;
std::ofstream g_log;
bool g_inited = false;
std::atomic<int> g_level{static_cast<int>(Logger::Level::Info)};

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

static void setLevelFromEnv() {
    const char* s = std::getenv("LOG_LEVEL");
    if (!s) return;
    std::string v(s);
    for (auto& c : v) c = (char)std::tolower(c);
    if (v == "debug") g_level.store((int)Logger::Level::Debug);
    else if (v == "info") g_level.store((int)Logger::Level::Info);
    else if (v == "warn" || v == "warning") g_level.store((int)Logger::Level::Warn);
    else if (v == "error") g_level.store((int)Logger::Level::Error);
    else if (v == "none" || v == "off") g_level.store((int)Logger::Level::None);
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
        setLevelFromEnv();
        g_log << "===== session start " << nowTs() << " =====" << '\n';
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

void Logger::logImpl(Level lvl, const std::string& msg) {
    std::lock_guard<std::mutex> lock(g_logMtx);
    if (!g_inited || !g_log.is_open()) return;
    if ((int)lvl < g_level.load()) return;
    std::ostringstream tid;
    tid << std::this_thread::get_id();
    const char* name = (lvl == Level::Debug ? "DEBUG" : lvl == Level::Info ? "INFO" : lvl == Level::Warn ? "WARN" : "ERROR");
    g_log << nowTs() << " [" << name << "] [t:" << tid.str() << "] " << msg << '\n';
    if (lvl >= Level::Warn) g_log.flush();
}

void Logger::info(const std::string& msg) { logImpl(Level::Info, msg); }
void Logger::warn(const std::string& msg) { logImpl(Level::Warn, msg); }
void Logger::error(const std::string& msg) { logImpl(Level::Error, msg); }
void Logger::debug(const std::string& msg) { logImpl(Level::Debug, msg); }

void Logger::logException(const std::string& where, const std::exception& e) {
    logImpl(Level::Error, where + ": " + e.what());
}

void Logger::logUnknownException(const std::string& where) {
    logImpl(Level::Error, where + ": unknown exception");
}

void Logger::setLevel(Level lvl) { g_level.store((int)lvl); }
Logger::Level Logger::level() { return (Level)g_level.load(); }
