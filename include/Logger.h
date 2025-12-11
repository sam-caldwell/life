/**
 * @file Logger.h
 * @brief Minimal thread-safe file logger writing to ./<command>.log
 */
#pragma once

#include <string>

class Logger {
public:
    // Initialize using argv[0] to derive <command>.log path.
    static void initFromArgv0(const char* argv0);
    // Initialize explicitly with a filename (relative or absolute).
    static void init(const std::string& filename);
    // Flush and close the log file; safe to call multiple times.
    static void shutdown();

    static void info(const std::string& msg);
    static void warn(const std::string& msg);
    static void error(const std::string& msg);
    static void debug(const std::string& msg);

    // Convenience helpers for exception logging
    static void logException(const std::string& where, const std::exception& e);
    static void logUnknownException(const std::string& where);

private:
    static void logImpl(const char* level, const std::string& msg);
};
