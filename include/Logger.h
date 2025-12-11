/**
 * @file Logger.h
 * @brief Minimal thread-safe file logger writing to ./<command>.log
 */
#pragma once

#include <string>

class Logger {
public:
    enum class Level { Debug = 0, Info = 1, Warn = 2, Error = 3, None = 4 };

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

    // Control log level (default: Info). Also honors LOG_LEVEL env (debug, info, warn, error, none)
    static void setLevel(Level lvl);
    static Level level();

private:
    static void logImpl(Level lvl, const std::string& msg);
};
