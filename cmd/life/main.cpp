/**
 * @file main.cpp
 * @brief Classic Life entry: initializes ncurses, installs signal handlers, runs UI loop and graceful shutdown.
 *
 * @copyright Copyright (c) 2025 Sam Caldwell. Released under the MIT License.
 */
#include <ncurses.h>
#include <chrono>
#include <thread>
#include <cstdlib>
#include <csignal>
#include <cstdlib>
#include <exception>
#include "Board.h"
#include "Logger.h"

static volatile sig_atomic_t g_stop = 0;
static void handle_signal(int) { g_stop = 1; }

// Forward decl for use in signal handler
static void init_colors();

static bool g_curses_inited = false;
static volatile sig_atomic_t g_needs_full_redraw = 0;
static void atexit_cleanup() {
    if (g_curses_inited) {
        endwin();
        g_curses_inited = false;
    }
}

// Handle terminal suspension (Ctrl+Z): restore tty before stopping.
static void handle_sigtstp(int) {
    if (g_curses_inited) {
        def_prog_mode(); // save current curses state
        endwin();        // restore tty modes for the shell
        g_curses_inited = false;
    }
    // Revert to default action and re-raise to actually stop the process
    struct sigaction sa{};
    sa.sa_handler = SIG_DFL;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGTSTP, &sa, nullptr);
    raise(SIGTSTP);
}

// Resume after suspension: restore curses program mode and redraw UI
static void handle_sigcont(int) {
    // Reinstall our SIGTSTP handler
    struct sigaction st{}; st.sa_handler = handle_sigtstp; sigemptyset(&st.sa_mask); st.sa_flags = 0; sigaction(SIGTSTP, &st, nullptr);

    reset_prog_mode(); // restore saved curses state
    refresh();         // refresh the screen
    cbreak();
    noecho();
    curs_set(0);
    keypad(stdscr, TRUE);
    nodelay(stdscr, TRUE);
    timeout(0);
    init_colors();
    clearok(stdscr, TRUE);
    refresh();
    g_curses_inited = true;
    g_needs_full_redraw = 1;
}

/** @brief Initialize ncurses color pairs used to render weight spectrum. */
static void init_colors() {
    if (!has_colors()) return;
    start_color();
    use_default_colors();
    // Weight-based color spectrum (darker to lighter)
    // 1: 0      -> black
    // 2: 10     -> blue
    // 3: 20     -> green
    // 4: 30     -> cyan
    // 5: 50     -> yellow
    // 6: 70     -> magenta
    // 7: 80     -> red
    // 8: 90/100 -> white (100 drawn in bold)
    init_pair(1, COLOR_BLACK, -1);
    init_pair(2, COLOR_BLUE, -1);
    init_pair(3, COLOR_GREEN, -1);
    init_pair(4, COLOR_CYAN, -1);
    init_pair(5, COLOR_YELLOW, -1);
    init_pair(6, COLOR_MAGENTA, -1);
    init_pair(7, COLOR_RED, -1);
    init_pair(8, COLOR_WHITE, -1);
}

/** @brief Program entry: sets up terminal UI, seeds the board, handles input, and exits cleanly on signals. */
int main(int argc, char** argv) {
    Logger::initFromArgv0((argc > 0) ? argv[0] : "life");
    Logger::info("life starting");
    std::set_terminate([]{
        try {
            auto ep = std::current_exception();
            if (ep) {
                try { std::rethrow_exception(ep); }
                catch (const std::exception& e) { Logger::logException("std::terminate (life)", e); }
                catch (...) { Logger::logUnknownException("std::terminate (life)"); }
            } else {
                Logger::error("std::terminate (life): no active exception");
            }
        } catch (...) {}
        if (g_curses_inited) { endwin(); }
        Logger::shutdown();
        std::_Exit(1);
    });
    try {
    struct sigaction sa{};
    sa.sa_handler = handle_signal;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGINT, &sa, nullptr);
    sigaction(SIGTERM, &sa, nullptr);
    // Suspend/resume handlers
    struct sigaction st{}; st.sa_handler = handle_sigtstp; sigemptyset(&st.sa_mask); st.sa_flags = 0; sigaction(SIGTSTP, &st, nullptr);
    struct sigaction sc{}; sc.sa_handler = handle_sigcont; sigemptyset(&sc.sa_mask); sc.sa_flags = 0; sigaction(SIGCONT, &sc, nullptr);

    initscr();
    g_curses_inited = true;
    std::atexit(atexit_cleanup);
    cbreak();
    noecho();
    curs_set(0);
    keypad(stdscr, TRUE);
    nodelay(stdscr, TRUE); // non-blocking getch
    timeout(0);
    init_colors();

    int rows, cols;
    getmaxyx(stdscr, rows, cols);
    int gridH = rows - 1;
    int gridW = cols;
    if (gridH < 1 || gridW < 1) {
        Logger::error("terminal too small");
        endwin();
        return 1;
    }

    Board board(gridW, gridH);

    // Seed with 10% of grid filled by default
    unsigned initial = (unsigned)((gridW * gridH) / 10);
    board.setWindow(stdscr);
    board.reseed(initial);
    board.draw(stdscr); // initial full draw only
    board.setRunning(false); // start paused
    Logger::info("board initialized: " + std::to_string(gridW) + "x" + std::to_string(gridH));

    bool done = false;
    while (!done) {
        if (g_stop) done = true;
        if (g_needs_full_redraw) {
            board.draw(stdscr);
            board.drawStatusLine(stdscr);
            g_needs_full_redraw = 0;
        }
        // Update status line only (incremental board updates happen per automaton)
        board.drawStatusLine(stdscr);

        int ch = getch();
        switch (ch) {
            case 'q':
            case 'Q':
                board.setStatusNote("terminating");
                board.drawStatusLine(stdscr);
                Logger::info("quit requested");
                done = true; break;
            case 's': case 'S':
                board.toggleRunning();
                Logger::info(std::string("running = ") + (board.isRunning()?"true":"false"));
                break;
            case 'p': case 'P':
                board.setRunning(false);
                Logger::info("paused");
                break;
            case 'r': case 'R': {
                unsigned cnt = (unsigned)((gridW * gridH) / 10);
                board.reseed(cnt);
                Logger::info("reseed: count=" + std::to_string(cnt));
                break; }
            case 'c': case 'C':
                board.clear();
                Logger::info("clear requested");
                break;
            case '+': {
                int ms = board.getStepDelayMs() - 10;
                board.setStepDelayMs(ms);
                Logger::info("delay set(ms): " + std::to_string(board.getStepDelayMs()));
                break; }
            case '-': {
                int ms = board.getStepDelayMs() + 10;
                board.setStepDelayMs(ms);
                Logger::info("delay set(ms): " + std::to_string(board.getStepDelayMs()));
                break; }
            default:
                break;
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(16)); // ~60 FPS UI
    }

    // Graceful shutdown: stop sim threads and cleanup
    board.setRunning(false);
    board.clear();
    endwin();
    g_curses_inited = false;
    Logger::info("life terminating");
    Logger::shutdown();
    return 0;
    } catch (const std::exception& e) {
        if (g_curses_inited) { endwin(); g_curses_inited = false; }
        Logger::logException("unhandled exception (life)", e);
        Logger::shutdown();
        return 2;
    } catch (...) {
        if (g_curses_inited) { endwin(); g_curses_inited = false; }
        Logger::logUnknownException("unhandled exception (life)");
        Logger::shutdown();
        return 2;
    }
}
