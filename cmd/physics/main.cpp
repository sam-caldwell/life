/**
 * @file main_physics.cpp
 * @brief Physics entry point for the Newtonian variant.
 */
#include <ncurses.h>
#include <chrono>
#include <thread>
#include <csignal>
#include <cstdlib>
#include <exception>
#include <cstring>
#include <string>
#include <cerrno>
#include <random>
#include "PhysicsWorld.h"
#include "Logger.h"

static volatile sig_atomic_t g_stop = 0;
static void handle_signal(int) { g_stop = 1; }

// Forward decl for signal handler
static void init_colors_physics();

static bool g_curses_inited = false;
static volatile sig_atomic_t g_needs_full_redraw = 0;
static void atexit_cleanup() {
    if (g_curses_inited) {
        endwin();
        g_curses_inited = false;
    }
}

// Suspend: restore tty, then stop process with default action
static void handle_sigtstp(int) {
    if (g_curses_inited) {
        def_prog_mode();
        endwin();
        g_curses_inited = false;
    }
    struct sigaction sa{}; sa.sa_handler = SIG_DFL; sigemptyset(&sa.sa_mask); sa.sa_flags = 0; sigaction(SIGTSTP, &sa, nullptr);
    raise(SIGTSTP);
}

// Resume: restore curses state and redraw
static void handle_sigcont(int) {
    struct sigaction st{}; st.sa_handler = handle_sigtstp; sigemptyset(&st.sa_mask); st.sa_flags = 0; sigaction(SIGTSTP, &st, nullptr);
    reset_prog_mode();
    refresh();
    cbreak();
    noecho();
    curs_set(0);
    keypad(stdscr, TRUE);
    nodelay(stdscr, TRUE);
    timeout(0);
    init_colors_physics();
    clearok(stdscr, TRUE);
    refresh();
    g_curses_inited = true;
    g_needs_full_redraw = 1;
}

static void init_colors_physics() {
    if (!has_colors()) return;
    start_color();
    use_default_colors();
    // Define 16 color pairs. We map 1..8 to standard colors and 9..16 to the same colors rendered with A_BOLD.
    // 1: BLACK, 2: RED, 3: GREEN, 4: YELLOW, 5: BLUE, 6: MAGENTA, 7: CYAN, 8: WHITE
    init_pair(1, COLOR_BLACK, -1);
    init_pair(2, COLOR_RED, -1);
    init_pair(3, COLOR_GREEN, -1);
    init_pair(4, COLOR_YELLOW, -1);
    init_pair(5, COLOR_BLUE, -1);
    init_pair(6, COLOR_MAGENTA, -1);
    init_pair(7, COLOR_CYAN, -1);
    init_pair(8, COLOR_WHITE, -1);
    // 9..16 reuse the same colors; brightness is achieved via A_BOLD when drawing
    init_pair(9,  COLOR_BLACK, -1);
    init_pair(10, COLOR_RED, -1);
    init_pair(11, COLOR_GREEN, -1);
    init_pair(12, COLOR_YELLOW, -1);
    init_pair(13, COLOR_BLUE, -1);
    init_pair(14, COLOR_MAGENTA, -1);
    init_pair(15, COLOR_CYAN, -1);
    init_pair(16, COLOR_WHITE, -1);
}

static bool parseFloat(const char* s, float& out) {
    if (!s) return false;
    char* end = nullptr;
    errno = 0;
    double v = std::strtod(s, &end);
    if (end == s || errno != 0) return false;
    out = static_cast<float>(v);
    return true;
}

int main(int argc, char** argv) {
    Logger::initFromArgv0((argc > 0) ? argv[0] : "physics");
    Logger::info("physics starting");
    std::set_terminate([]{
        try {
            auto ep = std::current_exception();
            if (ep) {
                try { std::rethrow_exception(ep); }
                catch (const std::exception& e) { Logger::logException("std::terminate (physics)", e); }
                catch (...) { Logger::logUnknownException("std::terminate (physics)"); }
            } else {
                Logger::error("std::terminate (physics): no active exception");
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
    struct sigaction st{}; st.sa_handler = handle_sigtstp; sigemptyset(&st.sa_mask); st.sa_flags = 0; sigaction(SIGTSTP, &st, nullptr);
    struct sigaction sc{}; sc.sa_handler = handle_sigcont; sigemptyset(&sc.sa_mask); sc.sa_flags = 0; sigaction(SIGCONT, &sc, nullptr);

    initscr();
    g_curses_inited = true;
    std::atexit(atexit_cleanup);
    cbreak();
    noecho();
    curs_set(0);
    keypad(stdscr, TRUE);
    nodelay(stdscr, TRUE);
    timeout(0);
    init_colors_physics();

    int rows, cols; getmaxyx(stdscr, rows, cols);
    int gridH = rows - 1;
    int gridW = cols;
    if (gridH < 1 || gridW < 1) { Logger::error("terminal too small"); endwin(); Logger::shutdown(); return 1; }

    // Defaults
    float gravity = 9.8f;
    float radius = 1.0f;
    float restitution = 0.98f;

    // Env overrides
    float tmp;
    const char* eg = std::getenv("PHYSICS_G");
    if (parseFloat(eg, tmp)) gravity = tmp;
    const char* er = std::getenv("PHYSICS_RADIUS");
    if (parseFloat(er, tmp)) radius = tmp;
    const char* ee = std::getenv("PHYSICS_RESTITUTION");
    if (parseFloat(ee, tmp)) restitution = tmp;

    // Arg overrides
    for (int i = 1; i < argc; ++i) {
        std::string a(argv[i]);
        auto read_next = [&](int& idx) -> const char* {
            if (idx + 1 < argc) return argv[++idx];
            return nullptr;
        };
        if (a == "-g" || a == "--gravity") {
            const char* v = read_next(i); if (parseFloat(v, tmp)) gravity = tmp;
        } else if (a.rfind("--gravity=", 0) == 0) {
            const char* v = a.c_str() + std::strlen("--gravity="); if (parseFloat(v, tmp)) gravity = tmp;
        } else if (a == "-r" || a == "--radius") {
            const char* v = read_next(i); if (parseFloat(v, tmp)) radius = tmp;
        } else if (a.rfind("--radius=", 0) == 0) {
            const char* v = a.c_str() + std::strlen("--radius="); if (parseFloat(v, tmp)) radius = tmp;
        } else if (a == "-e" || a == "--restitution") {
            const char* v = read_next(i); if (parseFloat(v, tmp)) restitution = tmp;
        } else if (a.rfind("--restitution=", 0) == 0) {
            const char* v = a.c_str() + std::strlen("--restitution="); if (parseFloat(v, tmp)) restitution = tmp;
        }
    }

    // Validation/clamping
    if (!(gravity >= 0.0f)) gravity = 9.8f;
    if (!(radius >= 0.0f)) radius = 1.0f;
    if (!(restitution >= 0.0f && restitution <= 1.0f)) restitution = 0.98f;

    PhysicsWorld world(gridW, gridH, gravity, radius, restitution);
    world.setWindow(stdscr);
    // Initial population: random count up to 20 (and within cap)
    {
        std::random_device rd{};
        std::mt19937 gen(rd());
        int cap = (int)std::min<size_t>(20, world.maxParticles());
        if (cap < 1) cap = 1;
        std::uniform_int_distribution<int> d(1, cap);
        unsigned n = (unsigned)d(gen);
        world.reseed(n);
        Logger::info("physics reseed: count=" + std::to_string(n));
    }
    world.drawAll(stdscr);
    world.drawStatusLine(stdscr);
    world.setRunning(false);

    bool done = false;
    using namespace std::chrono;
    auto lastStep = steady_clock::now();
    while (!done) {
        if (g_stop) done = true;
        if (g_needs_full_redraw) {
            world.drawAll(stdscr);
            world.drawStatusLine(stdscr);
            g_needs_full_redraw = 0;
        }
        // Run simulation step based on elapsed time and configured cadence
        if (world.isRunning()) {
            auto now = steady_clock::now();
            auto msDesired = world.getStepDelayMs();
            if (msDesired < 1) msDesired = 1;
            auto elapsed = duration_cast<milliseconds>(now - lastStep).count();
            if (elapsed >= msDesired) {
                float dt = static_cast<float>(elapsed) / 1000.0f;
                world.step(dt);
                lastStep = now;
            }
        } else {
            // keep time reference fresh while paused
            lastStep = steady_clock::now();
        }
        world.drawStatusLine(stdscr);
        int ch = getch();
        switch (ch) {
            case 'q': case 'Q':
                Logger::info("quit requested"); done = true; break;
            case 's': case 'S':
                world.toggleRunning(); Logger::info(std::string("running = ") + (world.isRunning()?"true":"false")); break;
            case 'p': case 'P':
                world.setRunning(false); Logger::info("paused"); break;
            case 'r': case 'R':
                world.reseedRandom(); world.drawAll(stdscr); Logger::info("reseed random"); break;
            case 'c': case 'C':
                world.setRunning(false); world.clear(); Logger::info("clear requested"); break;
            case '+':
                world.setStepDelayMs(world.getStepDelayMs() - 5); Logger::info("delay set(ms): " + std::to_string(world.getStepDelayMs())); break;
            case '-':
                world.setStepDelayMs(world.getStepDelayMs() + 5); Logger::info("delay set(ms): " + std::to_string(world.getStepDelayMs())); break;
            default:
                break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(16));
    }

    endwin();
    g_curses_inited = false;
    Logger::info("physics terminating");
    Logger::shutdown();
    return 0;
    } catch (const std::exception& e) {
        if (g_curses_inited) { endwin(); g_curses_inited = false; }
        Logger::logException("unhandled exception (physics)", e);
        Logger::shutdown();
        return 2;
    } catch (...) {
        if (g_curses_inited) { endwin(); g_curses_inited = false; }
        Logger::logUnknownException("unhandled exception (physics)");
        Logger::shutdown();
        return 2;
    }
}
