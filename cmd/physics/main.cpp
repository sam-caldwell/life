/**
 * @file main_physics.cpp
 * @brief Physics entry point for the Newtonian variant.
 */
#include <ncurses.h>
#include <chrono>
#include <thread>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cerrno>
#include <random>
#include "PhysicsWorld.h"

static volatile sig_atomic_t g_stop = 0;
static void handle_signal(int) { g_stop = 1; }

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
    struct sigaction sa{};
    sa.sa_handler = handle_signal;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGINT, &sa, nullptr);
    sigaction(SIGTERM, &sa, nullptr);

    initscr();
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
    if (gridH < 1 || gridW < 1) { endwin(); return 1; }

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
        world.reseed((unsigned)d(gen));
    }
    world.drawAll(stdscr);
    world.drawStatusLine(stdscr);
    world.setRunning(false);

    bool done = false;
    using namespace std::chrono;
    auto lastStep = steady_clock::now();
    while (!done) {
        if (g_stop) done = true;
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
                done = true; break;
            case 's': case 'S':
                world.toggleRunning(); break;
            case 'p': case 'P':
                world.setRunning(false); break;
            case 'r': case 'R':
                world.reseedRandom(); world.drawAll(stdscr); break;
            case 'c': case 'C':
                world.setRunning(false); world.clear(); break;
            case '+':
                world.setStepDelayMs(world.getStepDelayMs() - 5); break;
            case '-':
                world.setStepDelayMs(world.getStepDelayMs() + 5); break;
            default:
                break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(16));
    }

    endwin();
    return 0;
}
