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
#include "Board.h"

static volatile sig_atomic_t g_stop = 0;
static void handle_signal(int) { g_stop = 1; }

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
int main() {
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
    nodelay(stdscr, TRUE); // non-blocking getch
    timeout(0);
    init_colors();

    int rows, cols;
    getmaxyx(stdscr, rows, cols);
    int gridH = rows - 1;
    int gridW = cols;
    if (gridH < 1 || gridW < 1) {
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

    bool done = false;
    while (!done) {
        if (g_stop) done = true;
        // Update status line only (incremental board updates happen per automaton)
        board.drawStatusLine(stdscr);

        int ch = getch();
        switch (ch) {
            case 'q':
            case 'Q':
                board.setStatusNote("terminating");
                board.drawStatusLine(stdscr);
                done = true; break;
            case 's': case 'S':
                board.toggleRunning(); break;
            case 'p': case 'P':
                board.setRunning(false); break;
            case 'r': case 'R': {
                unsigned cnt = (unsigned)((gridW * gridH) / 10);
                board.reseed(cnt);
                break; }
            case 'c': case 'C':
                board.clear();
                break;
            case '+': {
                board.setStepDelayMs(board.getStepDelayMs() - 10);
                break; }
            case '-': {
                board.setStepDelayMs(board.getStepDelayMs() + 10);
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
    return 0;
}
