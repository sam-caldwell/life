/**
 * @file Board.cpp
 * @brief Board implementation: state management, interactions, and incremental ncurses rendering.
 *
 * @copyright Copyright (c) 2025 Sam Caldwell. Released under the MIT License.
 */
#include "Board.h"
#include "Automaton.h"
#include "Logger.h"

#include <algorithm>
#include <chrono>
#include <array>
#include <thread>

namespace {
/** @brief Compute flattened index into grid vector for coordinates (x,y) in a width w grid. */
inline int idx(int x, int y, int w) { return y * w + x; }
}

/** @copydoc Board::Board */
Board::Board(int width, int height)
    : w(width), h(height), grid(static_cast<size_t>(width * height)), prng(rd()) {
    unsigned int cores = std::thread::hardware_concurrency();
    if (cores == 0) cores = 1;
    maxAutomataCap = static_cast<size_t>(cores) * 200ULL;
}

/** @copydoc Board::~Board */
Board::~Board() {
    quitting.store(true);
    Logger::info("Board destructor: initiating clear and shutdown");
    clear();
}

/** @copydoc Board::setRunning */
void Board::setRunning(bool on) {
    running.store(on);
}

/** @copydoc Board::setStepDelayMs */
void Board::setStepDelayMs(int ms) {
    if (ms < 5) ms = 5;
    if (ms > 2000) ms = 2000;
    stepDelayMs.store(ms);
}

/** @copydoc Board::clear */
void Board::clear() {
    std::vector<std::shared_ptr<Automaton>> toJoin;
    {
        std::lock_guard<std::mutex> lock(mtx);
        for (auto& a : automata) {
            if (a) a->requestStop();
        }
        Logger::info("Board::clear: stopping automata count=" + std::to_string(automata.size()));
        toJoin = automata;
        automata.clear();
        for (auto& c : grid) c.occ.reset();
        positions.clear();
        if (win) {
            for (int y = 0; y < h; ++y) {
                for (int x = 0; x < w; ++x) mvwaddch(win, y, x, ' ');
            }
            wrefresh(win);
        }
    }
    for (auto& a : toJoin) if (a) a->join();
}

/** @copydoc Board::reseed */
void Board::reseed(unsigned count) {
    clear();
    placeInitial(count);
    Logger::info("Board::reseed: placed initial automata count=" + std::to_string(count));
}

/** @brief See Board::reseed; helper to place initial automata with random positions and traits. */
void Board::placeInitial(unsigned count) {
    std::uniform_int_distribution<int> weightDist(10, 50);
    // Weighted species distribution: A most likely, Z least likely
    auto pickSpeciesWeighted = [this]() -> char {
        static const std::array<int,26> weights = []{
            std::array<int,26> w{}; // A..Z => 26..1
            for (int i=0;i<26;++i) w[i] = 26 - i;
            return w;
        }();
        static const int total = []{
            int s=0; for (int i=0;i<26;++i) s += 26 - i; return s; // 351
        }();
        std::uniform_int_distribution<int> dist(1, total);
        int r = dist(prng);
        int acc = 0;
        for (int i=0;i<26;++i) { acc += weights[i]; if (r <= acc) return static_cast<char>('A' + i); }
        return 'Z';
    };

    std::lock_guard<std::mutex> lock(mtx);
    int total = w * h;
    if (count > static_cast<unsigned>(total)) count = static_cast<unsigned>(total);
    if (count > static_cast<unsigned>(maxAutomataCap)) count = static_cast<unsigned>(maxAutomataCap);

    // Random unique positions
    std::vector<int> posIdx(total);
    for (int i = 0; i < total; ++i) posIdx[i] = i;
    std::shuffle(posIdx.begin(), posIdx.end(), prng);

    // Ensure at least two Z automata if possible
    unsigned z_needed = (count >= 2) ? 2u : count;
    unsigned idxPos = 0;
    for (unsigned i = 0; i < z_needed; ++i) {
        int p = posIdx[idxPos++];
        int x = p % w;
        int y = p / w;
        char sym = 'Z';
        short color = static_cast<short>(1 + (sym - 'A') % 7);
        auto a = std::make_shared<Automaton>(*this, sym, color, weightDist(prng));
        automata.push_back(a);
        grid[static_cast<size_t>(p)].occ = a;
        positions[a.get()] = {x, y};
        if (win) { drawCellUnlocked(x, y); }
        try {
            a->start();
        } catch (const std::exception& e) {
            Logger::logException("Board::placeInitial start(Z) failed", e);
            grid[static_cast<size_t>(p)].occ.reset();
            positions.erase(a.get());
            automata.pop_back();
            if (win) { drawCellUnlocked(x, y); wrefresh(win); }
            continue;
        } catch (...) {
            Logger::logUnknownException("Board::placeInitial start(Z) failed");
            grid[static_cast<size_t>(p)].occ.reset();
            positions.erase(a.get());
            automata.pop_back();
            if (win) { drawCellUnlocked(x, y); wrefresh(win); }
            continue;
        }
    }

    for (; idxPos < count; ++idxPos) {
        int p = posIdx[idxPos];
        int x = p % w;
        int y = p / w;
        char sym = pickSpeciesWeighted();
        short color = static_cast<short>(1 + (sym - 'A') % 7);
        auto a = std::make_shared<Automaton>(*this, sym, color, weightDist(prng));
        automata.push_back(a);
        grid[static_cast<size_t>(p)].occ = a;
        positions[a.get()] = {x, y};
        if (win) { drawCellUnlocked(x, y); }
        try {
            a->start();
        } catch (const std::exception& e) {
            Logger::logException("Board::placeInitial start failed", e);
            grid[static_cast<size_t>(p)].occ.reset();
            positions.erase(a.get());
            automata.pop_back();
            if (win) { drawCellUnlocked(x, y); wrefresh(win); }
            continue;
        } catch (...) {
            Logger::logUnknownException("Board::placeInitial start failed");
            grid[static_cast<size_t>(p)].occ.reset();
            positions.erase(a.get());
            automata.pop_back();
            if (win) { drawCellUnlocked(x, y); wrefresh(win); }
            continue;
        }
    }
}

/** @brief See Board::placeAt */
bool Board::placeAt(const std::shared_ptr<Automaton>& a, int x, int y) {
    if (!inBounds(x, y)) return false;
    auto& cell = grid[static_cast<size_t>(idx(x, y, w))];
    if (!cell.occ.expired()) return false;
    cell.occ = a;
    automata.push_back(a);
    return true;
}

/** @brief See Board::removeAutomaton; internal locked removal. */
void Board::removeLocked(const std::shared_ptr<Automaton>& a) {
    // Remove from grid
    auto it = positions.find(a.get());
    if (it != positions.end()) {
        int x = it->second.first;
        int y = it->second.second;
        auto& c = grid[static_cast<size_t>(idx(x,y,w))];
        c.occ.reset();
        positions.erase(it);
        if (win) { drawCellUnlocked(x, y); wrefresh(win); }
    } else {
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                auto& c = grid[static_cast<size_t>(idx(x,y,w))];
                auto sp = c.occ.lock();
                if (sp == a) { c.occ.reset(); if (win) { drawCellUnlocked(x, y); wrefresh(win);} }
            }
        }
    }
    // Remove from list
    automata.erase(std::remove(automata.begin(), automata.end(), a), automata.end());
}

/** @brief Internal locked move to absolute destination; updates positions and repaints. */
bool Board::moveLocked(const std::shared_ptr<Automaton>& a, int nx, int ny) {
    if (!inBounds(nx, ny)) return false;
    auto& dest = grid[static_cast<size_t>(idx(nx, ny, w))];
    if (!dest.occ.expired()) return false;

    auto it = positions.find(a.get());
    if (it == positions.end()) return false;
    int x = it->second.first;
    int y = it->second.second;
    auto& src = grid[static_cast<size_t>(idx(x, y, w))];
    auto sp = src.occ.lock();
    if (sp != a) return false;

    // move
    src.occ.reset();
    dest.occ = a;
    positions[a.get()] = {nx, ny};
    if (win) {
        drawCellUnlocked(x, y);
        drawCellUnlocked(nx, ny);
        wrefresh(win);
    }
    return true;
}

/** @brief See Board::handleSpawn; picks a random empty neighbor cell around (cx,cy). */
bool Board::findRandomEmptyAdjacent(int cx, int cy, int& ox, int& oy) {
    std::vector<std::pair<int,int>> spots;
    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) continue;
            int nx = cx + dx;
            int ny = cy + dy;
            if (!inBounds(nx, ny)) continue;
            auto& cell = grid[static_cast<size_t>(idx(nx, ny, w))];
            if (cell.occ.expired()) spots.emplace_back(nx, ny);
        }
    }
    if (spots.empty()) return false;
    std::uniform_int_distribution<size_t> dist(0, spots.size() - 1);
    auto p = spots[dist(prng)];
    ox = p.first; oy = p.second;
    return true;
}

/** @copydoc Board::handlePush */
bool Board::handlePush(const std::shared_ptr<Automaton>& actor,
                       const std::shared_ptr<Automaton>& target) {
    std::lock_guard<std::mutex> lock(mtx);
    if (!actor || !target) return false;

    // Locate positions via map
    auto ita = positions.find(actor.get());
    auto itt = positions.find(target.get());
    if (ita == positions.end() || itt == positions.end()) return false;
    int ax = ita->second.first, ay = ita->second.second;
    int tx = itt->second.first, ty = itt->second.second;

    if (actor->weight() <= target->weight()) return false;
    int dx = tx - ax;
    int dy = ty - ay;
    if (dx == 0 && dy == 0) return false;
    // move target one cardinal step directly away from actor
    int stepx = 0, stepy = 0;
    int adx = std::abs(dx), ady = std::abs(dy);
    if (adx > ady) { stepx = (dx > 0 ? 1 : -1); stepy = 0; }
    else if (ady > adx) { stepx = 0; stepy = (dy > 0 ? 1 : -1); }
    else { if (randInt(0,1) == 0) { stepx = (dx > 0 ? 1 : -1); stepy = 0; } else { stepx = 0; stepy = (dy > 0 ? 1 : -1); } }

    int nx = tx + stepx;
    int ny = ty + stepy;
    if (!inBounds(nx, ny)) return false;
    auto& dest = grid[static_cast<size_t>(idx(nx, ny, w))];
    auto occ = dest.occ.lock();
    if (!occ) {
        // Simple case: move pushed target one step
        bool ok = moveLocked(target, nx, ny);
        if (ok) Logger::debug("push: " + std::string(1, actor->symbol()) + "->" + std::string(1, target->symbol()));
        return ok;
    }

    // Destination occupied -> compute recoil rule
    int sum = actor->weight() + occ->weight();
    int tw = target->weight();
    if (sum > tw) {
        int k = sum - tw; // max recoil steps
        int cx = tx;
        int cy = ty;
        bool moved = false;
        for (int i = 0; i < k; ++i) {
            int rx = cx - stepx;
            int ry = cy - stepy;
            if (!inBounds(rx, ry)) break; // stop at boundary
            auto& rcell = grid[static_cast<size_t>(idx(rx, ry, w))];
            if (!rcell.occ.expired()) break; // stop at obstacle
            // perform one-step move
            if (!moveLocked(target, rx, ry)) break;
            moved = true;
            cx = rx; cy = ry;
        }
        if (moved) Logger::debug("push-recoil: " + std::string(1, actor->symbol()) + "->" + std::string(1, target->symbol()));
        return moved;
    }

    // Sum not greater than target weight -> push fails
    return false;
}

/** @copydoc Board::handleEat */
bool Board::handleEat(const std::shared_ptr<Automaton>& actor,
                      const std::shared_ptr<Automaton>& target) {
    if (!actor || !target) return false;
    bool success = false;
    {
        std::lock_guard<std::mutex> lock(mtx);

        // Verify adjacency
        int ax=-1, ay=-1, tx=-1, ty=-1;
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                auto sp = grid[static_cast<size_t>(idx(x,y,w))].occ.lock();
                if (sp == actor) { ax = x; ay = y; }
                if (sp == target) { tx = x; ty = y; }
            }
        }
        if (ax < 0 || tx < 0) return false;
        if (std::max(std::abs(ax - tx), std::abs(ay - ty)) > 1) return false;

        // Rules via canEat() (including same-species <10 exception)
        if (!canEat(actor, target)) return false;

        int newW = actor->weight() + target->weight();
        if (newW > 100) newW = 100;
        actor->setWeight(newW);

        // Remove target from grid and list
        removeLocked(target);
        // Update actor cell color for new weight
        auto itp = positions.find(actor.get());
        if (itp != positions.end() && win) {
            drawCellUnlocked(itp->second.first, itp->second.second);
            wrefresh(win);
        }
        success = true;
    }
    if (success) {
        target->requestStop();
        if (target->threadId() != std::this_thread::get_id()) {
            target->join();
        }
        Logger::info(std::string("eat: ") + actor->symbol() + " ate " + target->symbol() +
                     ", newW=" + std::to_string(actor->weight()));
    }
    return success;
}

/** @copydoc Board::handleSpawn */
bool Board::handleSpawn(const std::shared_ptr<Automaton>& actor,
                        const std::shared_ptr<Automaton>& partner) {
    std::lock_guard<std::mutex> lock(mtx);
    if (!actor || !partner) return false;

    if (automata.size() >= maxAutomataCap) return false;
    double lf = std::min(1.0, static_cast<double>(automata.size()) / static_cast<double>(maxAutomataCap));
    // Higher spawn probability when population is low; lower as it approaches the cap
    double accept = 0.2 + 0.8 * (1.0 - lf);
    if (rand01() > accept) return false;

    // Find actor position
    int ax=-1, ay=-1;
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            auto sp = grid[static_cast<size_t>(idx(x,y,w))].occ.lock();
            if (sp == actor) { ax = x; ay = y; break; }
        }
    }
    if (ax < 0) return false;

    int sx, sy;
    if (!findRandomEmptyAdjacent(ax, ay, sx, sy)) return false;

    // New automaton inherits actor species/color; random weight
    std::uniform_int_distribution<int> weightDist(10, 50);
    auto a = std::make_shared<Automaton>(*this, actor->symbol(), actor->color(), weightDist(prng));
    automata.push_back(a);
    grid[static_cast<size_t>(idx(sx, sy, w))].occ = a;
    positions[a.get()] = {sx, sy};
    if (win) { drawCellUnlocked(sx, sy); wrefresh(win); }
    try {
        a->start();
    } catch (const std::exception& e) {
        Logger::logException("Board::handleSpawn start failed", e);
        grid[static_cast<size_t>(idx(sx, sy, w))].occ.reset();
        positions.erase(a.get());
        automata.pop_back();
        if (win) { drawCellUnlocked(sx, sy); wrefresh(win); }
        return false;
    } catch (...) {
        Logger::logUnknownException("Board::handleSpawn start failed");
        grid[static_cast<size_t>(idx(sx, sy, w))].occ.reset();
        positions.erase(a.get());
        automata.pop_back();
        if (win) { drawCellUnlocked(sx, sy); wrefresh(win); }
        return false;
    }
    Logger::info(std::string("spawn: ") + actor->symbol() + " + " + partner->symbol() +
                 " -> " + a->symbol() + " at (" + std::to_string(sx) + "," + std::to_string(sy) + ")");
    return true;
}

/** @copydoc Board::getNeighbors */
void Board::getNeighbors(int x, int y, std::vector<std::shared_ptr<Automaton>>& out) {
    std::lock_guard<std::mutex> lock(mtx);
    out.clear();
    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) continue;
            int nx = x + dx;
            int ny = y + dy;
            if (!inBounds(nx, ny)) continue;
            auto sp = grid[static_cast<size_t>(idx(nx, ny, w))].occ.lock();
            if (sp) out.push_back(sp);
        }
    }
}

/** @copydoc Board::draw */
void Board::draw(WINDOW* win) {
    std::lock_guard<std::mutex> lock(mtx);
    werase(win);
    // Draw occupants
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            auto sp = grid[static_cast<size_t>(idx(x,y,w))].occ.lock();
            if (sp) {
                int wpair = colorPairForWeight(sp->weight());
                bool bold = (sp->weight() >= 100);
                if (bold) wattron(win, A_BOLD);
                wattron(win, COLOR_PAIR(wpair));
                mvwaddch(win, y, x, sp->symbol());
                wattroff(win, COLOR_PAIR(wpair));
                if (bold) wattroff(win, A_BOLD);
            } else {
                // empty
                mvwaddch(win, y, x, ' ');
            }
        }
    }
    wrefresh(win);
}

/** @copydoc Board::drawStatusLine */
void Board::drawStatusLine(WINDOW* win) {
    std::lock_guard<std::mutex> lock(mtx);
    int rows, cols; getmaxyx(win, rows, cols);
    // Ensure cursor at bottom-left and clear the status line only
    wmove(win, rows - 1, 0);
    wclrtoeol(win);
    // Legend for weight spectrum
    const int thresholds[9] = {0,10,20,30,50,70,80,90,100};
    int x = 0;
    mvwprintw(win, rows - 1, x, "W:"); x += 2;
    for (int i = 0; i < 9; ++i) {
        int t = thresholds[i];
        int pair = colorPairForWeight(t);
        bool bold = (t >= 100);
        if (bold) wattron(win, A_BOLD);
        wattron(win, COLOR_PAIR(pair));
        // Print the label with padding
        char buf[8];
        snprintf(buf, sizeof(buf), "%d", t);
        mvwprintw(win, rows - 1, x, "%s", buf);
        x += (int)strlen(buf);
        wattroff(win, COLOR_PAIR(pair));
        if (bold) wattroff(win, A_BOLD);
        if (i != 8) { mvwaddch(win, rows - 1, x, ' '); ++x; }
    }

    // Controls and status
    // Compute live thread count while holding lock
    size_t live = 0; for (auto& a : automata) if (a && a->isAlive()) ++live;
    std::string status = "  | [s]tart/[p]ause  [r]eseed  [c]lear  speed[-/+]  [q]uit  | Threads: ";
    status += std::to_string(live);
    status += "/";
    status += std::to_string(maxAutomataCap);
    status += "  Delay(ms): ";
    status += std::to_string(getStepDelayMs());
    status += isRunning() ? "  | RUNNING" : "  | PAUSED";
    if (!statusNote.empty()) { status += "  | "; status += statusNote; }
    if ((int)status.size() + x < cols) status.append((size_t)(cols - x - (int)status.size()), ' ');
    mvwprintw(win, rows - 1, x, "%s", status.c_str());
    wrefresh(win);
}

/** @copydoc Board::tryGetPosition */
bool Board::tryGetPosition(const std::shared_ptr<Automaton>& who, int& x, int& y) {
    std::lock_guard<std::mutex> lock(mtx);
    auto it = positions.find(who.get());
    if (it == positions.end()) return false;
    x = it->second.first;
    y = it->second.second;
    return true;
}

/** @copydoc Board::threadCount */
size_t Board::threadCount() const {
    std::lock_guard<std::mutex> lock(mtx);
    size_t live = 0;
    for (auto& a : automata) if (a && a->isAlive()) ++live;
    return live;
}

/** @copydoc Board::randInt */
int Board::randInt(int lo, int hi) {
    std::uniform_int_distribution<int> d(lo, hi);
    return d(prng);
}

/** @copydoc Board::moveAutomaton */
bool Board::moveAutomaton(const std::shared_ptr<Automaton>& who, int dx, int dy) {
    std::lock_guard<std::mutex> lock(mtx);
    auto it = positions.find(who.get());
    if (it == positions.end()) return false;
    int nx = it->second.first + dx;
    int ny = it->second.second + dy;
    return moveLocked(who, nx, ny);
}

/** @copydoc Board::neighborQuery */
bool Board::neighborQuery(const std::shared_ptr<Automaton>& actor,
                          std::vector<NeighborInfo>& out,
                          int& ax, int& ay) {
    std::lock_guard<std::mutex> lock(mtx);
    auto it = positions.find(actor.get());
    if (it == positions.end()) return false;
    ax = it->second.first;
    ay = it->second.second;
    out.clear();
    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) continue;
            int nx = ax + dx;
            int ny = ay + dy;
            if (!inBounds(nx, ny)) continue;
            NeighborInfo info;
            info.nx = nx; info.ny = ny; info.dx = dx; info.dy = dy;
            auto sp = grid[static_cast<size_t>(idx(nx, ny, w))].occ.lock();
            info.occupied = (bool)sp;
            info.ref = sp;
            if (sp) {
                info.weight = sp->weight();
                info.symbol = sp->symbol();
                info.actorCanEat = canEat(actor, sp);
                info.canEatActor = canEat(sp, actor);
            }
            out.push_back(std::move(info));
        }
    }
    return true;
}

int Board::colorPairForWeight(int wgt) const {
    if (wgt <= 0) return 1;      // black
    if (wgt <= 10) return 2;     // blue
    if (wgt <= 20) return 3;     // green
    if (wgt <= 30) return 4;     // cyan
    if (wgt <= 50) return 5;     // yellow
    if (wgt <= 70) return 6;     // magenta
    if (wgt <= 80) return 7;     // red
    return 8;                    // white (90..100)
}

/** @copydoc Board::rand01 */
double Board::rand01() {
    std::uniform_real_distribution<double> d(0.0, 1.0);
    return d(prng);
}

/** @copydoc Board::loadFactor */
double Board::loadFactor() const {
    std::lock_guard<std::mutex> lock(mtx);
    double denom = static_cast<double>(maxAutomataCap > 0 ? maxAutomataCap : 1);
    double lf = static_cast<double>(automata.size()) / denom;
    if (lf < 0.0) lf = 0.0; if (lf > 1.0) lf = 1.0; return lf;
}

void Board::drawCellUnlocked(int x, int y) {
    if (!win) return;
    auto sp = grid[static_cast<size_t>(idx(x,y,w))].occ.lock();
    if (sp) {
        int wpair = colorPairForWeight(sp->weight());
        bool bold = (sp->weight() >= 100);
        if (bold) wattron(win, A_BOLD);
        wattron(win, COLOR_PAIR(wpair));
        mvwaddch(win, y, x, sp->symbol());
        wattroff(win, COLOR_PAIR(wpair));
        if (bold) wattroff(win, A_BOLD);
    } else {
        mvwaddch(win, y, x, ' ');
    }
}

/** @copydoc Board::attemptFlee */
bool Board::attemptFlee(const std::shared_ptr<Automaton>& a) {
    std::lock_guard<std::mutex> lock(mtx);
    if (!a) return false;
    int ax=-1, ay=-1;
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            auto sp = grid[static_cast<size_t>(idx(x,y,w))].occ.lock();
            if (sp == a) { ax = x; ay = y; break; }
        }
    }
    if (ax < 0) return false;

    // Identify predator neighbors (8-neighborhood) and candidate empty cardinal cells
    std::vector<std::pair<int,int>> heavy;
    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) continue;
            int nx = ax + dx, ny = ay + dy;
            if (!inBounds(nx, ny)) continue;
            auto sp = grid[static_cast<size_t>(idx(nx,ny,w))].occ.lock();
            if (sp && canEat(sp, a)) heavy.emplace_back(nx, ny);
        }
    }
    if (heavy.empty()) return false;

    std::vector<std::pair<int,int>> candidates;
    const std::pair<int,int> dirs[4]{{1,0},{-1,0},{0,1},{0,-1}};
    for (auto [dx,dy] : dirs) {
        int nx = ax + dx, ny = ay + dy;
        if (!inBounds(nx, ny)) continue;
        auto& cell = grid[static_cast<size_t>(idx(nx,ny,w))];
        if (cell.occ.expired()) candidates.emplace_back(nx, ny);
    }
    if (candidates.empty()) return false;

    // Choose candidate that maximizes min Chebyshev distance to heavier neighbors
    auto cheb = [](int x1,int y1,int x2,int y2){ return std::max(std::abs(x1-x2), std::abs(y1-y2)); };
    int bestIdx = -1; int bestScore = -1;
    for (size_t i = 0; i < candidates.size(); ++i) {
        auto [cx, cy] = candidates[i];
        int minDist = 999999;
        for (auto& hn : heavy) {
            minDist = std::min(minDist, cheb(cx, cy, hn.first, hn.second));
        }
        if (minDist > bestScore) { bestScore = minDist; bestIdx = (int)i; }
    }
    if (bestIdx < 0) return false;
    auto [mx, my] = candidates[(size_t)bestIdx];
    return moveLocked(a, mx, my);
}

/** @copydoc Board::removeAutomaton */
void Board::removeAutomaton(const std::shared_ptr<Automaton>& a) {
    if (!a) return;
    {
        std::lock_guard<std::mutex> lock(mtx);
        removeLocked(a);
    }
    // Stop and join outside the lock when possible
    a->requestStop();
    if (a->threadId() != std::this_thread::get_id()) {
        a->join();
    }
}

/** @copydoc Board::attemptPursue */
bool Board::attemptPursue(const std::shared_ptr<Automaton>& a) {
    std::lock_guard<std::mutex> lock(mtx);
    if (!a) return false;
    int ax=-1, ay=-1;
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            auto sp = grid[static_cast<size_t>(idx(x,y,w))].occ.lock();
            if (sp == a) { ax = x; ay = y; break; }
        }
        if (ax >= 0) break;
    }
    if (ax < 0) return false;

    // Identify prey neighbors (8-neighborhood)
    std::vector<std::pair<int,int>> prey;
    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) continue;
            int nx = ax + dx, ny = ay + dy;
            if (!inBounds(nx, ny)) continue;
            auto sp = grid[static_cast<size_t>(idx(nx,ny,w))].occ.lock();
            if (sp && canEat(a, sp)) prey.emplace_back(nx, ny);
        }
    }
    if (prey.empty()) return false;

    auto manh = [](int x1,int y1,int x2,int y2){ return std::abs(x1-x2)+std::abs(y1-y2); };
    int currentMin = INT_MAX;
    for (auto& p : prey) currentMin = std::min(currentMin, manh(ax, ay, p.first, p.second));

    // Candidate moves: cardinal empty cells that reduce distance to prey
    const std::pair<int,int> dirs[4]{{1,0},{-1,0},{0,1},{0,-1}};
    int bestx = ax, besty = ay, bestDist = currentMin;
    bool found = false;
    for (auto [dx,dy] : dirs) {
        int nx = ax + dx, ny = ay + dy;
        if (!inBounds(nx, ny)) continue;
        auto& cell = grid[static_cast<size_t>(idx(nx,ny,w))];
        if (!cell.occ.expired()) continue;
        int minD = INT_MAX;
        for (auto& p : prey) minD = std::min(minD, manh(nx, ny, p.first, p.second));
        if (minD < bestDist) { bestDist = minD; bestx = nx; besty = ny; found = true; }
    }
    if (!found) return false;
    return moveLocked(a, bestx, besty);
}

/** @copydoc Board::canEat */
bool Board::canEat(const std::shared_ptr<Automaton>& actor,
                   const std::shared_ptr<Automaton>& target) const {
    if (!actor || !target) return false;
    char as = actor->symbol();
    char ts = target->symbol();
    if (as > ts) return true; // higher letter eats lower letter
    if (as == ts) {
        // Same-species cannibalism: allowed only if species A–R and actor's weight < 10
        if (as <= 'R' && actor->weight() < 10) return true;
        return false;
    }
    return false;
}
void Board::drawCell(int x, int y) {
    std::lock_guard<std::mutex> lock(mtx);
    if (!win) return;
    drawCellUnlocked(x, y);
    wrefresh(win);
}

void Board::refresh() {
    std::lock_guard<std::mutex> lock(mtx);
    if (win) wrefresh(win);
}

/** @copydoc Board::setStatusNote */
void Board::setStatusNote(const std::string& note) {
    std::lock_guard<std::mutex> lock(mtx);
    statusNote = note;
}

/** @copydoc Board::clearStatusNote */
void Board::clearStatusNote() {
    std::lock_guard<std::mutex> lock(mtx);
    statusNote.clear();
}
