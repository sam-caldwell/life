/**
 * @file Automaton.cpp
 * @brief Automaton implementation: threaded decision loop and interaction requests to Board.
 *
 * @copyright Copyright (c) 2025 Sam Caldwell. Released under the MIT License.
 */
#include "Automaton.h"
#include "Board.h"
#include "Logger.h"

#include <vector>
#include <chrono>
#include <thread>
#include <string>

/** @copydoc Automaton::Automaton */
Automaton::Automaton(Board& b, char symbol, short colorPair_, int weight)
    : board(b), sym(symbol), colorPair(colorPair_), wt(weight) {}

/** @copydoc Automaton::~Automaton */
Automaton::~Automaton() {
    requestStop();
    if (worker.joinable()) {
        // Avoid joining self; detach if needed
        if (std::this_thread::get_id() != worker.get_id()) worker.join();
        else worker.detach();
    }
}

/** @copydoc Automaton::start */
void Automaton::start() {
    worker = std::thread(&Automaton::run, this);
}

/** @copydoc Automaton::requestStop */
void Automaton::requestStop() {
    alive.store(false);
}

/** @copydoc Automaton::join */
void Automaton::join() {
    if (worker.joinable()) worker.join();
}

/** @copydoc Automaton::push */
bool Automaton::push(const std::shared_ptr<Automaton>& actor) {
    // actor pushes this automaton
    return board.handlePush(actor, shared_from_this());
}

/** @copydoc Automaton::eat */
bool Automaton::eat(const std::shared_ptr<Automaton>& actor) {
    return board.handleEat(actor, shared_from_this());
}

/** @copydoc Automaton::spawn */
bool Automaton::spawn(const std::shared_ptr<Automaton>& actor) {
    return board.handleSpawn(actor, shared_from_this());
}

/** @brief Worker loop: sense -> decide -> act -> sleep, until stopped or removed. */
void Automaton::run() {
    using namespace std::chrono;
    std::shared_ptr<Automaton> self;
    try {
        Logger::debug("automaton thread starting: sym=" + std::string(1, sym));
        std::vector<std::shared_ptr<Automaton>> neighbors;
        self = shared_from_this();

        while (isAlive()) {
            if (!board.isRunning()) {
                std::this_thread::sleep_for(10ms);
                continue;
            }

            int x=-1,y=-1;
            if (!board.tryGetPosition(shared_from_this(), x, y)) {
                requestStop();
                break;
            }

            std::vector<Board::NeighborInfo> ninfo;
            if (!board.neighborQuery(self, ninfo, x, y)) {
                requestStop();
                break;
            }

            bool fled = false;
            std::vector<Board::NeighborInfo> predators;
            for (auto& ni : ninfo) if (ni.occupied && ni.symbol != sym && ni.symbol > sym) predators.push_back(ni);
            double moveProb = std::max(0.0, std::min(1.0, (100.0 - (double)wt.load()) / 100.0));
            if (!predators.empty() && board.rand01() < moveProb) {
                auto& pred = predators[(size_t)board.randInt(0, (int)predators.size()-1)];
                int dx = -pred.dx;
                int dy = -pred.dy;
                int sdx = dx == 0 ? 0 : (dx > 0 ? 1 : -1);
                int sdy = dy == 0 ? 0 : (dy > 0 ? 1 : -1);
                int mdx = 0, mdy = 0;
                if (std::abs(dx) > std::abs(dy)) { mdx = sdx; mdy = 0; }
                else if (std::abs(dy) > std::abs(dx)) { mdx = 0; mdy = sdy; }
                else { if (board.randInt(0,1)==0) { mdx = sdx; mdy = 0; } else { mdx = 0; mdy = sdy; } }
                fled = board.moveAutomaton(self, mdx, mdy);
            }
            if (fled) {
                double lf = board.loadFactor();
                int period = std::max(1, 10 - static_cast<int>(8 * lf));
                int dec = 1 + static_cast<int>(2 * lf);
                cyclesSinceEat++;
                if (cyclesSinceEat >= period) {
                    wt.fetch_sub(dec);
                    cyclesSinceEat = 0;
                }
                if (wt.load() <= 0) {
                    board.removeAutomaton(self);
                    requestStop();
                    break;
                }
                std::this_thread::sleep_for(milliseconds(board.getStepDelayMs()));
                continue;
            }

            bool ateThisCycle = false;
            if (!ninfo.empty()) {
                int myW = wt.load();
            bool ateRecently = (cyclesSinceEat < 10);

            // Prefer eating valid prey
            std::vector<Board::NeighborInfo> prey;
            for (auto& ni : ninfo) if (ni.occupied && ni.actorCanEat) prey.push_back(ni);
            if (!prey.empty()) {
                auto n = prey[(size_t)board.randInt(0, (int)prey.size()-1)].ref;
                if (n && n->eat(self)) {
                    ateThisCycle = true;
                } else if (board.rand01() < moveProb) {
                    auto chosen = prey.back();
                    int dx = chosen.dx; int dy = chosen.dy; // toward prey
                    int sdx = dx == 0 ? 0 : (dx > 0 ? 1 : -1);
                    int sdy = dy == 0 ? 0 : (dy > 0 ? 1 : -1);
                    int mdx = 0, mdy = 0;
                    if (std::abs(dx) > std::abs(dy)) { mdx = sdx; mdy = 0; }
                    else if (std::abs(dy) > std::abs(dx)) { mdx = 0; mdy = sdy; }
                    else { if (board.randInt(0,1)==0) { mdx = sdx; mdy = 0; } else { mdx = 0; mdy = sdy; } }
                    (void)board.moveAutomaton(self, mdx, mdy);
                }
            } else {
                // Seek same species to spawn
                std::vector<Board::NeighborInfo> sameKind;
                for (auto& ni : ninfo) if (ni.occupied && ni.symbol == sym) sameKind.push_back(ni);
                // If adjacent same-kind and spawn-eligible, spawn (with backoff if no space)
                bool spawned = false;
                if (ateRecently) {
                    for (auto& ni : sameKind) {
                        if (std::abs(myW - ni.weight) <= 1 && spawnBackoffTTL == 0) {
                            bool hasEmpty = false;
                            for (auto& sj : ninfo) { if (!sj.occupied) { hasEmpty = true; break; } }
                            if (hasEmpty) {
                                if (ni.ref->spawn(self)) { spawned = true; break; }
                                else { spawnBackoffTTL = SpawnBackoffMaxTTL; }
                            } else {
                                spawnBackoffTTL = SpawnBackoffMaxTTL;
                            }
                        }
                    }
                }
                // Update mate memory from immediate neighbors if none stored or to refresh TTL
                if (!sameKind.empty()) {
                    auto chosen = sameKind[(size_t)board.randInt(0, (int)sameKind.size()-1)].ref;
                    mateTarget = chosen;
                    mateMemoryTTL = MateMemoryMaxTTL;
                } else if (mateMemoryTTL > 0) {
                    // Pursue remembered mate target across cycles
                    auto sp = mateTarget.lock();
                    if (sp) {
                        int mx, my;
                        if (board.tryGetPosition(sp, mx, my)) {
                            int dx = mx - x; int dy = my - y;
                            int sdx = dx == 0 ? 0 : (dx > 0 ? 1 : -1);
                            int sdy = dy == 0 ? 0 : (dy > 0 ? 1 : -1);
                            int mdx = 0, mdy = 0;
                            if (std::abs(dx) > std::abs(dy)) { mdx = sdx; mdy = 0; }
                            else if (std::abs(dy) > std::abs(dx)) { mdx = 0; mdy = sdy; }
                            else { if (board.randInt(0,1)==0) { mdx = sdx; mdy = 0; } else { mdx = 0; mdy = sdy; } }
                            if (board.rand01() < moveProb) (void)board.moveAutomaton(self, mdx, mdy);
                        } else {
                            // target no longer on board
                            mateMemoryTTL = 0; mateTarget.reset();
                        }
                    } else {
                        mateMemoryTTL = 0;
                    }
                    if (mateMemoryTTL > 0) mateMemoryTTL--;
                }
                if (!spawned) {
                    // Move toward same-kind to meet up
                    if (!sameKind.empty() && board.rand01() < moveProb) {
                        auto chosen = sameKind[(size_t)board.randInt(0, (int)sameKind.size()-1)];
                        int dx = chosen.dx; int dy = chosen.dy; // toward same species
                        int sdx = dx == 0 ? 0 : (dx > 0 ? 1 : -1);
                        int sdy = dy == 0 ? 0 : (dy > 0 ? 1 : -1);
                        int mdx = 0, mdy = 0;
                        if (std::abs(dx) > std::abs(dy)) { mdx = sdx; mdy = 0; }
                        else if (std::abs(dy) > std::abs(dx)) { mdx = 0; mdy = sdy; }
                        else { if (board.randInt(0,1)==0) { mdx = sdx; mdy = 0; } else { mdx = 0; mdy = sdy; } }
                        (void)board.moveAutomaton(self, mdx, mdy);
                    } else if (!ninfo.empty() && spawnBackoffTTL == 0) {
                        // Fallback: attempt spawn/push with a random neighbor (legacy behavior)
                        int idx = board.randInt(0, static_cast<int>(ninfo.size()) - 1);
                        auto ni = ninfo[(size_t)idx];
                        if (ni.occupied) {
                            int nbW = ni.weight;
                            if (ateRecently && std::abs(myW - nbW) <= 1) {
                                bool hasEmpty = false;
                                for (auto& sj : ninfo) { if (!sj.occupied) { hasEmpty = true; break; } }
                                if (hasEmpty) {
                                    if (!ni.ref->spawn(self)) spawnBackoffTTL = SpawnBackoffMaxTTL;
                                } else {
                                    spawnBackoffTTL = SpawnBackoffMaxTTL;
                                }
                            } else if (myW > nbW) {
                                ni.ref->push(self);
                            }
                        }
                    }
                }
            }
        }

        // If still hungry and didn't flee or eat, attempt to wander using a heuristic
        if (!fled && !ateThisCycle) {
            double moveProb2 = std::max(0.0, std::min(1.0, (100.0 - (double)wt.load()) / 100.0));
            if (board.rand01() < moveProb2) {
                struct Cand { int dx; int dy; double score; };
                std::vector<Cand> cands;
                for (auto& ni : ninfo) {
                    if (!ni.occupied && ((ni.dx == 0) ^ (ni.dy == 0))) {
                        cands.push_back({ni.dx, ni.dy, 0.0});
                    }
                }
                if (!cands.empty()) {
                    for (auto& c : cands) {
                        double s = 0.0;
                        int edibleCount = 0;
                int sameCount = 0;
                for (auto& ni : ninfo) {
                    int dot = c.dx * ni.dx + c.dy * ni.dy; // alignment with candidate direction
                    int align = (dot > 0 ? 2 : (dot == 0 ? 1 : 0));
                    if (ni.actorCanEat) { s += 2.0 * align; ++edibleCount; }
                    if (ni.canEatActor) { s -= 3.0 * align; }
                    else if (ni.occupied && ni.weight > wt.load()) { s -= 1.0 * align; }
                    if (ni.occupied && ni.symbol == sym) { s += 1.0 * align; ++sameCount; }
                }
                if (cyclesSinceEat >= 7) s += edibleCount * 0.5; // get bolder when hungry
                if (cyclesSinceEat < 10) s += sameCount * 1.0;   // bias toward same species to meet and spawn
                s += (board.rand01() - 0.5) * 0.2; // small jitter to break ties
                c.score = s;
                    }
                    // pick best by score
                    size_t best = 0;
                    for (size_t i = 1; i < cands.size(); ++i) if (cands[i].score > cands[best].score) best = i;
                    (void)board.moveAutomaton(self, cands[best].dx, cands[best].dy);
                }
            }
        }

        if (ateThisCycle) {
            cyclesSinceEat = 0;
        } else {
            double lf = board.loadFactor();
            int period = std::max(1, 10 - static_cast<int>(8 * lf));
            int dec = 1 + static_cast<int>(2 * lf);
            cyclesSinceEat++;
            if (cyclesSinceEat >= period) {
                wt.fetch_sub(dec);
                cyclesSinceEat = 0;
            }
        }

        if (wt.load() <= 0) {
            board.removeAutomaton(self);
            requestStop();
            break;
        }

        // Track time spent at max weight and die if > 100 cycles at 100
        if (wt.load() >= 100) {
            cyclesAtMaxWeight++;
        } else {
            cyclesAtMaxWeight = 0;
        }
        if (cyclesAtMaxWeight > 100) {
            board.removeAutomaton(self);
            requestStop();
            break;
        }

        // Decrement spawn backoff TTL
        if (spawnBackoffTTL > 0) --spawnBackoffTTL;

        std::this_thread::sleep_for(milliseconds(board.getStepDelayMs()));
        }
        Logger::debug("automaton thread exiting normally: sym=" + std::string(1, sym));
    } catch (const std::exception& e) {
        alive.store(false);
        try { if (self) board.removeAutomaton(self); } catch (...) {}
        try { board.setStatusNote(std::string("thread error: ") + e.what()); } catch (...) {}
        Logger::error(std::string("automaton exception: ") + e.what());
    } catch (...) {
        alive.store(false);
        try { if (self) board.removeAutomaton(self); } catch (...) {}
        try { board.setStatusNote("thread error: unknown"); } catch (...) {}
        Logger::error("automaton exception: unknown");
    }
}
