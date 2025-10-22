/**
 * @file Board.h
 * @brief Declares the Board class which owns the simulation grid, automata, rendering, and thread-safe interactions.
 *
 * The Board centralizes all mutable state for the simulation. All Automaton threads request actions through
 * these APIs; Board applies them atomically under an internal mutex and updates the ncurses display incrementally.
 *
 * @copyright Copyright (c) 2025 Sam Caldwell. Released under the MIT License.
 */
#pragma once

#include <vector>
#include <memory>
#include <mutex>
#include <atomic>
#include <random>
#include <unordered_map>
#include <ncurses.h>
#include <functional>

class Automaton;

/**
 * @class Board
 * @brief Thread-safe owner of the simulation world and incremental ncurses renderer.
 *
 * Responsibilities:
 * - Grid storage and occupancy checks
 * - Automata lifecycle (seeding, removal, capped population)
 * - Conflict-resolved interactions (push, eat, spawn)
 * - Neighbor queries and position lookups
 * - Incremental drawing (cells + status line)
 */
class Board {
public:
    /** @brief Construct a board of given width/height measured in character cells. */
    Board(int width, int height);
    /** @brief Destructor; stops and joins all automata and clears the grid. */
    ~Board();

    /** @brief Hard cap on the number of concurrent automata (1000 per CPU core, computed at runtime). */
    size_t maxAutomata() const { return maxAutomataCap; }

    // Simulation control
    /** @brief Set the simulation running (true) or paused (false). */
    void setRunning(bool on);
    /** @brief Query whether the simulation is running. */
    bool isRunning() const { return running.load(); }
    /** @brief Toggle running/paused. */
    void toggleRunning() { setRunning(!isRunning()); }

    /** @brief Set per-automaton sleep delay in milliseconds; clamped to [5,2000]. */
    void setStepDelayMs(int ms);
    /** @brief Get the current per-automaton sleep delay in milliseconds. */
    int getStepDelayMs() const { return stepDelayMs.load(); }

    // Automata management
    /** @brief Stop all automata, join threads, clear the grid and positions, and erase display. */
    void clear();
    /** @brief Reset the board and randomly place @p count automata with random species/weights. */
    void reseed(unsigned count);

    // Grid helpers
    /** @brief Grid width in characters. */
    int width() const { return w; }
    /** @brief Grid height in characters. */
    int height() const { return h; }
    /** @brief Check if coordinates are within the grid. */
    bool inBounds(int x, int y) const { return x >= 0 && y >= 0 && x < w && y < h; }

    // Called by Automaton threads for interactions (actor calls neighbor->X(actor))
    /** @brief Attempt to push @p target one step away from @p actor, applying recoil rules if blocked. */
    bool handlePush(const std::shared_ptr<Automaton>& actor,
                    const std::shared_ptr<Automaton>& target);
    /** @brief Attempt to have @p actor eat @p target if rules allow; removes target and transfers weight. */
    bool handleEat(const std::shared_ptr<Automaton>& actor,
                   const std::shared_ptr<Automaton>& target);
    /** @brief Attempt to spawn a new automaton adjacent to @p actor; respects cap and load-based probability. */
    bool handleSpawn(const std::shared_ptr<Automaton>& actor,
                     const std::shared_ptr<Automaton>& partner);

    // Called by Automaton to query neighbors/position atomically
    /** @brief Legacy neighbor fetch (8-neighborhood); prefer neighborQuery for richer info. */
    void getNeighbors(int x, int y, std::vector<std::shared_ptr<Automaton>>& out);

    // Rendering
    /** @brief Full-screen initial draw of the grid contents (expensive); used once at startup. */
    void draw(WINDOW* win);
    /** @brief Update the bottom status line (legend, controls, thread count); cursor is kept at bottom. */
    void drawStatusLine(WINDOW* win);
    /** @brief Map a weight to a color pair id for ncurses rendering. */
    int colorPairForWeight(int w) const;
    /** @brief Assign the window used for incremental drawing. */
    void setWindow(WINDOW* w) { win = w; }
    /** @brief Draw a single cell (thread-safe) and refresh. */
    void drawCell(int x, int y);
    /** @brief Refresh the window (thread-safe). */
    void refresh();

    // Position accessors used by Automaton; these are guarded internally
    /** @brief Look up the current coordinates of @p who. */
    bool tryGetPosition(const std::shared_ptr<Automaton>& who, int& x, int& y);
    /** @brief Request a one-step cardinal move by delta (dx,dy); succeeds only if destination is free/in-bounds. */
    bool moveAutomaton(const std::shared_ptr<Automaton>& who, int dx, int dy);

    /**
     * @struct NeighborInfo
     * @brief Rich neighbor metadata returned by neighborQuery.
     */
    struct NeighborInfo {
        int nx{0};            /**< neighbor x coordinate */
        int ny{0};            /**< neighbor y coordinate */
        int dx{0};            /**< neighbor offset x relative to actor (nx-ax) */
        int dy{0};            /**< neighbor offset y relative to actor (ny-ay) */
        bool occupied{false}; /**< whether a neighbor exists at (nx,ny) */
        std::shared_ptr<Automaton> ref{}; /**< shared_ptr to neighbor (if occupied) */
        bool actorCanEat{false}; /**< true if the actor could eat this neighbor under rules */
        bool canEatActor{false}; /**< true if this neighbor could eat the actor under rules */
        int weight{0};        /**< neighbor weight (if occupied) */
        char symbol{' '};     /**< neighbor species letter (if occupied) */
    };
    /** @brief Populate @p out with metadata for the 8 neighbors of @p actor; also returns actor position in ax,ay. */
    bool neighborQuery(const std::shared_ptr<Automaton>& actor,
                       std::vector<NeighborInfo>& out,
                       int& ax, int& ay);

    // Info
    /** @brief Count of live automaton threads (best-effort snapshot). */
    size_t threadCount() const;

    // Random helpers
    /** @brief Uniform integer in [lo,hi]. */
    int randInt(int lo, int hi);
    /** @brief Access the board PRNG (mt19937). */
    std::mt19937& rng() { return prng; }
    /** @brief Uniform real in [0,1). */
    double rand01();
    /** @brief Current load as fraction of maxAutomata() in [0,1]. */
    double loadFactor() const;

    // Movement helpers
    /** @brief Attempt to flee a threat by moving one step to the best empty cardinal cell (internal heuristic). */
    bool attemptFlee(const std::shared_ptr<Automaton>& a);
    /** @brief Attempt to move toward edible prey (internal heuristic). */
    bool attemptPursue(const std::shared_ptr<Automaton>& a);
    /** @brief Remove an automaton from the grid and from the registry (thread-safe). */
    void removeAutomaton(const std::shared_ptr<Automaton>& a);

    /** @brief Rule helper: whether @p actor can eat @p target under species/weight rules. */
    bool canEat(const std::shared_ptr<Automaton>& actor,
                const std::shared_ptr<Automaton>& target) const;

private:
    /** @brief Grid cell: holds a weak reference to the occupant (if any). */
    struct Cell {
        std::weak_ptr<Automaton> occ; /**< occupant weak_ptr */
    };

    /** @brief Randomly place up to @p count automata on empty cells, assigning random species/weights. */
    void placeInitial(unsigned count);
    /** @brief Place an automaton @p a at (x,y) if empty; also registers it. */
    bool placeAt(const std::shared_ptr<Automaton>& a, int x, int y);
    /** @brief Remove @p a from grid and lists; caller must hold mtx. */
    void removeLocked(const std::shared_ptr<Automaton>& a);
    /** @brief Move @p a to (nx,ny) if empty; updates drawing; caller must hold mtx. */
    bool moveLocked(const std::shared_ptr<Automaton>& a, int nx, int ny);
    /** @brief Find a random empty neighbor (8-neighborhood) around (cx,cy). */
    bool findRandomEmptyAdjacent(int cx, int cy, int& ox, int& oy);
    /** @brief Draw a single cell without locking (caller holds mtx). */
    void drawCellUnlocked(int x, int y);

    // Internals
    int w, h; /**< grid dimensions */
    std::vector<Cell> grid; /**< flattened grid storage of size w*h */
    std::unordered_map<Automaton*, std::pair<int,int>> positions; /**< fast position map */
    WINDOW* win{nullptr}; /**< ncurses window for drawing */
    std::vector<std::shared_ptr<Automaton>> automata; /**< owned automata */
    mutable std::mutex mtx; /**< global board mutex protecting state and drawing */
    std::atomic<bool> running{false}; /**< run/pause flag */
    std::atomic<bool> quitting{false}; /**< reserved for future use */
    std::atomic<int> stepDelayMs{150}; /**< per-automaton sleep delay (ms) */

    std::random_device rd; /**< entropy for PRNG seeding */
    std::mt19937 prng;     /**< board PRNG */

    size_t maxAutomataCap{1000}; /**< runtime cap: 1000 threads per CPU core */
};
