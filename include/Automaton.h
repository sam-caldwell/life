/**
 * @file Automaton.h
 * @brief Declares the Automaton class: an autonomous, threaded life form acting on the Board.
 *
 * @copyright Copyright (c) 2025 Sam Caldwell. Released under the MIT License.
 */
#pragma once

#include <thread>
#include <atomic>
#include <memory>

class Board;

/**
 * @class Automaton
 * @brief One agent in the simulation. Each instance runs in its own thread and makes local decisions.
 *
 * Responsibilities:
 * - Decide to flee/pursue/eat/push/spawn based on neighbor info and internal state
 * - Request movements and actions through Board
 * - Track hunger, weight, spawn backoff, and mate-seeking memory
 */
class Automaton : public std::enable_shared_from_this<Automaton> {
public:
    /** @brief Construct an automaton bound to @p board with species @p symbol, display color, and starting weight. */
    Automaton(Board& board, char symbol, short colorPair, int weight);
    /** @brief Destructor; requests stop and joins the worker thread if needed. */
    ~Automaton();

    // Lifecycle
    /** @brief Start the automaton worker thread. */
    void start();
    /** @brief Ask the worker thread to stop at the next opportunity. */
    void requestStop();
    /** @brief Whether the worker is still running. */
    bool isAlive() const { return alive.load(); }
    /** @brief Join the worker thread if joinable. */
    void join();
    /** @brief Return the std::thread::id of the worker thread (default-constructed if not started). */
    std::thread::id threadId() const { return worker.get_id(); }

    // Identity
    /** @brief Species letter (A–Z). */
    char symbol() const { return sym; }
    /** @brief Current color pair identifier (used by the renderer). */
    short color() const { return colorPair; }
    /** @brief Current weight [0,100]. */
    int weight() const { return wt.load(); }
    /** @brief Add @p delta to weight (clamped by caller). */
    void addWeight(int delta) { wt.fetch_add(delta); }

    // Interactions invoked by neighbors (actor is the invoker)
    /** @brief Be pushed by @p actor according to Board rules. */
    bool push(const std::shared_ptr<Automaton>& actor);
    /** @brief Be eaten by @p actor (if allowed) and removed from the board. */
    bool eat(const std::shared_ptr<Automaton>& actor);
    /** @brief Spawn with @p actor if eligible; new automaton appears adjacent to @p actor. */
    bool spawn(const std::shared_ptr<Automaton>& actor);

private:
    /** @brief Worker loop: senses neighbors, decides actions, sleeps between steps. */
    void run();

    Board& board;                  /**< reference to the owning Board */
    char sym;                      /**< species letter */
    short colorPair;               /**< color pair (weight rendering handled by Board) */
    std::atomic<int> wt;           /**< current weight */
    std::atomic<bool> alive{true}; /**< liveness flag for the worker loop */
    std::thread worker;            /**< worker thread */

    // Hunger/behavior state
    int cyclesSinceEat{0};        /**< cycles since last successful eat */
    int cyclesAtMaxWeight{0};     /**< consecutive cycles at weight==100 */

    // Mate-seeking memory
    std::weak_ptr<Automaton> mateTarget; /**< remembered mate candidate */
    int mateMemoryTTL{0};                /**< remaining cycles to pursue remembered mate */
    static constexpr int MateMemoryMaxTTL = 30; /**< default TTL for mate memory */

    // Spawn backoff to avoid repeated attempts when no space
    int spawnBackoffTTL{0};              /**< cycles to wait before reattempting spawn */
    static constexpr int SpawnBackoffMaxTTL = 10; /**< max spawn backoff cycles */

public:
    /** @brief Set absolute weight. */
    void setWeight(int v) { wt.store(v); }
};
