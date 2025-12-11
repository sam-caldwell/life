/**
 * @file PhysicsWorld.h
 * @brief Newtonian particle simulation rendered with ncurses (physics variant).
 */
#pragma once

#include <vector>
#include <mutex>
#include <atomic>
#include <random>
#include <utility>
#include <ncurses.h>

/**
 * @class PhysicsWorld
 * @brief Manages a set of particles with mass, velocity, gravity, collisions, and boundary bounces.
 */
class PhysicsWorld {
public:
    PhysicsWorld(int width, int height,
                 float gravity = 9.8f,
                 float radius = 1.0f,
                 float restitution = 0.98f);
    ~PhysicsWorld();

    // Control
    void setWindow(WINDOW* w);
    void setRunning(bool on);
    bool isRunning() const { return running.load(); }
    void toggleRunning() { setRunning(!isRunning()); }
    void setStepDelayMs(int ms);
    int  getStepDelayMs() const { return stepDelayMs.load(); }
    // Single-threaded stepping; caller invokes step(dt) from the UI loop.

    // World ops
    void clear();
    void reseed(unsigned count);
    void reseedRandom();

    // Rendering
    void drawAll(WINDOW* w);
    void drawStatusLine(WINDOW* w);

    // Color mapping for mass [0,100] -> color pair 1..16 (9..16 rendered with A_BOLD)
    int colorPairForMass(int m) const;

    // Info
    size_t particleCount() const;
    size_t maxParticles() const { return MaxParticles; }

    // Tunables
    void setGravity(float g);
    void setRadius(float r);
    void setRestitution(float e);
    float gravity() const { return gravityG; }
    float radius() const { return partRadius; }
    float restitution() const { return bounceRestitution; }

    // Advance simulation by dt seconds (single-threaded)
    void step(float dt);

private:
    struct Particle {
        unsigned long long id{0};
        float x{0}, y{0};
        float vx{0}, vy{0};
        int mass{1};           // [0,100]
        char symbol{'o'};      // display glyph
        int prev_ix{-1};       // previous drawn integer x
        int prev_iy{-1};       // previous drawn integer y
        bool alive{true};
        int elasticity10{10};  // 0..10 (0=inelastic, 10=elastic)
        int decayTicks{0};     // cycles counted for periodic decay
    };

    struct Cluster {
        double mass{0};
        double cx{0};
        double cy{0};
        std::vector<size_t> members; // indices of particles
    };

    // Thread loop
    void run(); // unused (no background thread)

    // Physics helpers
    void computeClusters(std::vector<Cluster>& out);
    void applyGravity(const std::vector<Cluster>& clusters, float dt);
    void handleBoundary(float padding);
    void handleCollisions();
    void updateDecay();
    void updateAdjacencyAndCombine();
    void updateParticle(struct Particle& p, float dts);

    // Drawing helpers (callers hold mtx)
    void drawParticleUnlocked(const Particle& p);
    void eraseParticleUnlocked(const Particle& p);

    int w, h; // drawable area (excluding status line)
    std::vector<Particle> particles;
    mutable std::mutex mtx;
    std::atomic<bool> running{false};
    std::atomic<int> stepDelayMs{25}; // ms per desired step cadence
    WINDOW* win{nullptr};

    // Randomness
    std::random_device rd;
    std::mt19937 prng;

    // Identity and adjacency
    unsigned long long nextId{1};
    struct PairKey { unsigned long long a; unsigned long long b; };
    struct PairKeyHash {
        size_t operator()(const PairKey& k) const noexcept {
            // 64-bit mix
            unsigned long long x = k.a * 0x9E3779B185EBCA87ULL ^ k.b;
            x ^= (x >> 33);
            x *= 0xff51afd7ed558ccdULL;
            x ^= (x >> 33);
            x *= 0xc4ceb9fe1a85ec53ULL;
            x ^= (x >> 33);
            return (size_t)x;
        }
    };
    struct PairKeyEq {
        bool operator()(const PairKey& x, const PairKey& y) const noexcept { return x.a==y.a && x.b==y.b; }
    };
    std::unordered_map<PairKey, int, PairKeyHash, PairKeyEq> adjacencyCounts; // consecutive frames adjacent

    // Constants
    float gravityG;               // gravitational constant
    float partRadius;             // particle radius (uniform)
    float bounceRestitution;      // near-elastic boundary bounces
    static constexpr size_t MaxParticles = 100; // hard cap per requirements

    // Runtime clock (ms spent while running)
    long long runAccumMs{0};
};
