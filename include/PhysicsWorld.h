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
#include <cstdint>
#include <functional>
#include <cstdint>
#include <thread>
#include <condition_variable>
#include <atomic>

/**
 * @class PhysicsWorld
 * @brief Manages a set of particles with mass, velocity, gravity, collisions, and boundary interactions
 *        (elastic bounces; inelastic wall collisions may split a particle into 2–4 lower‑symbol children
 *        with equal mass/energy shares when capacity allows).
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
    // Start/stop background physics thread (fixed-dt stepping)
    void startThread();
    void stopThread();
    void notifyThread();
    // Snapshot access for renderer (triple-buffered, lock-free read)
    struct RenderItem { int16_t ix; int16_t iy; char sym; uint8_t mass; uint8_t alive; };
    struct RenderSnapshot {
        std::vector<RenderItem> items;
        std::vector<uint8_t> hl; // per-cell highlight: 0 none, 1 split (yellow), 2 collision (red)
        int w{0}; int h{0};
        std::atomic<uint64_t> seq{0};
    };
    bool getLatestSnapshot(const std::vector<RenderItem>*& items, const std::vector<uint8_t>*& hl,
                           uint64_t& seq, int& outW, int& outH) const;

private:
    // Physics helpers (SoA variants)
    void applyGravityBarnesHutSoA(float dt);
    void handleBoundarySoA(float padding);
    void handleCollisionsSoA();
    void updateParticleSoA(float dts);
    // SoA Phase 2
    // SoA variants (Phase 2)
    void updateDecaySoA();
    void updateAdjacencyAndCombineSoA();
    void maybeSpawnFromPoolsSoA();
    void applyFragmentationsSoAFromSchedule();
    void soaCompact();
    void reseedSoA(unsigned count);
    void accumulateLostMass(double dm);
    void accumulateLostEnergy(double de);

    int w, h; // drawable area (excluding status line)
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
    // Flat open-addressing hash map specialized for PairKey -> int
    class FlatAdjMap {
    public:
        void clear() { _size = 0; std::fill(_state.begin(), _state.end(), (uint8_t)0); }
        void reserve(size_t n) {
            size_t need = (size_t)((double)n * 2.0); // load factor ~0.5
            size_t cap = _capacity;
            if (cap >= need) return;
            cap = 1;
            while (cap < need) cap <<= 1;
            rehash(cap);
        }
        void increment(const PairKey& k) {
            int idx = findOrInsert(k);
            _vals[(size_t)idx] += 1;
        }
        bool get(const PairKey& k, int& out) const {
            if (_capacity == 0) return false;
            size_t cap = _capacity;
            uint64_t h = hashOf(k);
            size_t i = (size_t)h & (cap - 1);
            for (;;) {
                if (_state[i] == 0) return false;
                if (_state[i] == 1 && _keys[i].a == k.a && _keys[i].b == k.b) { out = _vals[i]; return true; }
                i = (i + 1) & (cap - 1);
            }
        }
        void set(const PairKey& k, int v) {
            int idx = findOrInsert(k);
            _vals[(size_t)idx] = v;
        }
        bool erase(const PairKey& k) {
            if (_capacity == 0) return false;
            size_t cap = _capacity;
            uint64_t h = hashOf(k);
            size_t i = (size_t)h & (cap - 1);
            for (;;) {
                if (_state[i] == 0) return false;
                if (_state[i] == 1 && _keys[i].a == k.a && _keys[i].b == k.b) {
                    _state[i] = 2; // tombstone
                    _keys[i] = {};
                    _vals[i] = 0;
                    --_size;
                    return true;
                }
                i = (i + 1) & (cap - 1);
            }
        }
        template <class F>
        void forEach(F&& f) {
            for (size_t i = 0; i < _capacity; ++i) if (_state[i] == 1) f(_keys[i], _vals[i]);
        }
        size_t size() const { return _size; }
    private:
        size_t _size{0};
        size_t _capacity{0};
        std::vector<PairKey> _keys;
        std::vector<int> _vals;
        std::vector<uint8_t> _state; // 0 empty, 1 used, 2 tomb
        static inline uint64_t hashOf(const PairKey& k) {
            unsigned long long x = k.a * 0x9E3779B185EBCA87ULL ^ k.b;
            x ^= (x >> 33);
            x *= 0xff51afd7ed558ccdULL;
            x ^= (x >> 33);
            x *= 0xc4ceb9fe1a85ec53ULL;
            x ^= (x >> 33);
            return x;
        }
        void rehash(size_t newCap) {
            std::vector<PairKey> oldK = std::move(_keys);
            std::vector<int> oldV = std::move(_vals);
            std::vector<uint8_t> oldS = std::move(_state);
            size_t oldCap = _capacity;
            _capacity = newCap; _size = 0;
            _keys.assign(_capacity, PairKey{});
            _vals.assign(_capacity, 0);
            _state.assign(_capacity, 0);
            if (oldCap == 0) return;
            for (size_t i = 0; i < oldCap; ++i) if (oldS[i] == 1) {
                PairKey k = oldK[i]; int v = oldV[i];
                int idx = findOrInsert(k);
                _vals[(size_t)idx] = v;
            }
        }
        int findOrInsert(const PairKey& k) {
            if ((_size + 1) * 2 >= _capacity) {
                size_t newCap = _capacity ? _capacity * 2 : 64;
                rehash(newCap);
            }
            size_t cap = _capacity;
            uint64_t h = hashOf(k);
            size_t i = (size_t)h & (cap - 1);
            for (;;) {
                if (_state[i] != 1) {
                    _keys[i] = k; _vals[i] = 0; _state[i] = 1; ++_size; return (int)i;
                }
                if (_keys[i].a == k.a && _keys[i].b == k.b) return (int)i;
                i = (i + 1) & (cap - 1);
            }
        }
    };
    FlatAdjMap adjacencyCounts; // consecutive frames adjacent

    // Constants
    float gravityG;               // gravitational constant
    float partRadius;             // particle radius (uniform)
    float bounceRestitution;      // boundary coefficient of restitution for elastic bounces
    static constexpr size_t MaxParticles = 300; // hard cap per requirements
    // Recycling pools: accumulated losses due to constraints (mass) and dissipation (energy)
    double massPool{0.0};
    double energyPool{0.0};
    static constexpr double MassSpawnUnit = 100.0;   // mass units to spawn a Z
    static constexpr double EnergySpawnUnit = 200.0; // energy units to spawn a Z

    // Runtime clock (ms spent while running)
    long long runAccumMs{0};

    // Cached legend text
    int legendCacheCols{-1};
    std::string legendCache;

    // Reusable buffers to reduce allocations
    struct Contact { size_t i; size_t j; float nx; float ny; float overlap; float rel; float vn1; float vn2; };
    std::vector<Contact> contacts;
    std::vector<std::vector<size_t>> broadGrid;
    int gridWCells{0};
    int gridHCells{0};
    float gridCellSize{0};
    
    // Barnes–Hut quadtree nodes
    struct QuadNode {
        float minx, miny, maxx, maxy;
        double mass; double cx, cy;
        int child[4];
        int particleIndex; // -1 if internal node
    };
    std::vector<QuadNode> quadTree;
    float bhTheta{0.6f};
    std::mutex qtMtx;

    // LBVH nodes (binary tree)
    struct BVHNode {
        float minx, miny, maxx, maxy;
        double mass; double cx, cy;
        int left;  // child index or -1
        int right; // child index or -1
        int leaf;  // particle index or -1
    };
    std::vector<BVHNode> bvhNodes;
    std::vector<std::pair<size_t,size_t>> fragSchedule;

    // SoA storage (groundwork for refactor)
    std::vector<unsigned long long> soa_id;
    std::vector<float> soa_x, soa_y, soa_vx, soa_vy;
    std::vector<uint8_t> soa_mass, soa_elast10, soa_decay, soa_alive;
    std::vector<char> soa_sym;
    std::vector<int16_t> soa_pix, soa_piy;
    void soaSyncFromVector();
    void soaSyncToVector();

    // Parallel helpers and CSR broad-phase
    int numWorkers{1};
    void initWorkers();
    void parallelFor(size_t n, const std::function<void(size_t,size_t,int)>& fn);
    // CSR grid: counts, offsets (size cells+1), indices of particle ids, and per-particle cell id
    std::vector<int> csrCounts;
    std::vector<int> csrOffsets;
    std::vector<size_t> csrIndex;
    std::vector<int> csrPcell;

    // Background stepping
    std::thread physThread;
    std::atomic<bool> threadExit{false};
    std::condition_variable cv;
    mutable std::mutex stateMtx; // for cv and start/stop signaling
    // Triple-buffer snapshots
    mutable std::mutex snapMtxDummy; // placeholder; readers don't lock; writer exclusive
    RenderSnapshot snaps[3];
    int snapWriteIdx{0};
    std::atomic<int> snapPublishedIdx{-1};
    std::atomic<uint64_t> snapSeq{0};
    // One-tick highlight grid (size w*h): 0 none, 1 split, 2 collision
    std::vector<uint8_t> highlightGrid;
    inline void markSplitAtCell(int ix, int iy) {
        if (ix<0||ix>=w||iy<0||iy>=h) return; size_t idx=(size_t)iy*(size_t)w+(size_t)ix; if (idx<highlightGrid.size()) highlightGrid[idx]=1;
    }
    inline void markCollisionAtCell(int ix, int iy) {
        if (ix<0||ix>=w||iy<0||iy>=h) return; size_t idx=(size_t)iy*(size_t)w+(size_t)ix; if (idx<highlightGrid.size()) highlightGrid[idx]=2;
    }
    inline void markSplit3x3AtCell(int ix, int iy) {
        for (int dy=-1; dy<=1; ++dy) for (int dx=-1; dx<=1; ++dx) markSplitAtCell(ix+dx, iy+dy);
    }
};
