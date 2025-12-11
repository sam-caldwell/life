/**
 * @file PhysicsWorld.cpp
 * @brief Newtonian particle simulation rendered with ncurses (physics variant).
 */
#include "PhysicsWorld.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <queue>
#include <unordered_set>

namespace {
inline int clampi(int v, int lo, int hi) { return std::max(lo, std::min(hi, v)); }
}

PhysicsWorld::PhysicsWorld(int width, int height,
                           float gravity,
                           float radius,
                           float restitution)
    : w(width), h(height), prng(rd()), gravityG(gravity), partRadius(radius), bounceRestitution(restitution) {}

PhysicsWorld::~PhysicsWorld() {}

void PhysicsWorld::setWindow(WINDOW* w_) {
    std::lock_guard<std::mutex> lock(mtx);
    win = w_;
}

void PhysicsWorld::setRunning(bool on) {
    running.store(on);
}

void PhysicsWorld::setStepDelayMs(int ms) {
    if (ms < 5) ms = 5;
    if (ms > 200) ms = 200; // tighter upper bound for smoother physics
    stepDelayMs.store(ms);
}

void PhysicsWorld::setGravity(float g) { gravityG = g; }
void PhysicsWorld::setRadius(float r) { partRadius = r; }
void PhysicsWorld::setRestitution(float e) { bounceRestitution = e; }

// No background thread; stepping is driven by the caller (UI loop)

void PhysicsWorld::clear() {
    std::lock_guard<std::mutex> lock(mtx);
    for (auto& p : particles) {
        if (win && p.prev_ix >= 0 && p.prev_iy >= 0) {
            mvwaddch(win, p.prev_iy, p.prev_ix, ' ');
        }
    }
    particles.clear();
    if (win) wrefresh(win);
}

void PhysicsWorld::reseed(unsigned count) {
    if (count > MaxParticles) count = MaxParticles;
    std::uniform_real_distribution<float> vxDist(-3.5f, 3.5f);
    std::uniform_real_distribution<float> vyDist(-2.0f, 2.0f);
    std::uniform_int_distribution<int> massDist(1, 100);
    std::uniform_int_distribution<int> symDist(0, 25);
    std::uniform_int_distribution<int> elastDist(0, 10);
    std::uniform_real_distribution<float> xDist(0.5f, std::max(0.5f, (float)w - 1.5f));
    std::uniform_real_distribution<float> yDist(0.5f, std::max(0.5f, (float)h - 1.5f));

    {
        std::lock_guard<std::mutex> lock(mtx);
        // clear existing
        for (auto& p : particles) {
            if (win && p.prev_ix >= 0 && p.prev_iy >= 0) mvwaddch(win, p.prev_iy, p.prev_ix, ' ');
        }
        particles.clear();
        particles.reserve(count);
        for (unsigned i = 0; i < count; ++i) {
            Particle p;
            p.id = nextId++;
            p.x = xDist(prng);
            p.y = yDist(prng);
            p.vx = vxDist(prng);
            p.vy = vyDist(prng);
            p.mass = massDist(prng);
            p.symbol = static_cast<char>('A' + symDist(prng));
            p.prev_ix = -1; p.prev_iy = -1;
            p.alive = true;
            p.elasticity10 = elastDist(prng);
            particles.push_back(p);
        }
        if (win) wrefresh(win);
    }
}

void PhysicsWorld::reseedRandom() {
    // Random count 5%..15% of drawable cells, capped at MaxParticles
    int maxCount = std::max(10, (w * h) / 7);
    int minCount = std::max(5, (w * h) / 30);
    if (maxCount > (int)MaxParticles) maxCount = (int)MaxParticles;
    if (minCount > (int)MaxParticles) minCount = (int)MaxParticles;
    if (minCount > maxCount) minCount = maxCount;
    std::uniform_int_distribution<int> cnt(minCount, maxCount);
    reseed((unsigned)cnt(prng));
}

void PhysicsWorld::drawAll(WINDOW* w_) {
    std::lock_guard<std::mutex> lock(mtx);
    WINDOW* ww = w_ ? w_ : win;
    if (!ww) return;
    werase(ww);
    for (auto& p : particles) {
        if (!p.alive) continue;
        int ix = clampi((int)std::round(p.x), 0, w - 1);
        int iy = clampi((int)std::round(p.y), 0, h - 1);
        int pair = colorPairForMass(p.mass);
        bool bold = (pair >= 9);
        if (bold) wattron(ww, A_BOLD);
        wattron(ww, COLOR_PAIR(pair));
        mvwaddch(ww, iy, ix, p.symbol);
        wattroff(ww, COLOR_PAIR(pair));
        if (bold) wattroff(ww, A_BOLD);
    }
    wrefresh(ww);
}

size_t PhysicsWorld::particleCount() const {
    std::lock_guard<std::mutex> lock(mtx);
    return particles.size();
}

int PhysicsWorld::colorPairForMass(int m) const {
    // Map [0,100] -> 16 buckets: 0..15 then shift to pairs 2..16 to avoid black (pair 1)
    if (m < 0) m = 0; if (m > 100) m = 100;
    int bucket = (int)std::floor((m / 100.0) * 16.0);
    if (bucket < 0) bucket = 0; if (bucket > 15) bucket = 15;
    int pair = bucket + 1; // 1..16
    if (pair < 2) pair = 2; // avoid pair 1 (black)
    if (pair > 16) pair = 16;
    return pair;
}

void PhysicsWorld::drawStatusLine(WINDOW* w_) {
    std::lock_guard<std::mutex> lock(mtx);
    WINDOW* ww = w_ ? w_ : win;
    if (!ww) return;
    int rows, cols; getmaxyx(ww, rows, cols);
    int y = rows - 1;
    wmove(ww, y, 0); wclrtoeol(ww);
    // Legend for mass buckets
    int x = 0;
    mvwprintw(ww, y, x, "M:"); x += 2;
    for (int i = 0; i < 16; ++i) {
        int m = (int)std::round((i / 15.0) * 100.0);
        int pair = std::max(2, i + 1); // skip black in legend
        bool bold = (i >= 8);
        if (bold) wattron(ww, A_BOLD);
        wattron(ww, COLOR_PAIR(pair));
        char buf[8]; snprintf(buf, sizeof(buf), "%d", m);
        mvwprintw(ww, y, x, "%s", buf);
        x += (int)strlen(buf);
        wattroff(ww, COLOR_PAIR(pair));
        if (bold) wattroff(ww, A_BOLD);
        if (i != 15) { mvwaddch(ww, y, x, ' '); ++x; }
    }
    // Build runtime clock HH:MM:SS.nn from accumulated milliseconds
    long long ms = runAccumMs;
    if (ms < 0) ms = 0;
    int hh = (int)(ms / (1000LL * 60LL * 60LL));
    ms %= (1000LL * 60LL * 60LL);
    int mm = (int)(ms / (1000LL * 60LL));
    ms %= (1000LL * 60LL);
    int ss = (int)(ms / 1000LL);
    int nn = (int)((ms % 1000LL) / 10LL);

    char status[256];
    snprintf(status, sizeof(status),
             "  | [s]tart/[p]ause  [r]eseed  [c]lear  speed[-/+]  [q]uit  | Particles: %zu/%zu  Delay(ms): %d  %s  | Time: %02d:%02d:%02d.%02d",
             particles.size(), MaxParticles, getStepDelayMs(), (running.load() ? "RUNNING" : "PAUSED"), hh, mm, ss, nn);
    if (x + (int)strlen(status) < cols) {
        mvwprintw(ww, y, x, "%s%*s", status, cols - x - (int)strlen(status), "");
    } else {
        mvwprintw(ww, y, x, "%s", status);
    }
    wrefresh(ww);
}

// No run loop; step(dt) is invoked by main

void PhysicsWorld::step(float dt) {
    // Advance running clock in milliseconds
    if (dt > 0) {
        long long inc = (long long)std::llround(dt * 1000.0f);
        if (inc > 0) runAccumMs += inc;
    }
    // Compute clusters for gravity aggregation
    std::vector<Cluster> clusters;
    computeClusters(clusters);

    // Apply gravity
    applyGravity(clusters, dt);

    // Substep integration to prevent tunneling through other particles.
    float maxDisp = 0.0f;
    {
        std::lock_guard<std::mutex> lock(mtx);
        for (auto& p : particles) {
            if (!p.alive) continue;
            float d = std::sqrt(p.vx*p.vx + p.vy*p.vy) * dt;
            if (d > maxDisp) maxDisp = d;
        }
    }
    const float maxStepLen = 0.5f; // max displacement per substep in cell units
    int substeps = (int)std::ceil(maxDisp / std::max(1e-6f, maxStepLen));
    if (substeps < 1) substeps = 1;
    float dts = dt / (float)substeps;

    for (int s = 0; s < substeps; ++s) {
        // Integrate one substep per particle (processor time per particle)
        {
            std::lock_guard<std::mutex> lock(mtx);
            for (auto& p : particles) {
                if (!p.alive) continue;
                updateParticle(p, dts);
            }
        }
        // Boundary and collisions per substep
        handleBoundary(partRadius);
        handleCollisions();
    }

    // Decay (Z -> 2xY, T..Y -> next lower) on 10-cycle intervals
    updateDecay();

    // After decay, apply adjacency-based combining
    updateAdjacencyAndCombine();

    // Incremental draw
    {
        std::lock_guard<std::mutex> lock(mtx);
        if (!win) return;
        for (auto& p : particles) {
            int ix = clampi((int)std::round(p.x), 0, w - 1);
            int iy = clampi((int)std::round(p.y), 0, h - 1);
            if (ix != p.prev_ix || iy != p.prev_iy) {
                if (p.prev_ix >= 0 && p.prev_iy >= 0) eraseParticleUnlocked(p);
                p.prev_ix = ix; p.prev_iy = iy;
                drawParticleUnlocked(p);
            } else {
                // Update color/attr in place in case mass mapping changed (mass constant here)
                drawParticleUnlocked(p);
            }
        }
        wrefresh(win);
    }
}

void PhysicsWorld::computeClusters(std::vector<Cluster>& out) {
    std::lock_guard<std::mutex> lock(mtx);
    out.clear();
    const int n = (int)particles.size();
    if (n == 0) return;
    // Map integer grid position -> list of particle indices
    std::unordered_map<long long, std::vector<int>> cellMap;
    cellMap.reserve((size_t)n * 2);
    auto keyOf = [](int x, int y) -> long long {
        return (static_cast<long long>(x) << 32) ^ static_cast<unsigned long long>(y);
    };
    std::vector<std::pair<int,int>> cellOf(n);
    for (int i = 0; i < n; ++i) {
        int ix = clampi((int)std::round(particles[i].x), 0, w - 1);
        int iy = clampi((int)std::round(particles[i].y), 0, h - 1);
        cellOf[i] = {ix, iy};
        cellMap[keyOf(ix, iy)].push_back(i);
    }
    // Build adjacency graph where nodes are particle indices; edges if their cells are 8-neighbors
    std::vector<std::vector<int>> adj(n);
    for (int i = 0; i < n; ++i) {
        auto [xi, yi] = cellOf[i];
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0) continue;
                int nx = xi + dx, ny = yi + dy;
                auto it = cellMap.find(keyOf(nx, ny));
                if (it == cellMap.end()) continue;
                for (int j : it->second) adj[i].push_back(j);
            }
        }
    }
    // Find connected components
    std::vector<int> comp(n, -1);
    int cid = 0;
    for (int i = 0; i < n; ++i) if (comp[i] < 0) {
        std::queue<int> q; q.push(i); comp[i] = cid;
        Cluster c{};
        while (!q.empty()) {
            int u = q.front(); q.pop();
            c.members.push_back((size_t)u);
            double m = std::max(0, particles[u].mass);
            c.mass += m;
            c.cx += m * particles[u].x;
            c.cy += m * particles[u].y;
            for (int v : adj[u]) if (comp[v] < 0) { comp[v] = cid; q.push(v); }
        }
        if (c.mass > 0) { c.cx /= c.mass; c.cy /= c.mass; }
        out.push_back(c);
        ++cid;
    }
}

void PhysicsWorld::applyGravity(const std::vector<Cluster>& clusters, float dt) {
    std::lock_guard<std::mutex> lock(mtx);
    if (particles.empty()) return;
    // For each particle, sum acceleration from all clusters except the one containing it (approximation)
    // Build quick lookup: particle index -> cluster idx
    std::unordered_map<size_t,int> p2c;
    p2c.reserve(particles.size());
    for (size_t ci = 0; ci < clusters.size(); ++ci) {
        for (size_t pi : clusters[ci].members) p2c[pi] = (int)ci;
    }
    for (size_t i = 0; i < particles.size(); ++i) {
        auto& p = particles[i];
        if (!p.alive) continue;
        float ax = 0.0f, ay = 0.0f;
        for (size_t ci = 0; ci < clusters.size(); ++ci) {
            auto& c = clusters[ci];
            auto it = p2c.find(i);
            if (it != p2c.end() && it->second == (int)ci) {
                // Apply influence of the rest of the cluster excluding self
                double massExcl = c.mass - std::max(0, p.mass);
                if (massExcl <= 0.0) continue;
                double cxExcl = (c.mass * c.cx - p.mass * p.x) / massExcl;
                double cyExcl = (c.mass * c.cy - p.mass * p.y) / massExcl;
                double dx = cxExcl - p.x;
                double dy = cyExcl - p.y;
                double r2 = dx*dx + dy*dy + 0.0001; // softening epsilon
                double invr = 1.0 / std::sqrt(r2);
                double f = (gravityG * massExcl) / r2;
                ax += (float)(f * dx * invr);
                ay += (float)(f * dy * invr);
            } else {
                double dx = c.cx - p.x;
                double dy = c.cy - p.y;
                double r2 = dx*dx + dy*dy + 0.0001;
                double invr = 1.0 / std::sqrt(r2);
                double f = (gravityG * c.mass) / r2;
                ax += (float)(f * dx * invr);
                ay += (float)(f * dy * invr);
            }
        }
        p.vx += ax * dt;
        p.vy += ay * dt;
    }
}

void PhysicsWorld::handleBoundary(float padding) {
    std::lock_guard<std::mutex> lock(mtx);
    float minx = 0.0f + padding, miny = 0.0f + padding;
    float maxx = (float)(w - 1) - padding;
    float maxy = (float)(h - 1) - padding;
    for (auto& p : particles) {
        if (!p.alive) continue;
        if (p.x < minx) { p.x = minx; p.vx = -p.vx * bounceRestitution; }
        if (p.x > maxx) { p.x = maxx; p.vx = -p.vx * bounceRestitution; }
        if (p.y < miny) { p.y = miny; p.vy = -p.vy * bounceRestitution; }
        if (p.y > maxy) { p.y = maxy; p.vy = -p.vy * bounceRestitution; }
    }
}

void PhysicsWorld::handleCollisions() {
    std::lock_guard<std::mutex> lock(mtx);
    const float R = partRadius;
    const float R2 = (2*R)*(2*R);
    const size_t n = particles.size();

    struct Pair { size_t i; size_t j; float nx; float ny; float overlap; float rel; float vn1; float vn2; };
    std::vector<Pair> contacts;
    contacts.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            auto& a = particles[i];
            auto& b = particles[j];
            if (!a.alive || !b.alive) continue;
            float dx = a.x - b.x;
            float dy = a.y - b.y;
            float dist2 = dx*dx + dy*dy;
            if (dist2 <= R2) {
                float dist = std::sqrt(std::max(dist2, 1e-6f));
                float nx = dx / dist;
                float ny = dy / dist;
                float overlap = (2*R) - dist;
                float vn1 = a.vx * nx + a.vy * ny;
                float vn2 = b.vx * nx + b.vy * ny;
                float rel = (a.vx - b.vx) * nx + (a.vy - b.vy) * ny;
                contacts.push_back({i,j,nx,ny,overlap,rel,vn1,vn2});
            }
        }
    }

    // Schedule fragmentations for high-momentum impacts
    // Threshold defined as 10x maximum mass of a single particle (100)
    const float MOMENTUM_THRESHOLD = 10.0f * 100.0f; // 10x max mass -> 1000
    std::vector<std::pair<size_t,size_t>> toFragment;
    std::unordered_set<size_t> deadIdx;
    for (auto& c : contacts) {
        auto& a = particles[c.i];
        auto& b = particles[c.j];
        float p_rel = std::fabs((float)a.mass * c.vn1 - (float)b.mass * c.vn2);
        if (p_rel > MOMENTUM_THRESHOLD) {
            toFragment.emplace_back(c.i, c.j);
        }
    }

    // Apply fragmentations
    auto spawnChild = [&](const Particle& parent, float offx, float offy) {
        Particle child;
        child.id = nextId++;
        child.symbol = (parent.symbol > 'A') ? (char)(parent.symbol - 1) : 'A';
        int m1 = parent.mass / 2;
        // Alternate to distribute remainder across two calls: caller will call twice
        child.mass = m1;
        child.elasticity10 = parent.elasticity10;
        child.x = parent.x + offx;
        child.y = parent.y + offy;
        child.vx = parent.vx;
        child.vy = parent.vy;
        child.prev_ix = -1; child.prev_iy = -1;
        child.alive = true;
        return child;
    };

    std::vector<Particle> newParticles;
    newParticles.reserve(toFragment.size() * 4);
    // Track live count to enforce MaxParticles capacity during fragmentation
    size_t liveCount = 0; for (auto& p : particles) if (p.alive) ++liveCount;
    for (auto [i,j] : toFragment) {
        if (deadIdx.count(i) || deadIdx.count(j)) continue;
        auto& a = particles[i];
        auto& b = particles[j];
        if (!a.alive || !b.alive) continue;

        // Erase originals from screen
        if (win) {
            if (a.prev_ix >= 0 && a.prev_iy >= 0) mvwaddch(win, a.prev_iy, a.prev_ix, ' ');
            if (b.prev_ix >= 0 && b.prev_iy >= 0) mvwaddch(win, b.prev_iy, b.prev_ix, ' ');
        }
        a.alive = false; b.alive = false;
        deadIdx.insert(i); deadIdx.insert(j);

        // Create 2 children from each parent, slightly offset orthogonal to their mutual normal
        // Find contact normal from contacts list
        float nx = 0, ny = 0;
        for (auto& c : contacts) if ((c.i==i && c.j==j) || (c.i==j && c.j==i)) { nx = c.nx; ny = c.ny; break; }
        // Tangent vector
        float tx = -ny, ty = nx;
        float sep = std::max(0.1f, partRadius * 0.6f);

        // For a
        Particle a1 = spawnChild(a, tx * sep, ty * sep);
        Particle a2 = spawnChild(a, -tx * sep, -ty * sep);
        // Mass split with minimum 1 per child when possible
        a1.mass = a.mass / 2;
        a2.mass = a.mass - a1.mass;
        if (a.mass >= 2) {
            if (a1.mass < 1) { a1.mass = 1; a2.mass = a.mass - 1; }
            if (a2.mass < 1) { a2.mass = 1; a1.mass = a.mass - 1; }
        }
        // For b
        Particle b1 = spawnChild(b, tx * sep, ty * sep);
        Particle b2 = spawnChild(b, -tx * sep, -ty * sep);
        b1.mass = b.mass / 2;
        b2.mass = b.mass - b1.mass;
        if (b.mass >= 2) {
            if (b1.mass < 1) { b1.mass = 1; b2.mass = b.mass - 1; }
            if (b2.mass < 1) { b2.mass = 1; b1.mass = b.mass - 1; }
        }

        // Enforce capacity: remove parents (net -2), then add children up to MaxParticles
        liveCount -= 2;
        auto tryAdd = [&](const Particle& c){ if (c.mass >= 1 && liveCount < MaxParticles) { newParticles.push_back(c); ++liveCount; } };
        // Only add non-zero mass children
        tryAdd(a1);
        tryAdd(a2);
        tryAdd(b1);
        tryAdd(b2);
    }

    // Resolve remaining contacts with collision response (skip dead pairs)
    for (auto& c : contacts) {
        if (deadIdx.count(c.i) || deadIdx.count(c.j)) continue;
        auto& a = particles[c.i];
        auto& b = particles[c.j];
        if (!a.alive || !b.alive) continue;

        // Position correction to resolve overlap
        if (c.overlap > 0) {
            float totalMass = (float)a.mass + (float)b.mass + 1e-4f;
            float pushA = (float)b.mass / totalMass;
            float pushB = (float)a.mass / totalMass;
            a.x += c.nx * c.overlap * pushA;
            a.y += c.ny * c.overlap * pushA;
            b.x -= c.nx * c.overlap * pushB;
            b.y -= c.ny * c.overlap * pushB;
        }

        // Normal impulse
        float rel = c.rel;
        if (rel > 0) continue; // separating
        float m1 = std::max(1, a.mass);
        float m2 = std::max(1, b.mass);
        float ea = std::max(0, a.elasticity10) / 10.0f;
        float eb = std::max(0, b.elasticity10) / 10.0f;
        float e = std::max(0.0f, std::min(1.0f, std::min(ea, eb)));
        float jimp = -(1 + e) * rel / (1.0f/m1 + 1.0f/m2);
        float jx = jimp * c.nx;
        float jy = jimp * c.ny;
        a.vx += jx / m1; a.vy += jy / m1;
        b.vx -= jx / m2; b.vy -= jy / m2;
    }

    // Remove dead and append new particles
    if (!deadIdx.empty() || !newParticles.empty()) {
        std::vector<Particle> out;
        out.reserve(particles.size() - deadIdx.size() + newParticles.size());
        for (size_t i = 0; i < particles.size(); ++i) {
            if (deadIdx.count(i)) continue;
            out.push_back(particles[i]);
        }
        for (auto& np : newParticles) out.push_back(np);
        particles.swap(out);
    }
}

void PhysicsWorld::updateDecay() {
    std::lock_guard<std::mutex> lock(mtx);
    if (particles.empty()) return;

    std::vector<size_t> toRemove;
    std::vector<Particle> toAdd;

    auto makeChild = [&](const Particle& parent, char sym, int mass, float offx, float offy) {
        Particle c;
        c.id = nextId++;
        c.symbol = sym;
        c.mass = std::max(1, std::min(100, mass));
        c.elasticity10 = parent.elasticity10;
        c.x = parent.x + offx;
        c.y = parent.y + offy;
        c.vx = parent.vx;
        c.vy = parent.vy;
        c.prev_ix = -1; c.prev_iy = -1;
        c.alive = true;
        c.decayTicks = 0;
        return c;
    };

    // Count current live size
    size_t liveCount = 0; for (auto& p : particles) if (p.alive) ++liveCount;

    for (size_t i = 0; i < particles.size(); ++i) {
        auto& p = particles[i];
        if (!p.alive) continue;
        p.decayTicks += 1;
        if (p.decayTicks < 10) continue;
        // Time to decay
        p.decayTicks = 0;
        if (p.symbol == 'Z') {
            // Z -> 2xY (if capacity allows)
            if (liveCount >= MaxParticles) continue; // cannot add more
            if (p.mass < 2) { p.symbol = 'Y'; continue; }
            char sym = 'Y';
            int m1 = p.mass / 2;
            int m2 = p.mass - m1;
            // Choose a random outward direction and place children on opposite sides; give them unit speed away.
            std::uniform_real_distribution<float> ang(0.0f, (float)M_PI * 2.0f);
            float theta = ang(prng);
            float ux = std::cos(theta);
            float uy = std::sin(theta);
            float sep = std::max(0.1f, partRadius * 0.6f);

            size_t addCap = (MaxParticles > liveCount) ? (MaxParticles - liveCount) : 0;
            // We aim to replace 1 particle with 2 -> net +1
            if (addCap == 0) continue;
            // remove parent
            if (win && p.prev_ix >= 0 && p.prev_iy >= 0) mvwaddch(win, p.prev_iy, p.prev_ix, ' ');
            p.alive = false; toRemove.push_back(i); --liveCount;
            // create children within cap
            Particle c1 = makeChild(p, sym, m1, ux*sep, uy*sep);
            c1.vx = ux * 1.0f; c1.vy = uy * 1.0f;
            toAdd.push_back(c1); ++liveCount;
            if (liveCount < MaxParticles) {
                Particle c2 = makeChild(p, sym, m2, -ux*sep, -uy*sep);
                c2.vx = -ux * 1.0f; c2.vy = -uy * 1.0f;
                toAdd.push_back(c2); ++liveCount;
            }
        }
    }

    if (!toRemove.empty() || !toAdd.empty()) {
        std::vector<Particle> out;
        out.reserve(particles.size() - toRemove.size() + toAdd.size());
        std::unordered_set<size_t> rem(toRemove.begin(), toRemove.end());
        for (size_t i = 0; i < particles.size(); ++i) if (!rem.count(i)) out.push_back(particles[i]);
        for (auto& c : toAdd) out.push_back(c);
        particles.swap(out);
    }
}

void PhysicsWorld::updateAdjacencyAndCombine() {
    std::lock_guard<std::mutex> lock(mtx);
    if (particles.empty()) return;
    // Build cell map id -> (ix,iy)
    std::unordered_map<long long, std::vector<size_t>> cellMap;
    cellMap.reserve(particles.size()*2);
    auto keyOf = [](int x, int y) -> long long { return (static_cast<long long>(x) << 32) ^ (unsigned long long)y; };
    std::vector<std::pair<int,int>> cellOf(particles.size());
    for (size_t i = 0; i < particles.size(); ++i) {
        int ix = clampi((int)std::round(particles[i].x), 0, w - 1);
        int iy = clampi((int)std::round(particles[i].y), 0, h - 1);
        cellOf[i] = {ix, iy};
        cellMap[keyOf(ix, iy)].push_back(i);
    }

    // Compute set of currently adjacent pairs (8-neighborhood)
    std::unordered_set<unsigned long long> present; // hashed PairKey via 128->64 fold
    auto pkhash = [&](unsigned long long a, unsigned long long b){ unsigned long long x = (a<<1) ^ (b + 0x9e3779b97f4a7c15ULL + (a<<6) + (a>>2)); return x; };
    auto mkpair = [&](unsigned long long a, unsigned long long b){ return a<b ? std::make_pair(a,b) : std::make_pair(b,a); };
    std::vector<PairKey> presentKeys;

    for (size_t i = 0; i < particles.size(); ++i) {
        auto [xi, yi] = cellOf[i];
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0) continue;
                int nx = xi + dx, ny = yi + dy;
                auto it = cellMap.find(keyOf(nx, ny));
                if (it == cellMap.end()) continue;
                for (size_t j : it->second) {
                    if (j <= i) continue;
                    auto pr = mkpair(particles[i].id, particles[j].id);
                    unsigned long long h = pkhash(pr.first, pr.second);
                    if (present.insert(h).second) {
                        presentKeys.push_back({pr.first, pr.second});
                    }
                }
            }
        }
    }

    // Increment counts; remove stale entries
    std::unordered_set<unsigned long long> keep;
    for (auto& k : presentKeys) {
        PairKey key{k.a, k.b};
        adjacencyCounts[key] += 1;
        unsigned long long hh = pkhash(k.a, k.b);
        keep.insert(hh);
    }
    // Collect stale
    std::vector<PairKey> toErase;
    toErase.reserve(adjacencyCounts.size());
    for (auto& kv : adjacencyCounts) {
        unsigned long long hh = pkhash(kv.first.a, kv.first.b);
        if (!keep.count(hh)) toErase.push_back(kv.first);
    }
    for (auto& k : toErase) adjacencyCounts.erase(k);

    // Attempt combinations where count >= 10
    std::unordered_set<unsigned long long> consumed; // particle ids that already combined
    std::vector<size_t> removeIdx;
    std::vector<Particle> addList;
    std::uniform_int_distribution<int> coin(0,1);
    for (auto it = adjacencyCounts.begin(); it != adjacencyCounts.end(); ) {
        if (it->second >= 10) {
            // Find indices for ids
            unsigned long long ida = it->first.a, idb = it->first.b;
            if (consumed.count(ida) || consumed.count(idb)) { it = adjacencyCounts.erase(it); continue; }
            // 50/50 chance
            if (coin(prng) == 1) {
                // locate particles
                size_t ia = SIZE_MAX, ib = SIZE_MAX;
                for (size_t i = 0; i < particles.size(); ++i) {
                    if (particles[i].id == ida) ia = i; else if (particles[i].id == idb) ib = i;
                }
                if (ia != SIZE_MAX && ib != SIZE_MAX && particles[ia].alive && particles[ib].alive) {
                    auto& A = particles[ia]; auto& B = particles[ib];
                    int iaL = (A.symbol - 'A' + 1);
                    int ibL = (B.symbol - 'A' + 1);
                    int sumL = iaL + ibL;
                    if (sumL > 26) sumL = 26;
                    int massSum = A.mass + B.mass;
                    if (massSum > 100) massSum = 100;
                    char sym = (char)('A' + (sumL - 1));
                    // Momentum-conserving velocity
                    float m1 = std::max(1, A.mass);
                    float m2 = std::max(1, B.mass);
                    float vx = (m1*A.vx + m2*B.vx) / (m1 + m2);
                    float vy = (m1*A.vy + m2*B.vy) / (m1 + m2);
                    // New particle at first's position (choose smaller id as first)
                    size_t firstIdx = (A.id < B.id) ? ia : ib;
                    auto& F = particles[firstIdx];
                    Particle C;
                    C.id = nextId++;
                    C.symbol = sym;
                    C.mass = massSum;
                    C.elasticity10 = (A.elasticity10 + B.elasticity10) / 2;
                    C.x = F.x; C.y = F.y;
                    C.vx = vx; C.vy = vy;
                    C.prev_ix = -1; C.prev_iy = -1;
                    C.alive = true;
                    addList.push_back(C);
                    // Mark originals for removal and clear from screen
                    if (win) {
                        if (A.prev_ix >= 0 && A.prev_iy >= 0) mvwaddch(win, A.prev_iy, A.prev_ix, ' ');
                        if (B.prev_ix >= 0 && B.prev_iy >= 0) mvwaddch(win, B.prev_iy, B.prev_ix, ' ');
                    }
                    particles[ia].alive = false; particles[ib].alive = false;
                    removeIdx.push_back(ia); removeIdx.push_back(ib);
                    consumed.insert(ida); consumed.insert(idb);
                    it = adjacencyCounts.erase(it);
                    continue;
                } else {
                    it = adjacencyCounts.erase(it);
                    continue;
                }
            } else {
                it->second = 0; // reset counter if chance failed
                ++it;
            }
        } else {
            ++it;
        }
    }

    if (!removeIdx.empty() || !addList.empty()) {
        // Compact particle list and append additions
        std::vector<Particle> out;
        out.reserve(particles.size() - removeIdx.size() + addList.size());
        std::unordered_set<size_t> remset(removeIdx.begin(), removeIdx.end());
        for (size_t i = 0; i < particles.size(); ++i) if (!remset.count(i)) out.push_back(particles[i]);
        for (auto& np : addList) out.push_back(np);
        particles.swap(out);
        // Clean adjacencyCounts entries involving consumed ids
        std::vector<PairKey> toErase2;
        for (auto& kv : adjacencyCounts) if (consumed.count(kv.first.a) || consumed.count(kv.first.b)) toErase2.push_back(kv.first);
        for (auto& k : toErase2) adjacencyCounts.erase(k);
    }
}

void PhysicsWorld::drawParticleUnlocked(const Particle& p) {
    if (!win) return;
    int ix = clampi((int)std::round(p.x), 0, w - 1);
    int iy = clampi((int)std::round(p.y), 0, h - 1);
    int pair = colorPairForMass(p.mass);
    bool bold = (pair >= 9);
    if (bold) wattron(win, A_BOLD);
    wattron(win, COLOR_PAIR(pair));
    mvwaddch(win, iy, ix, p.symbol);
    wattroff(win, COLOR_PAIR(pair));
    if (bold) wattroff(win, A_BOLD);
}

void PhysicsWorld::eraseParticleUnlocked(const Particle& p) {
    if (!win) return;
    int ix = clampi((int)std::round(p.prev_ix), 0, w - 1);
    int iy = clampi((int)std::round(p.prev_iy), 0, h - 1);
    mvwaddch(win, iy, ix, ' ');
}
void PhysicsWorld::updateParticle(Particle& p, float dts) {
    p.x += p.vx * dts;
    p.y += p.vy * dts;
}
