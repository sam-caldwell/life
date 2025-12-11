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
    : w(width), h(height), prng(rd()), gravityG(gravity), partRadius(radius), bounceRestitution(restitution) {
    initWorkers();
}

PhysicsWorld::~PhysicsWorld() {}
 
bool PhysicsWorld::getLatestSnapshot(const std::vector<RenderItem>*& items, uint64_t& seq, int& outW, int& outH) const {
    int idx = snapPublishedIdx.load(std::memory_order_acquire);
    if (idx < 0) return false;
    const RenderSnapshot& s = snaps[idx];
    items = &s.items;
    seq = s.seq.load(std::memory_order_relaxed);
    outW = s.w; outH = s.h;
    return true;
}

void PhysicsWorld::startThread() {
    if (physThread.joinable()) return;
    threadExit.store(false, std::memory_order_relaxed);
    physThread = std::thread([this]{
        using clock = std::chrono::steady_clock;
        auto last = clock::now();
        double accumulator = 0.0;
        while (!threadExit.load(std::memory_order_relaxed)) {
            // Pause while not running
            if (!running.load(std::memory_order_acquire)) {
                std::unique_lock<std::mutex> lk(stateMtx);
                cv.wait_for(lk, std::chrono::milliseconds(10));
                last = clock::now();
                continue;
            }
            auto now = clock::now();
            double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(now - last).count();
            last = now;
            double fixedDt = std::max(0.005, std::min(0.2, getStepDelayMs() / 1000.0));
            accumulator += elapsed;
            int maxSteps = 8;
            bool didStep = false;
            while (accumulator >= fixedDt && maxSteps-- > 0) {
                step((float)fixedDt); // headless: step does not draw
                accumulator -= fixedDt;
                didStep = true;
            }
            if (didStep) {
                // Publish snapshot
                RenderSnapshot& dst = snaps[snapWriteIdx];
                dst.w = w; dst.h = h;
                // Build items from SoA
                size_t n = soa_id.size();
                dst.items.clear();
                dst.items.reserve(n);
                for (size_t i=0;i<n;++i) {
                    if (!soa_alive[i]) continue;
                    int ix = clampi((int)std::round(soa_x[i]), 0, w - 1);
                    int iy = clampi((int)std::round(soa_y[i]), 0, h - 1);
                    RenderItem it{(int16_t)ix, (int16_t)iy, soa_sym[i], soa_mass[i], (uint8_t)1};
                    dst.items.push_back(it);
                }
                uint64_t s = snapSeq.fetch_add(1, std::memory_order_acq_rel) + 1;
                dst.seq.store(s, std::memory_order_release);
                snapPublishedIdx.store(snapWriteIdx, std::memory_order_release);
                // Advance write index to a buffer different from published
                int next = (snapWriteIdx + 1) % 3;
                if (next == snapPublishedIdx.load(std::memory_order_acquire)) next = (next + 1) % 3;
                snapWriteIdx = next;
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    });
}

void PhysicsWorld::stopThread() {
    if (!physThread.joinable()) return;
    threadExit.store(true, std::memory_order_relaxed);
    notifyThread();
    physThread.join();
}

void PhysicsWorld::notifyThread() {
    cv.notify_all();
}

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
    if (win) { werase(win); wrefresh(win); }
    soa_id.clear(); soa_x.clear(); soa_y.clear(); soa_vx.clear(); soa_vy.clear();
    soa_mass.clear(); soa_elast10.clear(); soa_decay.clear(); soa_alive.clear();
    soa_sym.clear(); soa_pix.clear(); soa_piy.clear();
}

void PhysicsWorld::reseed(unsigned count) {
    reseedSoA(count);
}

void PhysicsWorld::reseedSoA(unsigned count) {
    size_t cellCap = (size_t)std::max(0, w) * (size_t)std::max(0, h);
    if (count > MaxParticles) count = MaxParticles;
    if (count > cellCap) count = (unsigned)cellCap;
    // Clear SoA and vector glyphs
    {
        std::lock_guard<std::mutex> lock(mtx);
        if (win) werase(win);
    }
    soa_id.clear(); soa_x.clear(); soa_y.clear(); soa_vx.clear(); soa_vy.clear(); soa_mass.clear(); soa_elast10.clear(); soa_decay.clear(); soa_alive.clear(); soa_sym.clear(); soa_pix.clear(); soa_piy.clear();
    soa_id.reserve(count); soa_x.reserve(count); soa_y.reserve(count); soa_vx.reserve(count); soa_vy.reserve(count); soa_mass.reserve(count); soa_elast10.reserve(count); soa_decay.reserve(count); soa_alive.reserve(count); soa_sym.reserve(count); soa_pix.reserve(count); soa_piy.reserve(count);
    std::uniform_real_distribution<float> vxDist(-3.5f, 3.5f);
    std::uniform_real_distribution<float> vyDist(-2.0f, 2.0f);
    std::uniform_int_distribution<int> massDist(1, 100);
    std::uniform_int_distribution<int> symDist(0, 25);
    std::uniform_int_distribution<int> elastDist(0, 10);
    std::uniform_int_distribution<int> ixDist(0, std::max(0, w - 1));
    std::uniform_int_distribution<int> iyDist(0, std::max(0, h - 1));
    std::unordered_set<long long> occ; occ.reserve(count*2+16);
    unsigned placed=0; int attempts=0; const int maxAttempts = (int)count * 50 + 500;
    while (placed < count && attempts < maxAttempts) {
        ++attempts; int ix = ixDist(prng); int iy = iyDist(prng); long long key = ((long long)iy << 32) ^ (unsigned long long)ix; if (occ.count(key)) continue; occ.insert(key);
        soa_id.push_back(nextId++); soa_x.push_back((float)ix); soa_y.push_back((float)iy); soa_vx.push_back(vxDist(prng)); soa_vy.push_back(vyDist(prng)); soa_mass.push_back((uint8_t)massDist(prng)); soa_elast10.push_back((uint8_t)elastDist(prng)); soa_decay.push_back(0); soa_alive.push_back(1); soa_sym.push_back((char)('A'+symDist(prng))); soa_pix.push_back(-1); soa_piy.push_back(-1); ++placed;
    }
}

void PhysicsWorld::reseedRandom() {
    // Random count 5%..15% of drawable cells, capped at MaxParticles
    int maxCount = std::max(10, (w * h) / 7);
    int minCount = std::max(5, (w * h) / 30);
    int cellCapI = std::max(0, w) * std::max(0, h);
    if (maxCount > (int)MaxParticles) maxCount = (int)MaxParticles;
    if (minCount > (int)MaxParticles) minCount = (int)MaxParticles;
    if (maxCount > cellCapI) maxCount = cellCapI;
    if (minCount > cellCapI) minCount = cellCapI;
    if (minCount > maxCount) minCount = maxCount;
    std::uniform_int_distribution<int> cnt(minCount, maxCount);
    reseedSoA((unsigned)cnt(prng));
}

void PhysicsWorld::drawAll(WINDOW* w_) {
    std::lock_guard<std::mutex> lock(mtx);
    WINDOW* ww = w_ ? w_ : win;
    if (!ww) return;
    werase(ww);
    for (size_t i=0;i<soa_id.size();++i) {
        if (!soa_alive[i]) continue;
        int ix = clampi((int)std::round(soa_x[i]), 0, w - 1);
        int iy = clampi((int)std::round(soa_y[i]), 0, h - 1);
        int pair = colorPairForMass(soa_mass[i]);
        bool bold = (pair >= 9);
        if (bold) wattron(ww, A_BOLD);
        wattron(ww, COLOR_PAIR(pair));
        mvwaddch(ww, iy, ix, soa_sym[i]);
        wattroff(ww, COLOR_PAIR(pair));
        if (bold) wattroff(ww, A_BOLD);
    }
    wnoutrefresh(ww);
}

size_t PhysicsWorld::particleCount() const {
    std::lock_guard<std::mutex> lock(mtx);
    return soa_id.size();
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
    // Cache legend across frames for given width
    int x = 0;
    if (legendCacheCols != cols) {
        legendCache.clear();
        legendCache += "M:";
        for (int i = 0; i < 16; ++i) {
            int m = (int)std::round((i / 15.0) * 100.0);
            char buf[8]; snprintf(buf, sizeof(buf), "%d", m);
            legendCache += buf;
            if (i != 15) legendCache += ' ';
        }
        legendCacheCols = cols;
    }
    // Print monochrome cached text, then color over it for readability
    mvwprintw(ww, y, x, "%s", legendCache.c_str());
    x += (int)legendCache.size();
    // Color overlay for mass buckets
    int lx = 2; // after "M:"
    for (int i = 0; i < 16; ++i) {
        int m = (int)std::round((i / 15.0) * 100.0);
        int pair = std::max(2, i + 1); // skip black
        bool bold = (i >= 8);
        if (bold) wattron(ww, A_BOLD);
        wattron(ww, COLOR_PAIR(pair));
        char buf[8]; snprintf(buf, sizeof(buf), "%d", m);
        mvwprintw(ww, y, lx, "%s", buf);
        lx += (int)strlen(buf) + (i==15?0:1);
        wattroff(ww, COLOR_PAIR(pair));
        if (bold) wattroff(ww, A_BOLD);
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
             soa_id.size(), MaxParticles, getStepDelayMs(), (running.load() ? "RUNNING" : "PAUSED"), hh, mm, ss, nn);
    if (x + (int)strlen(status) < cols) {
        mvwprintw(ww, y, x, "%s%*s", status, cols - x - (int)strlen(status), "");
    } else {
        mvwprintw(ww, y, x, "%s", status);
    }
    wnoutrefresh(ww);
}

// No run loop; step(dt) is invoked by main

void PhysicsWorld::step(float dt) {
    // Advance running clock in milliseconds
    if (dt > 0) {
        long long inc = (long long)std::llround(dt * 1000.0f);
        if (inc > 0) runAccumMs += inc;
    }
    // Apply gravity via Barnes–Hut directly on SoA
    applyGravityBarnesHutSoA(dt);

    // Substep integration to prevent tunneling through other particles.
    float maxDisp = 0.0f;
    for (size_t i = 0; i < soa_x.size(); ++i) {
        if (!soa_alive[i]) continue;
        float d = std::sqrt(soa_vx[i]*soa_vx[i] + soa_vy[i]*soa_vy[i]) * dt;
        if (d > maxDisp) maxDisp = d;
    }
    const float maxStepLen = 0.5f; // max displacement per substep in cell units
    int substeps = (int)std::ceil(maxDisp / std::max(1e-6f, maxStepLen));
    if (substeps < 1) substeps = 1;
    float dts = dt / (float)substeps;

    fragSchedule.clear();
    for (int s = 0; s < substeps; ++s) {
        updateParticleSoA(dts);
        handleBoundarySoA(partRadius);
        handleCollisionsSoA();
    }

    // Apply scheduled fragmentations, decay, combining, and pool spawns in SoA
    applyFragmentationsSoAFromSchedule();
    updateDecaySoA();
    updateAdjacencyAndCombineSoA();
    maybeSpawnFromPoolsSoA();

    // Incremental draw
    {
        std::lock_guard<std::mutex> lock(mtx);
        if (!win) return;
        for (size_t i=0;i<soa_id.size();++i) {
            if (!soa_alive[i]) continue;
            int ix = clampi((int)std::round(soa_x[i]), 0, w - 1);
            int iy = clampi((int)std::round(soa_y[i]), 0, h - 1);
            if (ix != soa_pix[i] || iy != soa_piy[i]) {
                if (soa_pix[i] >= 0 && soa_piy[i] >= 0) {
                    int ex = clampi((int)soa_pix[i], 0, w - 1);
                    int ey = clampi((int)soa_piy[i], 0, h - 1);
                    mvwaddch(win, ey, ex, ' ');
                }
                soa_pix[i] = (int16_t)ix; soa_piy[i] = (int16_t)iy;
                int pair = colorPairForMass(soa_mass[i]); bool bold = (pair >= 9);
                if (bold) wattron(win, A_BOLD);
                wattron(win, COLOR_PAIR(pair));
                mvwaddch(win, iy, ix, soa_sym[i]);
                wattroff(win, COLOR_PAIR(pair)); if (bold) wattroff(win, A_BOLD);
            }
        }
        wnoutrefresh(win);
    }
}

void PhysicsWorld::applyGravityBarnesHutSoA(float dt) {
    // Build LBVH over alive particles (binary tree) using Morton ordering and median splits
    bvhNodes.clear(); bvhNodes.reserve(std::max<size_t>(64, soa_x.size() * 2));
    auto makeNode = [&](){ BVHNode n{}; n.minx=1e9f; n.miny=1e9f; n.maxx=-1e9f; n.maxy=-1e9f; n.mass=0; n.cx=0; n.cy=0; n.left=-1; n.right=-1; n.leaf=-1; bvhNodes.push_back(n); return (int)bvhNodes.size()-1; };
    // Build index list of alive particles and sort by Morton code for locality
    auto morton2d = [&](uint32_t x, uint32_t y){
        auto part1by1 = [](uint32_t n){ n = (n | (n << 8)) & 0x00FF00FFu; n = (n | (n << 4)) & 0x0F0F0F0Fu; n = (n | (n << 2)) & 0x33333333u; n = (n | (n << 1)) & 0x55555555u; return n; };
        return (part1by1(y) << 1) | part1by1(x);
    };
    std::vector<size_t> aliveIdx; aliveIdx.reserve(soa_id.size());
    for (size_t i=0;i<soa_id.size();++i) if (soa_alive[i]) aliveIdx.push_back(i);
    std::vector<std::pair<uint32_t,size_t>> mort; mort.reserve(aliveIdx.size());
    for (size_t i : aliveIdx) {
        int ix = clampi((int)std::round(soa_x[i]), 0, w - 1);
        int iy = clampi((int)std::round(soa_y[i]), 0, h - 1);
        mort.emplace_back(morton2d((uint32_t)ix, (uint32_t)iy), i);
    }
    std::sort(mort.begin(), mort.end(), [](auto& a, auto& b){ return a.first < b.first; });
    aliveIdx.clear(); aliveIdx.reserve(mort.size()); for (auto& p : mort) aliveIdx.push_back(p.second);
    std::function<int(int,int)> build = [&](int lo, int hi)->int {
        int node = makeNode();
        if (hi - lo == 1) {
            size_t pi = aliveIdx[lo];
            BVHNode& nd = bvhNodes[node]; nd.leaf = (int)pi; nd.minx = nd.maxx = soa_x[pi]; nd.miny = nd.maxy = soa_y[pi]; nd.mass = std::max(0,(int)soa_mass[pi]); nd.cx = soa_x[pi]; nd.cy = soa_y[pi]; return node;
        }
        int mid = (lo + hi) >> 1;
        int L = build(lo, mid);
        int R = build(mid, hi);
        BVHNode& nd = bvhNodes[node]; nd.left = L; nd.right = R;
        auto& a = bvhNodes[L]; auto& b = bvhNodes[R];
        nd.minx = std::min(a.minx, b.minx); nd.miny = std::min(a.miny, b.miny);
        nd.maxx = std::max(a.maxx, b.maxx); nd.maxy = std::max(a.maxy, b.maxy);
        nd.mass = a.mass + b.mass; if (nd.mass > 0) { nd.cx = (a.mass*a.cx + b.mass*b.cx)/nd.mass; nd.cy = (a.mass*a.cy + b.mass*b.cy)/nd.mass; } else { nd.cx=nd.cy=0; }
        return node;
    };
    int root = -1; if (!aliveIdx.empty()) root = build(0, (int)aliveIdx.size());
    auto sizeOf = [](const BVHNode& n){ return std::max(n.maxx - n.minx, n.maxy - n.miny); };
    parallelFor(soa_x.size(), [&](size_t begin, size_t end, int){
        for (size_t i=begin;i<end;++i) {
            if (!soa_alive[i]) continue; float ax=0.0f, ay=0.0f; std::function<void(int)> forceRec = [&](int idx){ if (idx<0) return; const BVHNode& nd = bvhNodes[idx]; if (nd.mass<=0) return; if (nd.leaf >= 0) { if (nd.leaf == (int)i) return; double dx=nd.cx-soa_x[i], dy=nd.cy-soa_y[i]; double r2=dx*dx+dy*dy+0.0001; double invr=1.0/std::sqrt(r2); double f=(gravityG*nd.mass)/r2; ax+=(float)(f*dx*invr); ay+=(float)(f*dy*invr); return; } double dx=nd.cx-soa_x[i], dy=nd.cy-soa_y[i]; double r2=dx*dx+dy*dy+0.0001; double s=sizeOf(nd); if ((s/std::sqrt(r2))<bhTheta){ double invr=1.0/std::sqrt(r2); double f=(gravityG*nd.mass)/r2; ax+=(float)(f*dx*invr); ay+=(float)(f*dy*invr);} else { forceRec(nd.left); forceRec(nd.right); } }; forceRec(root); soa_vx[i]+=ax*dt; soa_vy[i]+=ay*dt; }
    });
}

void PhysicsWorld::updateParticleSoA(float dts) {
    parallelFor(soa_x.size(), [&](size_t b, size_t e, int){ for (size_t i=b;i<e;++i) if (soa_alive[i]) { soa_x[i]+=soa_vx[i]*dts; soa_y[i]+=soa_vy[i]*dts; } });
}

void PhysicsWorld::handleBoundarySoA(float padding) {
    float minx = 0.0f + padding, miny = 0.0f + padding; float maxx = (float)(w - 1) - padding; float maxy = (float)(h - 1) - padding;
    // Thread-local energy accumulation to avoid data races
    std::vector<double> energyLocal((size_t)numWorkers, 0.0);
    parallelFor(soa_x.size(), [&](size_t b, size_t e, int tid){
        double eacc = 0.0;
        for (size_t i=b;i<e;++i) {
            if (!soa_alive[i]) continue;
            if (soa_x[i] < minx) { float vxb=soa_vx[i]; soa_x[i]=minx; soa_vx[i]=-soa_vx[i]*bounceRestitution; double dE=0.5*std::max(1,(int)soa_mass[i])*(vxb*vxb - soa_vx[i]*soa_vx[i]); if (dE>0) eacc += dE; }
            if (soa_x[i] > maxx) { float vxb=soa_vx[i]; soa_x[i]=maxx; soa_vx[i]=-soa_vx[i]*bounceRestitution; double dE=0.5*std::max(1,(int)soa_mass[i])*(vxb*vxb - soa_vx[i]*soa_vx[i]); if (dE>0) eacc += dE; }
            if (soa_y[i] < miny) { float vyb=soa_vy[i]; soa_y[i]=miny; soa_vy[i]=-soa_vy[i]*bounceRestitution; double dE=0.5*std::max(1,(int)soa_mass[i])*(vyb*vyb - soa_vy[i]*soa_vy[i]); if (dE>0) eacc += dE; }
            if (soa_y[i] > maxy) { float vyb=soa_vy[i]; soa_y[i]=maxy; soa_vy[i]=-soa_vy[i]*bounceRestitution; double dE=0.5*std::max(1,(int)soa_mass[i])*(vyb*vyb - soa_vy[i]*soa_vy[i]); if (dE>0) eacc += dE; }
        }
        energyLocal[(size_t)tid] += eacc;
    });
    double eSum = 0.0; for (double v : energyLocal) eSum += v; if (eSum > 0) accumulateLostEnergy(eSum);
}

void PhysicsWorld::handleCollisionsSoA() {
    // Grid (colorOfGridCell): partitions simulation space into 4 disjoint color classes so cells of the
    // same color (and their owned neighbor interactions) can be processed in parallel without touching
    // the same particles. Each cell owns its within-cell pairs and its right/down/down-right neighbor pairs.
    const float R = partRadius; const float R2 = (2*R)*(2*R);
    const float cellSize = std::max(1.0f, 2.0f * R);
    int gridW = std::max(1, (int)std::ceil((float)w / cellSize)); int gridH = std::max(1, (int)std::ceil((float)h / cellSize));
    int cells = gridW * gridH;
    if (gridW != gridWCells || gridH != gridHCells || std::fabs(cellSize - gridCellSize) > 1e-6f) { gridWCells = gridW; gridHCells = gridH; gridCellSize = cellSize; }
    auto cellIndex = [&](float x, float y) -> int { int cx=(int)std::floor(x/cellSize); int cy=(int)std::floor(y/cellSize); if (cx<0) cx=0; if (cx>=gridW) cx=gridW-1; if (cy<0) cy=0; if (cy>=gridH) cy=gridH-1; return cy*gridW + cx; };
    // Build CSR for alive particles
    csrPcell.resize(soa_x.size());
    parallelFor(soa_x.size(), [&](size_t b, size_t e, int){ for (size_t i=b;i<e;++i) csrPcell[i] = (!soa_alive[i] ? -1 : cellIndex(soa_x[i], soa_y[i])); });
    csrCounts.assign(cells, 0);
    std::vector<std::vector<int>> local; local.resize((size_t)numWorkers); for (int t=0;t<numWorkers;++t) local[(size_t)t].assign(cells, 0);
    parallelFor(soa_x.size(), [&](size_t b, size_t e, int tid){ auto& lc = local[(size_t)tid]; for (size_t i=b;i<e;++i) { int c = csrPcell[i]; if (c>=0) lc[(size_t)c]++; } });
    for (int t=0;t<numWorkers;++t) for (int c=0;c<cells;++c) csrCounts[(size_t)c] += local[(size_t)t][(size_t)c];
    csrOffsets.resize(cells+1); csrOffsets[0]=0; for (int c=0;c<cells;++c) csrOffsets[(size_t)c+1] = csrOffsets[(size_t)c] + csrCounts[(size_t)c]; csrIndex.resize((size_t)csrOffsets[(size_t)cells]);
    std::vector<int> cursor(cells); for (int c=0;c<cells;++c) cursor[(size_t)c] = csrOffsets[(size_t)c];
    for (size_t i=0;i<soa_x.size();++i) { int c = csrPcell[i]; if (c>=0) { int pos = cursor[(size_t)c]++; csrIndex[(size_t)pos] = i; } }
    // 4-color contact resolution
    auto colorOfGridCell = [&](int x,int y){ return ((x&1)<<1) | (y&1); };
    std::vector<double> energyLocal((size_t)numWorkers, 0.0);
    std::vector<std::vector<std::pair<size_t,size_t>>> fragLocal((size_t)numWorkers);
    const int neigh[3][2] = { {1,0}, {0,1}, {1,1} };
    auto within_pair = [&](size_t i, size_t j, int tid){ if (!soa_alive[i]||!soa_alive[j]) return; float dx=soa_x[i]-soa_x[j], dy=soa_y[i]-soa_y[j]; float dist2=dx*dx+dy*dy; if (dist2>R2) return; float dist=std::sqrt(std::max(dist2,1e-6f)); float nx=dx/dist, ny=dy/dist; float overlap=(2*R)-dist; if (overlap>0){ float tm = (float)std::max(1,(int)soa_mass[i]) + (float)std::max(1,(int)soa_mass[j]) + 1e-4f; float pushA = (float)std::max(1,(int)soa_mass[j]) / tm; float pushB = (float)std::max(1,(int)soa_mass[i]) / tm; soa_x[i] += nx*overlap*pushA; soa_y[i] += ny*overlap*pushA; soa_x[j] -= nx*overlap*pushB; soa_y[j] -= ny*overlap*pushB; } float vn1 = soa_vx[i]*nx + soa_vy[i]*ny; float vn2 = soa_vx[j]*nx + soa_vy[j]*ny; float rel=(soa_vx[i]-soa_vx[j])*nx + (soa_vy[i]-soa_vy[j])*ny; const float MOMENTUM_THRESHOLD = 10.0f * 100.0f; float p_rel = std::fabs((float)std::max(1,(int)soa_mass[i]) * vn1 - (float)std::max(1,(int)soa_mass[j]) * vn2); if (p_rel > MOMENTUM_THRESHOLD) fragLocal[(size_t)tid].emplace_back(i,j); if (rel > 0) return; float m1 = std::max(1,(int)soa_mass[i]); float m2 = std::max(1,(int)soa_mass[j]); float ea = std::max(0,(int)soa_elast10[i]) / 10.0f; float eb = std::max(0,(int)soa_elast10[j]) / 10.0f; float e = std::max(0.0f, std::min(1.0f, std::min(ea, eb))); float jimp = -(1 + e) * rel / (1.0f/m1 + 1.0f/m2); float jx = jimp * nx; float jy = jimp * ny; soa_vx[i] += jx / m1; soa_vy[i] += jy / m1; soa_vx[j] -= jx / m2; soa_vy[j] -= jy / m2; double mu = (m1*m2)/(m1+m2); double dE = 0.5 * mu * (1.0 - (double)e * (double)e) * (double)rel * (double)rel; if (dE>0) energyLocal[(size_t)tid] += dE; };
    for (int colorPass=0; colorPass<4; ++colorPass) {
        std::vector<int> colorCells; colorCells.reserve((size_t)cells/4+1);
        for (int y=0;y<gridH;++y) for (int x=0;x<gridW;++x) if (colorOfGridCell(x,y)==colorPass) colorCells.push_back(y*gridW + x);
        parallelFor(colorCells.size(), [&](size_t b, size_t e, int tid){ for (size_t ui=b; ui<e; ++ui) { int c = colorCells[ui]; int cx = c % gridW; int cy = c / gridW; int a0 = csrOffsets[(size_t)c], a1 = csrOffsets[(size_t)c+1]; for (int a=a0;a<a1;++a) for (int b1=a+1;b1<a1;++b1) within_pair(csrIndex[(size_t)a], csrIndex[(size_t)b1], tid); for (int k=0;k<3;++k){ int nx = cx + neigh[k][0], ny = cy + neigh[k][1]; if (nx<0||nx>=gridW||ny<0||ny>=gridH) continue; int nb = ny*gridW + nx; int b0 = csrOffsets[(size_t)nb], bE = csrOffsets[(size_t)nb+1]; for (int a=a0;a<a1;++a) for (int bb=b0; bb<bE; ++bb) within_pair(csrIndex[(size_t)a], csrIndex[(size_t)bb], tid);} } });
    }
    double eSum = 0.0; for (double v : energyLocal) eSum += v; if (eSum>0) accumulateLostEnergy(eSum);
    for (auto& v : fragLocal) for (auto& p : v) fragSchedule.push_back(p);
}

void PhysicsWorld::applyFragmentationsSoAFromSchedule() {
    if (fragSchedule.empty()) return;
    // Capacity-aware, conservative splits (2 parents -> 4 children, net +2)
    size_t liveCount = 0; for (size_t i=0;i<soa_alive.size();++i) if (soa_alive[i]) ++liveCount;
    std::vector<size_t> parentsToKill;
    parentsToKill.reserve(fragSchedule.size()*2);
    // Children buffers
    std::vector<unsigned long long> c_id;
    std::vector<float> c_x, c_y, c_vx, c_vy;
    std::vector<uint8_t> c_mass, c_elast, c_decay, c_alive;
    std::vector<char> c_sym; std::vector<int16_t> c_pix, c_piy;
    auto pushChild = [&](size_t pi, float offx, float offy, int m){
        c_id.push_back(nextId++);
        c_x.push_back(soa_x[pi] + offx);
        c_y.push_back(soa_y[pi] + offy);
        c_vx.push_back(soa_vx[pi]); c_vy.push_back(soa_vy[pi]);
        c_mass.push_back((uint8_t)std::max(1, std::min(100, m)));
        c_elast.push_back(soa_elast10[pi]);
        c_decay.push_back(0);
        c_alive.push_back(1);
        c_sym.push_back((soa_sym[pi] > 'A') ? (char)(soa_sym[pi] - 1) : 'A');
        c_pix.push_back(-1); c_piy.push_back(-1);
    };
    for (auto& pr : fragSchedule) {
        size_t i = pr.first, j = pr.second;
        if (i>=soa_alive.size()||j>=soa_alive.size()) continue;
        if (!soa_alive[i] || !soa_alive[j]) continue;
        if (soa_mass[i] < 2 || soa_mass[j] < 2) continue;
        if (liveCount + 2 > MaxParticles) continue;
        // compute ortho offset
        float dx = soa_x[i] - soa_x[j], dy = soa_y[i] - soa_y[j]; float len = std::sqrt(std::max(1e-6f, dx*dx+dy*dy));
        float nx = dx/len, ny = dy/len; float tx=-ny, ty=nx; float sep = std::max(0.1f, partRadius*0.6f);
        int am1 = soa_mass[i]/2, am2 = (int)soa_mass[i] - am1;
        int bm1 = soa_mass[j]/2, bm2 = (int)soa_mass[j] - bm1;
        parentsToKill.push_back(i); parentsToKill.push_back(j); liveCount -= 2;
        pushChild(i, tx*sep, ty*sep, am1); ++liveCount;
        pushChild(i,-tx*sep,-ty*sep, am2); ++liveCount;
        pushChild(j, tx*sep, ty*sep, bm1); ++liveCount;
        pushChild(j,-tx*sep,-ty*sep, bm2); ++liveCount;
        if (liveCount >= MaxParticles) break;
    }
    // Kill parents
    for (size_t idx : parentsToKill) soa_alive[idx] = 0;
    // Compact and append children
    if (!parentsToKill.empty() || !c_id.empty()) {
        soaCompact();
        size_t addN = c_id.size();
        size_t base = soa_id.size();
        soa_id.resize(base+addN); soa_x.resize(base+addN); soa_y.resize(base+addN);
        soa_vx.resize(base+addN); soa_vy.resize(base+addN); soa_mass.resize(base+addN);
        soa_elast10.resize(base+addN); soa_decay.resize(base+addN); soa_alive.resize(base+addN);
        soa_sym.resize(base+addN); soa_pix.resize(base+addN); soa_piy.resize(base+addN);
        for (size_t k=0;k<addN;++k) {
            soa_id[base+k]=c_id[k]; soa_x[base+k]=c_x[k]; soa_y[base+k]=c_y[k];
            soa_vx[base+k]=c_vx[k]; soa_vy[base+k]=c_vy[k]; soa_mass[base+k]=c_mass[k];
            soa_elast10[base+k]=c_elast[k]; soa_decay[base+k]=c_decay[k]; soa_alive[base+k]=c_alive[k];
            soa_sym[base+k]=c_sym[k]; soa_pix[base+k]=c_pix[k]; soa_piy[base+k]=c_piy[k];
        }
    }
    fragSchedule.clear();
}

void PhysicsWorld::updateDecaySoA() {
    // Parallel detection; centralized apply for capacity and conservation.
    // Z: 90% → 2x 'Y' (50/50); 10% → 4x 'X' (25% each) with unit velocities up/down/left/right.
    size_t liveCount = 0; for (size_t i=0;i<soa_alive.size();++i) if (soa_alive[i]) ++liveCount;
    struct ChildDec { size_t pi; float offx, offy; float vx, vy; int m; char sym; };
    std::vector<std::vector<size_t>> toKillLocal((size_t)numWorkers);
    std::vector<std::vector<ChildDec>> childrenLocal((size_t)numWorkers);
    parallelFor(soa_alive.size(), [&](size_t b, size_t e, int tid){
        auto& kill = toKillLocal[(size_t)tid]; auto& ch = childrenLocal[(size_t)tid];
        std::uniform_real_distribution<float> ang(0.0f, (float)M_PI*2.0f);
        std::uniform_int_distribution<int> coin10(0,9);
        for (size_t i=b;i<e;++i) {
            if (!soa_alive[i]) continue;
            uint8_t d = (uint8_t)(soa_decay[i] + 1); soa_decay[i] = d; if (d < 10) continue; soa_decay[i] = 0;
            if (soa_sym[i] != 'Z') continue;
            int massZ = (int)soa_mass[i];
            if (massZ < 2) { soa_sym[i] = 'Y'; continue; }
            bool fourWay = (coin10(prng) == 0) && (massZ >= 4);
            float sep = std::max(0.1f, partRadius*0.6f);
            kill.push_back(i);
            if (fourWay) {
                int q = massZ / 4; int r = massZ - 3*q; int mA=q, mB=q, mC=q, mD=q; if (r>0){mA++;r--; } if (r>0){mB++;r--; } if (r>0){mC++;r--; }
                ch.push_back({i, 0.0f, -sep, 0.0f, -1.0f, mA, 'X'});
                ch.push_back({i, 0.0f,  sep, 0.0f,  1.0f, mB, 'X'});
                ch.push_back({i, -sep, 0.0f, -1.0f, 0.0f, mC, 'X'});
                ch.push_back({i,  sep, 0.0f,  1.0f, 0.0f, mD, 'X'});
            } else {
                int m1 = massZ / 2; int m2 = massZ - m1;
                float theta = ang(prng); float ux = std::cos(theta), uy = std::sin(theta);
                ch.push_back({i, ux*sep,  uy*sep,  soa_vx[i], soa_vy[i], m1, 'Y'});
                ch.push_back({i,-ux*sep, -uy*sep, soa_vx[i], soa_vy[i], m2, 'Y'});
            }
        }
    });
    std::vector<size_t> toKill; for (auto& v : toKillLocal) for (auto x : v) toKill.push_back(x);
    std::vector<ChildDec> kids; for (auto& v : childrenLocal) { kids.insert(kids.end(), v.begin(), v.end()); }
    if (!toKill.empty() || !kids.empty()) {
        size_t allowNew = (MaxParticles > liveCount) ? (MaxParticles - liveCount) : 0;
        std::unordered_map<size_t, std::vector<ChildDec>> byParent; byParent.reserve(toKill.size()*2);
        for (auto& c : kids) byParent[c.pi].push_back(c);
        std::vector<size_t> acceptedKills; acceptedKills.reserve(toKill.size());
        std::vector<ChildDec> acceptedKids; acceptedKids.reserve(kids.size());
        for (size_t pi : toKill) {
            auto it = byParent.find(pi); if (it == byParent.end()) continue;
            size_t childCount = it->second.size(); size_t need = (childCount >= 1) ? (childCount - 1) : 0;
            if (allowNew >= need) {
                allowNew -= need; acceptedKills.push_back(pi);
                auto& vec = it->second; acceptedKids.insert(acceptedKids.end(), vec.begin(), vec.end());
            }
        }
        for (size_t idx : acceptedKills) soa_alive[idx] = 0;
        soaCompact(); size_t base = soa_id.size(); size_t addN = acceptedKids.size();
        soa_id.resize(base+addN); soa_x.resize(base+addN); soa_y.resize(base+addN); soa_vx.resize(base+addN); soa_vy.resize(base+addN); soa_mass.resize(base+addN); soa_elast10.resize(base+addN); soa_decay.resize(base+addN); soa_alive.resize(base+addN); soa_sym.resize(base+addN); soa_pix.resize(base+addN); soa_piy.resize(base+addN);
        for (size_t k=0;k<acceptedKids.size();++k) { auto& c = acceptedKids[k]; size_t dst = base+k; soa_id[dst]=nextId++; soa_x[dst]=soa_x[c.pi]+c.offx; soa_y[dst]=soa_y[c.pi]+c.offy; soa_vx[dst]=c.vx; soa_vy[dst]=c.vy; soa_mass[dst]=(uint8_t)std::max(1,std::min(100,c.m)); soa_elast10[dst]=soa_elast10[c.pi]; soa_decay[dst]=0; soa_alive[dst]=1; soa_sym[dst]=c.sym; soa_pix[dst]=-1; soa_piy[dst]=-1; }
    }
}

void PhysicsWorld::updateAdjacencyAndCombineSoA() {
    // Owner cells: each candidate pair is assigned to an owner screen cell to avoid duplicate work and data races.
    // We use the cell of the particle with the smaller (row-major) cell id as the owner. Owner cells are further
    // partitioned into 4 color classes (colorOfOwnerCell(ix,iy)) to process non-neighboring cells in parallel.
    // Use CSR grid for adjacency detection in parallel
    int gridW = w, gridH = h; int cells = gridW*gridH;
    auto cellIndex = [&](int ix, int iy)->int { ix = clampi(ix,0,w-1); iy = clampi(iy,0,h-1); return iy*gridW + ix; };
    csrPcell.resize(soa_x.size());
    parallelFor(soa_x.size(), [&](size_t b, size_t e, int){ for (size_t i=b;i<e;++i) { if (!soa_alive[i]) { csrPcell[i]=-1; continue; } int ix=(int)std::round(soa_x[i]); int iy=(int)std::round(soa_y[i]); csrPcell[i] = cellIndex(ix,iy); } });
    csrCounts.assign(cells,0);
    std::vector<std::vector<int>> local((size_t)numWorkers); for (int t=0;t<numWorkers;++t) local[(size_t)t].assign(cells,0);
    parallelFor(soa_x.size(), [&](size_t b, size_t e, int tid){ auto& lc=local[(size_t)tid]; for (size_t i=b;i<e;++i){ int c=csrPcell[i]; if (c>=0) lc[(size_t)c]++; } });
    for (int t=0;t<numWorkers;++t) for (int c=0;c<cells;++c) csrCounts[(size_t)c] += local[(size_t)t][(size_t)c];
    csrOffsets.resize(cells+1); csrOffsets[0]=0; for (int c=0;c<cells;++c) csrOffsets[(size_t)c+1]=csrOffsets[(size_t)c]+csrCounts[(size_t)c]; csrIndex.resize((size_t)csrOffsets[(size_t)cells]);
    std::vector<int> cursor(cells); for (int c=0;c<cells;++c) cursor[(size_t)c]=csrOffsets[(size_t)c]; for (size_t i=0;i<soa_x.size();++i){ int c=csrPcell[i]; if (c>=0){ int pos=cursor[(size_t)c]++; csrIndex[(size_t)pos]=i; } }
    // parallel gather present keys
    auto pkhash = [&](unsigned long long a, unsigned long long b){ unsigned long long x = (a<<1) ^ (b + 0x9e3779b97f4a7c15ULL + (a<<6) + (a>>2)); return x; };
    auto mkpair = [&](unsigned long long a, unsigned long long b){ return a<b ? std::make_pair(a,b) : std::make_pair(b,a); };
    std::vector<std::vector<PairKey>> presentLocal((size_t)numWorkers);
    parallelFor(cells, [&](size_t b, size_t e, int tid){ auto& out = presentLocal[(size_t)tid]; out.reserve(256); for (size_t ci=b; ci<e; ++ci){ int a0=csrOffsets[ci], a1=csrOffsets[ci+1]; // within
            for (int a=a0;a<a1;++a) for (int b1=a+1;b1<a1;++b1){ size_t i=csrIndex[(size_t)a], j=csrIndex[(size_t)b1]; auto pr=mkpair(soa_id[i], soa_id[j]); out.push_back({pr.first, pr.second}); }
            int cx = (int)(ci % gridW), cy = (int)(ci / gridW);
            for (int dy=-1; dy<=1; ++dy) for (int dx=-1; dx<=1; ++dx){ if (dx==0 && dy==0) continue; int nx=cx+dx, ny=cy+dy; if (nx<0||nx>=gridW||ny<0||ny>=gridH) continue; int nb = ny*gridW + nx; int b0=csrOffsets[(size_t)nb], bE=csrOffsets[(size_t)nb+1]; for (int a=a0;a<a1;++a) for (int bb=b0; bb<bE; ++bb){ size_t i=csrIndex[(size_t)a], j=csrIndex[(size_t)bb]; if (i>=j) continue; auto pr=mkpair(soa_id[i], soa_id[j]); out.push_back({pr.first, pr.second}); } }
        }
    });
    // merge: filter duplicates with a set for this frame
    std::unordered_set<unsigned long long> present;
    present.reserve(soa_id.size()*4+64);
    std::vector<PairKey> presentKeys;
    for (auto& vec : presentLocal) for (auto& k : vec) { unsigned long long h = pkhash(k.a, k.b); if (present.insert(h).second) presentKeys.push_back(k); }
    // Increment counts; prune stale
    std::unordered_set<unsigned long long> keep;
    for (auto& k : presentKeys) { adjacencyCounts.increment(k); unsigned long long hh=pkhash(k.a,k.b); keep.insert(hh); }
    std::vector<PairKey> toErase;
    adjacencyCounts.forEach([&](const PairKey& pk, int){ unsigned long long hh = pkhash(pk.a, pk.b); if (!keep.count(hh)) toErase.push_back(pk); });
    for (auto& k : toErase) adjacencyCounts.erase(k);
    // Attempt combinations
    // Build id->index map for quick lookups
    std::unordered_map<unsigned long long, size_t> id2idx; id2idx.reserve(soa_id.size()*2);
    for (size_t i=0;i<soa_id.size();++i) if (soa_alive[i]) id2idx[soa_id[i]] = i;
    // Partition keys by owner cell (cell of smaller id's index); then process per color in parallel
    auto cellOfIdx = [&](size_t i){ int ix = clampi((int)std::round(soa_x[i]), 0, w - 1); int iy = clampi((int)std::round(soa_y[i]), 0, h - 1); return iy*w + ix; };
    std::unordered_map<int, std::vector<PairKey>> byCell; byCell.reserve(presentKeys.size());
    for (auto& K : presentKeys) {
        int cnt=0; if (!adjacencyCounts.get(K, cnt)) continue; if (cnt < 10) continue;
        auto ita = id2idx.find(K.a); auto itb = id2idx.find(K.b); if (ita==id2idx.end()||itb==id2idx.end()) continue;
        size_t ia = ita->second, ib = itb->second; int ca = cellOfIdx(ia), cb = cellOfIdx(ib); int owner = (ca < cb) ? ca : cb; byCell[owner].push_back(K);
    }
    auto colorOfCell = [&](int cid){ int ix = cid % w; int iy = cid / w; return ((ix&1)<<1) | (iy&1); };
    struct Child { size_t firstIdx; float vx, vy; int mass; char sym; };
    std::uniform_int_distribution<int> coin(0,1);
    for (int pass=0; pass<4; ++pass) {
        std::vector<int> cellsList; cellsList.reserve(byCell.size());
        for (auto& kv : byCell) if (colorOfCell(kv.first) == pass) cellsList.push_back(kv.first);
        // Per-thread results
        std::vector<std::vector<size_t>> killLocal((size_t)numWorkers);
        std::vector<std::vector<Child>> childLocal((size_t)numWorkers);
        std::vector<std::vector<PairKey>> eraseLocal((size_t)numWorkers), resetLocal((size_t)numWorkers);
        parallelFor(cellsList.size(), [&](size_t b, size_t e, int tid){ auto& kill = killLocal[(size_t)tid]; auto& kids = childLocal[(size_t)tid]; auto& er = eraseLocal[(size_t)tid]; auto& rs = resetLocal[(size_t)tid]; std::unordered_set<size_t> consumed; for (size_t ui=b; ui<e; ++ui) { int cell = cellsList[ui]; auto& vec = byCell[cell]; for (auto& K : vec) {
                    int cnt=0; if (!adjacencyCounts.get(K, cnt)) continue; if (cnt < 10) continue;
                    auto ita = id2idx.find(K.a); auto itb = id2idx.find(K.b); if (ita==id2idx.end()||itb==id2idx.end()) { er.push_back(K); continue; }
                    size_t ia = ita->second, ib = itb->second; if (!soa_alive[ia]||!soa_alive[ib]) { er.push_back(K); continue; }
                    if (consumed.count(ia) || consumed.count(ib)) { rs.push_back(K); continue; }
                    if (coin(prng) == 0) { rs.push_back(K); continue; }
                    int iaL = (soa_sym[ia]-'A'+1); int ibL = (soa_sym[ib]-'A'+1); int sumL = iaL + ibL; if (sumL > 26) sumL = 26; char sym=(char)('A'+(sumL-1));
                    int rawMassSum = (int)soa_mass[ia] + (int)soa_mass[ib]; int massSum = rawMassSum; if (massSum > 100) { // record mass loss; applied centrally
                        massSum = 100; }
                    float m1 = std::max(1,(int)soa_mass[ia]); float m2 = std::max(1,(int)soa_mass[ib]); float vx = (m1*soa_vx[ia] + m2*soa_vx[ib])/(m1+m2); float vy = (m1*soa_vy[ia] + m2*soa_vy[ib])/(m1+m2);
                    size_t first = (soa_id[ia] < soa_id[ib]) ? ia : ib;
                    kill.push_back(ia); kill.push_back(ib); kids.push_back({first, vx, vy, massSum, sym}); er.push_back(K); consumed.insert(ia); consumed.insert(ib);
                } }
        });
        // Merge results and apply centrally
        std::vector<size_t> killIdx; for (auto& v : killLocal) for (auto x : v) killIdx.push_back(x);
        std::vector<Child> children; for (auto& v : childLocal) { children.insert(children.end(), v.begin(), v.end()); }
        std::vector<PairKey> toErase; for (auto& v : eraseLocal) for (auto& k : v) toErase.push_back(k);
        std::vector<PairKey> toReset; for (auto& v : resetLocal) for (auto& k : v) toReset.push_back(k);
        if (!killIdx.empty() || !children.empty()) {
            // Mark kills in parallel (no compaction yet to keep indices stable until all colors processed)
            parallelFor(killIdx.size(), [&](size_t b, size_t e, int){ for (size_t i=b;i<e;++i) { size_t idx = killIdx[i]; if (idx < soa_alive.size()) soa_alive[idx] = 0; } });
            // Append children in parallel
            size_t base = soa_id.size(); size_t addN = children.size();
            soa_id.resize(base+addN); soa_x.resize(base+addN); soa_y.resize(base+addN); soa_vx.resize(base+addN); soa_vy.resize(base+addN); soa_mass.resize(base+addN); soa_elast10.resize(base+addN); soa_decay.resize(base+addN); soa_alive.resize(base+addN); soa_sym.resize(base+addN); soa_pix.resize(base+addN); soa_piy.resize(base+addN);
            parallelFor(addN, [&](size_t b, size_t e, int){ for (size_t k=b;k<e;++k) { auto& c = children[k]; size_t dst = base + k; soa_id[dst]=nextId++; soa_x[dst]=soa_x[c.firstIdx]; soa_y[dst]=soa_y[c.firstIdx]; soa_vx[dst]=c.vx; soa_vy[dst]=c.vy; soa_mass[dst]=(uint8_t)c.mass; soa_elast10[dst]=8; soa_decay[dst]=0; soa_alive[dst]=1; soa_sym[dst]=c.sym; soa_pix[dst]=-1; soa_piy[dst]=-1; } });
        }
        for (auto& k : toErase) adjacencyCounts.erase(k);
        for (auto& k : toReset) adjacencyCounts.set(k, 0);
    }
    // Final compaction after processing all color batches
    soaCompact();
}

void PhysicsWorld::maybeSpawnFromPoolsSoA() {
    // Build occupancy set and check capacity
    size_t liveCount = 0; std::unordered_set<long long> occ; occ.reserve(soa_x.size()*2+16);
    for (size_t i=0;i<soa_x.size();++i) { if (!soa_alive[i]) continue; ++liveCount; int ix = clampi((int)std::round(soa_x[i]), 0, w - 1); int iy = clampi((int)std::round(soa_y[i]), 0, h - 1); long long key = ((long long)iy << 32) ^ (unsigned long long)ix; occ.insert(key); }
    const size_t totalCells = (size_t)w * (size_t)h; if (liveCount >= MaxParticles) return; if (occ.size() >= totalCells) return; if (massPool < MassSpawnUnit && energyPool < EnergySpawnUnit) return;
    std::uniform_int_distribution<int> ixDist(0, std::max(0, w-1)); std::uniform_int_distribution<int> iyDist(0, std::max(0, h-1)); std::uniform_real_distribution<float> ang(0.0f, (float)M_PI*2.0f);
    while (liveCount < MaxParticles && (massPool >= MassSpawnUnit || energyPool >= EnergySpawnUnit)) {
        bool found=false; int ix=0, iy=0; for (int attempt=0; attempt<200; ++attempt){ ix = ixDist(prng); iy = iyDist(prng); long long key=((long long)iy<<32) ^ (unsigned long long)ix; if (!occ.count(key)) { occ.insert(key); found=true; break; } }
        if (!found) break;
        int zMass = (int)std::floor(std::min(100.0, massPool)); if (zMass < 1) zMass = 1; double consumeM = std::min(massPool, (double)zMass); massPool -= consumeM;
        double useE = 0.0; if (energyPool >= EnergySpawnUnit) { useE = std::min(energyPool, EnergySpawnUnit); energyPool -= useE; }
        float theta = ang(prng); float ux = std::cos(theta), uy = std::sin(theta); float speed = useE>0.0 ? (float)std::sqrt(std::max(0.0, 2.0*useE / (double)std::max(1,zMass))) : 0.0f; if (speed > 5.0f) speed = 5.0f;
        // append SoA entry
        soa_id.push_back(nextId++); soa_x.push_back((float)ix); soa_y.push_back((float)iy); soa_vx.push_back(ux*speed); soa_vy.push_back(uy*speed); soa_mass.push_back((uint8_t)zMass); soa_elast10.push_back(8); soa_decay.push_back(0); soa_alive.push_back(1); soa_sym.push_back('Z'); soa_pix.push_back(-1); soa_piy.push_back(-1); ++liveCount;
        if (occ.size() >= totalCells) break;
    }
}

void PhysicsWorld::soaCompact() {
    size_t n = soa_id.size(); size_t write = 0; for (size_t i=0;i<n;++i) { if (!soa_alive[i]) continue; if (write != i) { soa_id[write]=soa_id[i]; soa_x[write]=soa_x[i]; soa_y[write]=soa_y[i]; soa_vx[write]=soa_vx[i]; soa_vy[write]=soa_vy[i]; soa_mass[write]=soa_mass[i]; soa_elast10[write]=soa_elast10[i]; soa_decay[write]=soa_decay[i]; soa_alive[write]=1; soa_sym[write]=soa_sym[i]; soa_pix[write]=soa_pix[i]; soa_piy[write]=soa_piy[i]; } ++write; }
    soa_id.resize(write); soa_x.resize(write); soa_y.resize(write); soa_vx.resize(write); soa_vy.resize(write); soa_mass.resize(write); soa_elast10.resize(write); soa_decay.resize(write); soa_alive.resize(write); soa_sym.resize(write); soa_pix.resize(write); soa_piy.resize(write);
}




#if 0
void PhysicsWorld::handleCollisions() {
    std::lock_guard<std::mutex> lock(mtx);
    const float R = partRadius;
    const float R2 = (2*R)*(2*R);
    const size_t n = particles.size();

    contacts.clear();
    contacts.reserve(n);

    // Uniform grid (linked-cell) broad-phase to prune pair checks (reused buffers)
    const float cellSize = std::max(1.0f, 2.0f * R);
    int gridW = std::max(1, (int)std::ceil((float)w / cellSize));
    int gridH = std::max(1, (int)std::ceil((float)h / cellSize));
    if (gridW != gridWCells || gridH != gridHCells || std::fabs(cellSize - gridCellSize) > 1e-6f) {
        gridWCells = gridW; gridHCells = gridH; gridCellSize = cellSize;
        broadGrid.clear(); broadGrid.resize((size_t)gridW * (size_t)gridH);
    } else {
        for (auto& v : broadGrid) v.clear();
    }
    auto cellIndex = [&](float x, float y) -> int {
        int cx = (int)std::floor(x / cellSize);
        int cy = (int)std::floor(y / cellSize);
        if (cx < 0) cx = 0; if (cx >= gridW) cx = gridW - 1;
        if (cy < 0) cy = 0; if (cy >= gridH) cy = gridH - 1;
        return cy * gridW + cx;
    };
    for (size_t i = 0; i < n; ++i) {
        auto& p = particles[i];
        if (!p.alive) continue;
        int idx = cellIndex(p.x, p.y);
        broadGrid[(size_t)idx].push_back(i);
    }

    static const int OFFS[4][2] = { {1,0}, {0,1}, {1,1}, {-1,1} };
    for (int cy = 0; cy < gridH; ++cy) {
        for (int cx = 0; cx < gridW; ++cx) {
            size_t base = (size_t)(cy * gridW + cx);
            auto& cell = broadGrid[base];
            // Pairs within the same cell
            for (size_t aidx = 0; aidx < cell.size(); ++aidx) {
                for (size_t bidx = aidx + 1; bidx < cell.size(); ++bidx) {
                    size_t i = cell[aidx], j = cell[bidx];
                    auto& a = particles[i]; auto& b = particles[j];
                    if (!a.alive || !b.alive) continue;
                    float dx = a.x - b.x; float dy = a.y - b.y;
                    float dist2 = dx*dx + dy*dy;
                    if (dist2 <= R2) {
                        float dist = std::sqrt(std::max(dist2, 1e-6f));
                        float nx = dx / dist; float ny = dy / dist;
                        float overlap = (2*R) - dist;
                        float vn1 = a.vx * nx + a.vy * ny;
                        float vn2 = b.vx * nx + b.vy * ny;
                        float rel = (a.vx - b.vx) * nx + (a.vy - b.vy) * ny;
                        contacts.push_back({i,j,nx,ny,overlap,rel,vn1,vn2});
                    }
                }
            }
            // Cross-cell pairs with 4 unique neighbor offsets to avoid duplicates
            for (auto& d : OFFS) {
                int nxCell = cx + d[0];
                int nyCell = cy + d[1];
                if (nxCell < 0 || nxCell >= gridW || nyCell < 0 || nyCell >= gridH) continue;
                size_t nb = (size_t)(nyCell * gridW + nxCell);
                auto& ncell = broadGrid[nb];
                for (size_t aidx = 0; aidx < cell.size(); ++aidx) {
                    for (size_t bidx = 0; bidx < ncell.size(); ++bidx) {
                        size_t i = cell[aidx], j = ncell[bidx];
                        auto& a = particles[i]; auto& b = particles[j];
                        if (!a.alive || !b.alive) continue;
                        float dx = a.x - b.x; float dy = a.y - b.y;
                        float dist2 = dx*dx + dy*dy;
                        if (dist2 <= R2) {
                            float dist = std::sqrt(std::max(dist2, 1e-6f));
                            float nx = dx / dist; float ny = dy / dist;
                            float overlap = (2*R) - dist;
                            float vn1 = a.vx * nx + a.vy * ny;
                            float vn2 = b.vx * nx + b.vy * ny;
                            float rel = (a.vx - b.vx) * nx + (a.vy - b.vy) * ny;
                            contacts.push_back({i,j,nx,ny,overlap,rel,vn1,vn2});
                        }
                    }
                }
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

        // Enforce strict conservation rules for collision-driven fragmentation:
        // - Both parents must have mass >= 2 (to avoid zero-mass children)
        // - Capacity must permit replacing 2 parents with 4 children (net +2)
        if (a.mass < 2 || b.mass < 2) continue; // cannot split conservatively
        if (liveCount + 2 > MaxParticles) continue; // insufficient capacity; defer fragmentation

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

        // For a: exact 50/50 split (integer), guaranteed >=1 per child
        Particle a1 = spawnChild(a, tx * sep, ty * sep);
        Particle a2 = spawnChild(a, -tx * sep, -ty * sep);
        a1.mass = a.mass / 2;
        a2.mass = a.mass - a1.mass;
        // For b: exact 50/50 split (integer)
        Particle b1 = spawnChild(b, tx * sep, ty * sep);
        Particle b2 = spawnChild(b, -tx * sep, -ty * sep);
        b1.mass = b.mass / 2;
        b2.mass = b.mass - b1.mass;

        // Capacity: we already ensured net +2 fits
        liveCount -= 2;
        newParticles.push_back(a1); ++liveCount;
        newParticles.push_back(a2); ++liveCount;
        newParticles.push_back(b1); ++liveCount;
        newParticles.push_back(b2); ++liveCount;
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
        float m1 = std::max(1, (int)a.mass);
        float m2 = std::max(1, (int)b.mass);
        float ea = std::max(0, (int)a.elasticity10) / 10.0f;
        float eb = std::max(0, (int)b.elasticity10) / 10.0f;
        float e = std::max(0.0f, std::min(1.0f, std::min(ea, eb)));
        float jimp = -(1 + e) * rel / (1.0f/m1 + 1.0f/m2);
        float jx = jimp * c.nx;
        float jy = jimp * c.ny;
        a.vx += jx / m1; a.vy += jy / m1;
        b.vx -= jx / m2; b.vy -= jy / m2;
        // Dissipated energy along normal due to inelasticity
        double mu = (m1 * m2) / (m1 + m2);
        double dE = 0.5 * mu * (1.0 - (double)e * (double)e) * (double)rel * (double)rel;
        if (dE > 0) accumulateLostEnergy(dE);
    }

    // In-place compaction of alive particles; then append new ones
    if (!deadIdx.empty() || !newParticles.empty()) {
        size_t orig = particles.size();
        size_t write = 0;
        for (size_t i = 0; i < orig; ++i) {
            if (!particles[i].alive) continue;
            if (write != i) particles[write] = particles[i];
            ++write;
        }
        particles.resize(write + newParticles.size());
        for (size_t k = 0; k < newParticles.size(); ++k) {
            particles[write + k] = newParticles[k];
        }
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
            // Z -> 2xY (exact mass conservation; 50/50 integer split)
            if (p.mass < 2) { p.symbol = 'Y'; continue; }
            // Need capacity to replace 1 with 2 (net +1)
            if (liveCount + 1 > MaxParticles) continue;
            char sym = 'Y';
            int m1 = p.mass / 2;
            int m2 = p.mass - m1;
            // Choose a random outward direction and place children on opposite sides; give them unit speed away.
            std::uniform_real_distribution<float> ang(0.0f, (float)M_PI * 2.0f);
            float theta = ang(prng);
            float ux = std::cos(theta);
            float uy = std::sin(theta);
            float sep = std::max(0.1f, partRadius * 0.6f);
            // remove parent
            if (win && p.prev_ix >= 0 && p.prev_iy >= 0) mvwaddch(win, p.prev_iy, p.prev_ix, ' ');
            p.alive = false; toRemove.push_back(i); --liveCount;
            // create children within cap
            Particle c1 = makeChild(p, sym, m1, ux*sep, uy*sep);
            c1.vx = ux * 1.0f; c1.vy = uy * 1.0f;
            toAdd.push_back(c1); ++liveCount;
            // We ensured capacity, so second child must fit
            Particle c2 = makeChild(p, sym, m2, -ux*sep, -uy*sep);
            c2.vx = -ux * 1.0f; c2.vy = -uy * 1.0f;
            toAdd.push_back(c2); ++liveCount;
        }
    }

    if (!toRemove.empty() || !toAdd.empty()) {
        // Mark removed as not alive
        std::unordered_set<size_t> rem(toRemove.begin(), toRemove.end());
        for (size_t i = 0; i < particles.size(); ++i) if (rem.count(i)) particles[i].alive = false;
        // Compact in place
        size_t orig = particles.size();
        size_t write = 0;
        for (size_t i = 0; i < orig; ++i) {
            if (!particles[i].alive) continue;
            if (write != i) particles[write] = particles[i];
            ++write;
        }
        particles.resize(write + toAdd.size());
        for (size_t k = 0; k < toAdd.size(); ++k) particles[write + k] = toAdd[k];
    }
}

void PhysicsWorld::updateAdjacencyAndCombine() {
    std::lock_guard<std::mutex> lock(mtx);
    if (particles.empty()) return;
    adjacencyCounts.reserve(particles.size() * 4);
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

    // Increment counts; remove stale entries using flat map
    std::unordered_set<unsigned long long> keep;
    for (auto& k : presentKeys) {
        PairKey key{k.a, k.b};
        adjacencyCounts.increment(key);
        unsigned long long hh = pkhash(k.a, k.b);
        keep.insert(hh);
    }
    // Collect stale and erase
    std::vector<PairKey> toErase;
    adjacencyCounts.forEach([&](const PairKey& pk, int){ unsigned long long hh = pkhash(pk.a, pk.b); if (!keep.count(hh)) toErase.push_back(pk); });
    for (auto& k : toErase) adjacencyCounts.erase(k);

    // Attempt combinations where count >= 10
    std::unordered_set<unsigned long long> consumed; // particle ids that already combined
    std::vector<size_t> removeIdx;
    std::vector<Particle> addList;
    std::uniform_int_distribution<int> coin(0,1);
    // Snapshot keys for iteration to avoid iterator invalidation
    std::vector<PairKey> keys;
    keys.reserve(adjacencyCounts.size());
    adjacencyCounts.forEach([&](const PairKey& pk, int){ keys.push_back(pk); });
    for (auto& K : keys) {
        int cnt=0; if (!adjacencyCounts.get(K, cnt)) continue;
        if (cnt >= 10) {
            unsigned long long ida = K.a, idb = K.b;
            if (consumed.count(ida) || consumed.count(idb)) { adjacencyCounts.erase(K); continue; }
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
                    int rawMassSum = A.mass + B.mass;
                    int massSum = rawMassSum;
                    if (massSum > 100) { accumulateLostMass((double)(massSum - 100)); massSum = 100; }
                    char sym = (char)('A' + (sumL - 1));
                    // Momentum-conserving velocity
                    float m1 = std::max(1, (int)A.mass);
                    float m2 = std::max(1, (int)B.mass);
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
                    adjacencyCounts.erase(K);
                    continue;
                } else {
                    adjacencyCounts.erase(K);
                    continue;
                }
            } else {
                adjacencyCounts.set(K, 0); // reset counter if chance failed
            }
        }
    }

    if (!removeIdx.empty() || !addList.empty()) {
        std::unordered_set<size_t> remset(removeIdx.begin(), removeIdx.end());
        for (size_t i = 0; i < particles.size(); ++i) if (remset.count(i)) particles[i].alive = false;
        size_t orig = particles.size();
        size_t write = 0;
        for (size_t i = 0; i < orig; ++i) {
            if (!particles[i].alive) continue;
            if (write != i) particles[write] = particles[i];
            ++write;
        }
        particles.resize(write + addList.size());
        for (size_t k = 0; k < addList.size(); ++k) particles[write + k] = addList[k];
        // Clean adjacencyCounts entries involving consumed ids
        std::vector<PairKey> toErase2;
        adjacencyCounts.forEach([&](const PairKey& pk, int& v){ if (consumed.count(pk.a) || consumed.count(pk.b)) toErase2.push_back(pk); });
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

// SoA sync helpers (for progressive refactor)
void PhysicsWorld::soaSyncFromVector() {
    std::lock_guard<std::mutex> lock(mtx);
    size_t n = particles.size();
    soa_id.resize(n);
    soa_x.resize(n); soa_y.resize(n); soa_vx.resize(n); soa_vy.resize(n);
    soa_mass.resize(n); soa_elast10.resize(n); soa_decay.resize(n); soa_alive.resize(n);
    soa_sym.resize(n); soa_pix.resize(n); soa_piy.resize(n);
    for (size_t i = 0; i < n; ++i) {
        const auto& p = particles[i];
        soa_id[i] = p.id; soa_x[i] = p.x; soa_y[i] = p.y; soa_vx[i] = p.vx; soa_vy[i] = p.vy;
        soa_mass[i] = (uint8_t)std::max(0, std::min(100, (int)p.mass));
        soa_elast10[i] = (uint8_t)std::max(0, std::min(10, (int)p.elasticity10));
        soa_decay[i] = (uint8_t)std::max(0, std::min(255, (int)p.decayTicks));
        soa_alive[i] = (uint8_t)(p.alive ? 1 : 0);
        soa_sym[i] = p.symbol;
        soa_pix[i] = (int16_t)p.prev_ix; soa_piy[i] = (int16_t)p.prev_iy;
    }
}

void PhysicsWorld::soaSyncToVector() {
    std::lock_guard<std::mutex> lock(mtx);
    size_t n = soa_id.size();
    particles.resize(n);
    for (size_t i = 0; i < n; ++i) {
        auto& p = particles[i];
        p.id = soa_id[i]; p.x = soa_x[i]; p.y = soa_y[i]; p.vx = soa_vx[i]; p.vy = soa_vy[i];
        p.mass = soa_mass[i]; p.elasticity10 = soa_elast10[i]; p.decayTicks = soa_decay[i]; p.alive = (soa_alive[i] != 0);
        p.symbol = soa_sym[i]; p.prev_ix = soa_pix[i]; p.prev_iy = soa_piy[i];
    }
}

void PhysicsWorld::accumulateLostMass(double dm) {
    if (dm > 0) massPool += dm;
}

void PhysicsWorld::accumulateLostEnergy(double de) {
    if (de > 0) energyPool += de;
}

void PhysicsWorld::maybeSpawnFromPools() {
    std::lock_guard<std::mutex> lock(mtx);
    // Count live and build occupancy of current screen cells
    size_t liveCount = 0; 
    std::unordered_set<long long> occ;
    occ.reserve(particles.size() * 2 + 16);
    for (auto& p : particles) {
        if (!p.alive) continue;
        ++liveCount;
        int ix = clampi((int)std::round(p.x), 0, w - 1);
        int iy = clampi((int)std::round(p.y), 0, h - 1);
        long long key = ((long long)iy << 32) ^ (unsigned long long)ix;
        occ.insert(key);
    }
    const size_t totalCells = (size_t)w * (size_t)h;
    if (liveCount >= MaxParticles) return;
    if (occ.size() >= totalCells) return; // screen full
    if (massPool < MassSpawnUnit && energyPool < EnergySpawnUnit) return;

    std::uniform_real_distribution<float> xDist(0.5f, std::max(0.5f, (float)w - 1.5f));
    std::uniform_real_distribution<float> yDist(0.5f, std::max(0.5f, (float)h - 1.5f));
    std::uniform_real_distribution<float> ang(0.0f, (float)M_PI * 2.0f);

    while (liveCount < MaxParticles && (massPool >= MassSpawnUnit || energyPool >= EnergySpawnUnit)) {
        // Find an empty screen cell for spawn (limited attempts)
        bool found = false;
        float x = 0.0f, y = 0.0f; int ix = 0, iy = 0;
        for (int attempt = 0; attempt < 200; ++attempt) {
            x = xDist(prng); y = yDist(prng);
            ix = clampi((int)std::round(x), 0, w - 1);
            iy = clampi((int)std::round(y), 0, h - 1);
            long long key = ((long long)iy << 32) ^ (unsigned long long)ix;
            if (!occ.count(key)) { occ.insert(key); found = true; break; }
        }
        if (!found) break; // couldn't find empty cell this pass

        int zMass = (int)std::floor(std::min(100.0, massPool));
        if (zMass < 1) zMass = 1;
        double consumeM = std::min(massPool, (double)zMass);
        massPool -= consumeM;

        double useE = 0.0;
        if (energyPool >= EnergySpawnUnit) { useE = std::min(energyPool, EnergySpawnUnit); energyPool -= useE; }

        float theta = ang(prng);
        float ux = std::cos(theta);
        float uy = std::sin(theta);
        float speed = 0.0f;
        if (useE > 0.0) {
            speed = (float)std::sqrt(std::max(0.0, 2.0 * useE / (double)std::max(1, zMass)));
            if (speed > 5.0f) speed = 5.0f;
        }

        Particle z;
        z.id = nextId++;
        z.symbol = 'Z';
        z.mass = zMass;
        z.elasticity10 = 8;
        z.x = (float)ix; z.y = (float)iy; // snap to chosen empty screen cell
        z.vx = ux * speed; z.vy = uy * speed;
        z.prev_ix = -1; z.prev_iy = -1;
        z.alive = true;
        z.decayTicks = 0;
        particles.push_back(z);
        ++liveCount;

        if (occ.size() >= totalCells) break; // no more empty cells
    }
}
#endif

void PhysicsWorld::accumulateLostMass(double dm) {
    if (dm > 0) massPool += dm;
}

void PhysicsWorld::accumulateLostEnergy(double de) {
    if (de > 0) energyPool += de;
}
void PhysicsWorld::initWorkers() {
    unsigned hc = std::thread::hardware_concurrency();
    if (hc == 0) hc = 2;
    // Leave 1 for UI thread; use at least 2 workers if possible
    numWorkers = (int)std::max(1u, hc - 1);
}

void PhysicsWorld::parallelFor(size_t n, const std::function<void(size_t,size_t,int)>& fn) {
    if (n == 0 || numWorkers <= 1) {
        fn(0, n, 0);
        return;
    }
    size_t block = (n + (size_t)numWorkers - 1) / (size_t)numWorkers;
    std::vector<std::thread> th;
    th.reserve((size_t)numWorkers);
    for (int t = 0; t < numWorkers; ++t) {
        size_t begin = (size_t)t * block;
        if (begin >= n) break;
        size_t end = std::min(n, begin + block);
        th.emplace_back([&, begin, end, t]{ fn(begin, end, t); });
    }
    for (auto& x : th) x.join();
}
