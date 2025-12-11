# Physics (Particles) — Details

This document describes the rules and inner workings of `physics`.

## Overview

2D Newtonian particle toy rendered with ncurses. Particles have mass, velocity, elasticity, collide with each other and
with boundaries, can fragment on high‑momentum impacts, combine after sustained adjacency, and occasionally decay.

No background thread: the main loop steps the world based on elapsed time.

## Controls

- `s` start/pause, `p` pause
- `r` reseed random, `c` clear
- `-`/`+` slower/faster
- `q` quit

## Rendering

- Mass: 1–100. Color buckets mapped to pairs 2..16 (avoids black); higher buckets render bold.
- Legend at bottom shows bucket thresholds.

## Physics and Rules

- Integration: 
  - simple Euler update per step (`x += vx*dt`, `y += vy*dt`).
- Boundary:
  - Elastic: particles bounce with coefficient of restitution `restitution`.
  - Inelastic wall splits: if a particle is sufficiently inelastic (elasticity ≤ 0.4) and hits a wall, it splits into
    2–4 children (next lower symbol) when capacity allows. Mass is divided equally (integer 50/50/25% shares with any
    remainder distributed by +1), and each child receives an equal share of the parent’s kinetic energy. Children are
    emitted in distinct directions consistent with the wall normal/tangent (away from the wall and spread along the
    tangent); at corners, children are emitted into different outward quadrants. If capacity is insufficient to
    instantiate all children (net +N−1 slots), the parent falls back to an inelastic bounce.
- Collisions (pairwise check):
  - Contact detected when centers within `2*radius`.
  - Resolve overlap by pushing apart based on relative masses.
  - Normal impulse uses minimum per‑particle elasticity (0..10 → 0..1).
- Fragmentation (high momentum):
  - If the relative impact momentum exceeds a threshold, each parent splits into two children (next lower symbol) with
    exact mass conservation using a 50/50 integer split (m1 = ⌊M/2⌋, m2 = M − m1). A parent with mass < 2 does not
    fragment.
  - Fragmentation proceeds only when capacity allows replacing two parents with four children (net +2). Otherwise, the
    pair resolves as a regular collision and fragmentation is deferred.
  - Parents are removed; children inherit velocity and elasticity with a slight spatial separation.
- Combination (adjacency over time):
  - Pairs that remain adjacent for 10 cycles have a 50% chance to combine into a heavier particle (symbol index adds,
    capped at `Z`; mass clamped to 100). Momentum is conserved.
- Decay:
  - Letters `W`, `X`, `Y`, `Z` decay periodically. For a parent letter `L ∈ {W,X,Y,Z}`:
    - 90%: split into two children of letter `(L−1)` with exact 50/50 integer mass split. If mass < 2, the parent simply
      converts to `(L−1)` (no split).
    - 10%: split into four children of letter `(L−2)` with integer 25% shares; children are emitted along unit cardinal
      directions. Requires mass ≥ 4; otherwise the two-child mode applies.
    - Capacity-aware: if a split would exceed capacity (needs net +1 or +3 slots), the decay is deferred for that parent.
- Recycling:
  - Mass that cannot be instantiated due to capacity or clamps, and energy dissipated by inelastic impacts/boundary
    bounces, accumulates in hidden pools. When the mass pool reaches 100 or the energy pool reaches 200 (arbitrary
    units), the world injects a new random `Z` particle using available mass and a velocity derived from the energy pool.

## Initialization & Reseed

- Startup seeding: populates exactly half of `MaxParticles` (cap-aware) at unique screen cells.
- Symbol range: initial and reseed populations use only upper-half letters `M`..`Z`.

## Tunables
- Command line:
  - `--gravity=<g>` or `-g` (non‑negative)
  - `--radius=<r>` or `-r` (non‑negative)
  - `--restitution=<e>` or `-e` (0..1)
- Environment:
  - `PHYSICS_G`, `PHYSICS_RADIUS`, `PHYSICS_RESTITUTION`
- Other constants:
  - `MaxParticles` = 300 (hard cap)
  - Color scale avoids pair 1 (black) and uses 2..16

## Performance Notes

- Collision detection uses a uniform grid (linked-cell) broad phase with cell size ≈ `2*radius`. Only pairs within the
  same or neighboring cells are tested, making the typical complexity near O(n + k) where k is the number of close
  pairs. In the worst case (all particles in one cell) it degrades toward O(n²), but this is rare with reasonable
  radii and distributions.

## Performance Improvements

- Build optimization: Release builds use `-O3 -DNDEBUG -march=native` with link-time optimization (LTO) and fast-math
  enabled when supported.
- Gravity: Replaced per-particle/cluster summation with a Barnes–Hut quadtree for O(n log n) far-field approximation.
- Memory reuse: Reuses the collision broad-phase grid and contact buffers each frame to avoid reallocations.
- In-place compaction: Removes dead particles by compacting the array in place and appends new particles without
  constructing throwaway vectors.
- Hot map reserves: `adjacencyCounts` pre-reserves capacity to reduce rehashing churn during combine checks.
- Fixed time step: Simulation advances with a fixed dt using an accumulator, improving numerical stability and CPU
  predictability under load.
- Rendering I/O: Skips redraw for particles that did not move screen cell this frame to reduce terminal writes.

## Safety and Cleanup

- Terminal is restored on SIGINT/SIGTERM and after suspend/resume. On resume, curses is reinitialized and a full
  redraw is triggered.
