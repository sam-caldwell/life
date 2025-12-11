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
  - bounce with coefficient of restitution `restitution`.
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
  - `Z` decays periodically into two `Y` children. Mass is conserved exactly with a 50/50 integer split (m1 = ⌊M/2⌋,
    m2 = M − m1). If a split would exceed capacity (needs net +1 slot), the decay is deferred. If mass < 2, `Z` converts
    to `Y` without splitting.
- Recycling:
  - Mass that cannot be instantiated due to capacity or clamps, and energy dissipated by inelastic impacts/boundary
    bounces, accumulates in hidden pools. When the mass pool reaches 100 or the energy pool reaches 200 (arbitrary
    units), the world injects a new random `Z` particle using available mass and a velocity derived from the energy pool.

## Tunables
- Command line:
  - `--gravity=<g>` or `-g` (non‑negative)
  - `--radius=<r>` or `-r` (non‑negative)
  - `--restitution=<e>` or `-e` (0..1)
- Environment:
  - `PHYSICS_G`, `PHYSICS_RADIUS`, `PHYSICS_RESTITUTION`
- Other constants:
  - `MaxParticles` = 200 (hard cap)
  - Color scale avoids pair 1 (black) and uses 2..16

## Performance Notes

- Collision detection uses a uniform grid (linked-cell) broad phase with cell size ≈ `2*radius`. Only pairs within the
  same or neighboring cells are tested, making the typical complexity near O(n + k) where k is the number of close
  pairs. In the worst case (all particles in one cell) it degrades toward O(n²), but this is rare with reasonable
  radii and distributions.

## Safety and Cleanup

- Terminal is restored on SIGINT/SIGTERM and after suspend/resume. On resume, curses is reinitialized and a full
  redraw is triggered.
