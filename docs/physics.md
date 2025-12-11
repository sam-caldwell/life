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
  - If the relative impact momentum exceeds a threshold, each parent splits into two children (next lower symbol),
    masses split roughly 50/50 with a minimum mass of 1.
  - Parents are removed; children inherit velocity and elasticity with a slight spatial separation.
  - Capacity enforced: new particles are only added up to `MaxParticles`.
- Combination (adjacency over time):
  - Pairs that remain adjacent for 10 cycles have a 50% chance to combine into a heavier particle (symbol index adds,
    capped at `Z`; mass clamped to 100). Momentum is conserved.
- Decay:
  - `Z` decays periodically into two `Y` children (mass split). Capacity enforced; if mass < 2, `Z` converts to `Y`
    without splitting.

## Tunables
- Command line:
  - `--gravity=<g>` or `-g` (non‑negative)
  - `--radius=<r>` or `-r` (non‑negative)
  - `--restitution=<e>` or `-e` (0..1)
- Environment:
  - `PHYSICS_G`, `PHYSICS_RADIUS`, `PHYSICS_RESTITUTION`
- Other constants:
  - `MaxParticles` = 100 (hard cap)
  - Color scale avoids pair 1 (black) and uses 2..16

## Performance Notes

- Collision detection is O(n²) at the current cap (100) which remains responsive on typical terminals. If you push the
  cap higher, consider adding a uniform grid / spatial hashing broad phase.

## Safety and Cleanup

- Terminal is restored on SIGINT/SIGTERM and after suspend/resume. On resume, curses is reinitialized and a full
  redraw is triggered.
