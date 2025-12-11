# Nuclear (Atoms) — Details

This document describes the rules and inner workings of `nuclear`.

## Overview

Terminal simulation rendering atoms with printable ASCII. Elements H (Z=1) through Pu (Z=94) map one‑to‑one to
printable ASCII characters `!`..`~`. Particles represent atoms with minimal properties. The engine models Newtonian
motion with simplified electromagnetism (Coulomb) and simplified special relativity (speed cap and proper‑time). 
Radioactive decay is approximate and capacity‑aware.

No background thread: the main loop steps the world based on elapsed time with a fixed‑dt accumulator.

## Controls

- `s` start/pause, `p` pause
- `r` reseed half‑capacity (unique cells), `c` clear
- `-`/`+` slower/faster
- `q` quit
- `e` toggle electromagnetism on/off
- `o` toggle overlay (Z/group vs charge)
- `h` toggle decay on/off
- `k`/`K` increase/decrease Coulomb strength (k)
- `]`/`[` increase/decrease LBVH theta (θ)
- `}`/`{` increase/decrease softening epsilon (ε)

## Rendering

- Mapping: `!`..`~` → H..Pu (Z=1..94). The status line shows “Map: !..~→H..Pu”.
- Overlays:
  - Z/group: color buckets by atomic number (periodic‑like grouping).
  - Charge: positive charges warm colors; negative charges cool; neutral white.

## Physics and Rules

- Integration:
  - Euler update per substep (`x += vx*dt`, `y += vy*dt`).
  - SR‑lite speed cap: `|v| < c_light` (scaled cells/s).
- Electromagnetism (simplified):
  - Coulomb force only: `F = k * q1 * q2 * r̂ / (r^2 + ε^2)` with softening `ε` to avoid singularities.
  - Acceleration is `a = F / m` using a scaled mass bucket.
  - Long‑range approximation: LBVH broad‑phase (Barnes–Hut style) with acceptance threshold `θ`.
  - Magnetism is intentionally omitted in this first version; a simplified v×B approximation may be added later.
- Boundaries:
  - Bounce from walls with coefficient of restitution `restitution`.
- Inter‑particle collisions:
  - Not modeled in this variant (Coulomb dominates); overlap resolution may be added in a later pass.
- Ionization/recombination:
  - Not modeled in this version; future work may allow threshold‑based `q±1` changes using local field/impact energy.
- Decay (approximate):
  - Each atom keeps a proper‑time accumulator `τ` advanced by `dt/γ` (`γ = 1/sqrt(1−β²)`, `β=|v|/c`).
  - Half‑life tiers by atomic number (toy model):
    - Z ≤ 82 (Pb): stable for our purposes.
    - 83..88: 60 s; 89..92: 30 s; 93..94: 10 s.
  - When `τ ≥ half_life`, atom decays to Z→Z−1 (min H), glyph updated, `τ` resets. No radiation products are modeled.
  - Decay can be toggled on/off; behavior is capacity‑agnostic (no new particles created).

## Initialization & Reseed

- Startup seeding populates exactly half of `MaxParticles` (cap‑aware) at unique screen cells.
- Elements are sampled uniformly from H...Pu; initial charges sampled from {−1,0,+1}; velocities are small random 
  values.

## Tunables
- Command line: none yet; use in‑app keys to adjust.
- In‑app keys:
  - `e` EM enable, `k/K` Coulomb strength (k)
  - `]`/`[` LBVH theta (θ)
  - `}`/`{` softening (ε)
  - `o` overlay mode
  - `h` decay enable
- Other constants:
  - `MaxParticles` = 300 (hard cap)
  - `c_light` ≈ 50 cells/s (velocity cap)
  - `restitution` for boundary bounces

## Performance Notes

- LBVH for Coulomb reduces all‑pairs O(n²) to roughly O(n log n) with adjustable `θ` trade‑off.
- Softening `ε` stabilizes near‑field interactions; tune alongside `k` for visual clarity.
- Parallel per‑particle traversal updates velocities across worker threads.
- SoA storage; buffer reuse; incremental redraws for responsiveness.

## Safety and Cleanup

- Terminal is restored on exit and exceptions. If needed in future, suspend/resume handlers can be added like in
  `physics`.

