# Life — OOP, Threaded, ncurses

## Description
This project is an object‑oriented, agent‑based take on Conway’s Game of Life. In the classic Game of Life, a 
two‑dimensional grid evolves in discrete time steps where each cell is either alive or dead. The next state of each 
cell depends only on its eight neighbors: a live cell survives if it has two or three live neighbors, otherwise it 
dies; a dead cell becomes alive if it has exactly three live neighbors (Gardner, 1970). This simulation extends the 
idea by representing each “life form” as an autonomous object and thread that can push, eat, or spawn based on local 
interactions. The simulation enforces a dynamic population cap of 200 automata per CPU core on the host.

## Requirements

- CMake 3.16+
- A C++17-compatible compiler (GCC/Clang/MSVC)
- `ncurses` (Linux/macOS). 
  - On macOS: 
    - `brew install ncurses` if missing; o
  - On Ubuntu:
    `sudo apt-get install libncurses5-dev`.
  - On Windows:
    - Install Linux Subsystem for Windows (WSL) and install `ncurses` there.
      I am not inclined to support the Redmond Mafia, and I'm not writing a bunch of .NET code.

## Build

- Configure: `make configure`
- Build: `make build`
- Run classic life: `./build/life`
- Run physics variant: `./build/physics`

## Controls

- `s`: start/pause toggle
- `p`: pause
- `r`: reseed (10% fill)
- `c`: clear board
- `-`/`+`: slow down / speed up automata
- `q`: quit

The physics variant uses the same controls and displays a mass legend (0–100) mapped to 16 color buckets.

## Physics Variant Rules

- Only Newtonian mechanics are applied (no Life/agent rules).
- Particles have mass [0–100], velocity, and collide elastically.
- Gravitational acceleration is computed from nearby mass clusters (adjacent particles combine as a single gravitational source).
- Particles bounce off screen boundaries.
- Population is capped at 100 particles.

### Particle Elasticity and Energy Loss

- Each particle has an `elasticity` in the range 0–10 (random at spawn).
- On particle–particle collision, the pair’s coefficient of restitution `e` is derived from their elasticities
  (0..10 → 0..1), combined conservatively via `e = min(e_a, e_b)`. Momentum is conserved and kinetic energy is
  reduced according to `e`, with the loss interpreted as heat.
- Boundary bounces use the global `restitution` setting.

### Fragmentation (High-Momentum Collisions)

- If two particles collide with sufficient momentum (greater than `10 × max_mass`, i.e., > 1,000 using max mass 100),
  each particle splits into two particles of the next lower letter (e.g., `Z+Z → 4×Y`, `Z+Y → 2×Y + 2×X`).
- Each child inherits half the parent’s mass (rounded), parent velocity (slightly separated), and the parent’s elasticity.
- The original two particles are removed; only the produced particles remain.

### Combination (Adjacency Over Time)

- If two particles remain immediately adjacent for 10 cycles, they have a 50% chance to combine.
- Letter combination follows index addition (e.g., `A+A=B`, `B+B=D`, `A+B=C`), capped at `Z`.
- Combined mass is clamped to the maximum mass (100). Momentum is conserved: the resulting velocity is the
  mass-weighted average of the two velocities, and the new particle appears at the first particle’s position.

### Decay

- Z decays every 10 cycles into two Y particles. The original Z disappears; only the two Ys remain. Each child
  inherits approximately half the mass and the parent’s elasticity. Each child is given unit speed (1) and heads in
  opposite random directions away from the original location. If the particle cap (100) would be exceeded, only one Y
  may be produced.

### Initial Population

- On launch, the physics variant seeds a random number of particles between 1 and 20 (subject to the 100-particle cap).

### Physics Tunables

- Defaults: `gravity=9.8`, `radius=1.0`, `restitution=0.98` (applies to collisions and boundary bounces).
- Flags (override env):
  - `./build/physics --gravity=12.0 --radius=0.8 --restitution=0.95`
  - Short forms: `-g 12.0 -r 0.8 -e 0.95`
- Environment variables:
  - `PHYSICS_G=12.0`
  - `PHYSICS_RADIUS=0.8`
  - `PHYSICS_RESTITUTION=0.95`
- Validation: `gravity>=0`, `radius>=0`, `0<=restitution<=1`. Flags take precedence over environment variables.

## Behavior

### Initial State
- Every automaton starts with a random weight in [10, 50], symbol `A`–`Z`, and color. Species are assigned with a bias
  where `A` is most likely and `Z` least likely; at least two `Z` automata are guaranteed on reseed if space allows.
  Weight then evolves over time (0–100) due to eating and starvation.
- Species is indicated by the letter (A–Z). 
- Weight is indicated by ANSI color: darker for low weight up to white at 100. A weight legend appears on the bottom 
  bar as 0-10-20-30-50-70-80-90-100 (100 is shown in bold white).

### General Movement
- Movement is restricted to one cardinal step (up/down/left/right) per move.
- Movement probability decreases as weight approaches 100 (heavier moves less).

### Actions
- Flee:
  - If an automaton is near an automaton which could eat it, it attempts to move to an adjacent empty cell maximizing 
    distance from the predator in order to survive.
- Eat: 
  - Allowed if the actor’s species letter is greater than the target’s (e.g., B can eat A; A cannot eat B).
  - Same-species eating is allowed only if the actor’s weight is < 10 and the species letter is A–R. Species S–Z cannot eat their own kind.
  - On success, the actor gains the target’s weight (capped at 100).
- Push: 
  - If heavier, 
    - push a neighbor one cell farther away from the actor (destination must be free).
  - If lighter, 
    - Pushing the neighbor will result in a recoil and the automata will be moved back in the opposite direction.
    - Pushing a neighbor into another neighbor will result in recoil equal to the difference between the actor's weight
      and the sum of the two neighbors' weight and the weight of the automata into which the neighbor was pushed.
    - If the recoil would cause the actor to collide with another automata, the actor will be destroyed.
- Spawn: 
  - permitted only if the invoker has eaten within the last 10 cycles and the peer’s weight is equal or within ±1 
    of the invoker.
  - The newborn automaton appears adjacent to the invoker.
- Hunger: 
  - If an automaton goes 10 cycles without eating, it loses 1 weight.
  - When weight reaches 0, it dies and is removed.
- Max weight: 
  - If an automaton remains at weight 100 for more than 100 cycles, it dies and is removed.
  - Each automaton runs in its own thread and sleeps according to the current delay shown on the status line.

## References

Gardner, M. (1970). Mathematical games: The fantastic combinations of John Conway’s new
     solitaire game “Life.” Scientific American, 223(4), 120–123.
