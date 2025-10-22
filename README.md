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
- Run: `./build/life`

## Controls

- `s`: start/pause toggle
- `p`: pause
- `r`: reseed (10% fill)
- `c`: clear board
- `-`/`+`: slow down / speed up automata
- `q`: quit

## Behavior

### Initial State
- Every automaton starts with a random weight in [10, 50], symbol `A`–`Z`, and color. Weight then evolves over 
  time (0–100) due to eating and starvation.
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
