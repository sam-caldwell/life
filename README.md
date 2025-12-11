# Life + Physics (ncurses)

Two small terminal simulations rendered with ncurses:

- `life`: An agent-based, threaded “life-like” simulation where autonomous actors move, eat, push, and spawn.
- `physics`: A Newtonian particle toy with collisions, fragmentation, combination, and simple gravity.

Both apps draw a legend and status bar at the bottom and handle Ctrl‑C/Ctrl‑Z gracefully (terminal is restored).

## Dedication

This project is dedicated to Claude Shannon, Ludwig Boltzmann, and all the other minds who saw science as a playground
rather than a career; who encouraged us to tinker and explore rather than memorize. May we all skip an occasional
class like Einstein, tinker like Shannon, and laugh like Maxwell.

## Requirements

- CMake 3.16+
- C++17 compiler (Clang/GCC/MSVC)
- ncurses (Linux/macOS)
  - macOS: `brew install ncurses` (if missing)
  - Ubuntu/Debian: `sudo apt-get install libncurses5-dev`
  - Windows: use WSL with ncurses

## Build

- Using the provided Makefile: `make clean configure build`
- Builds both `life` and `physics` in `./build/` directory.

## Run
- To Run Life (agent-based): `./build/life`
- To Run Physics (particles): `./build/physics`

## Tips

- Resize your terminal larger for more room; the grid uses all rows minus one (for the status bar).
- If rendering looks odd after suspend/resume, the app will force a full redraw on resume.
- Log verbosity is controlled by `LOG_LEVEL` (`debug|info|warn|error|none`). Default is `info`.

## Further Reading
- [Life Internals](docs/life.md) : Documentation about the `life` program and how it works.
- [Physics Internals](docs/physics.md) : Documentation about the `physics` program and how it works.

## License

MIT — see `LICENSE.txt`.
