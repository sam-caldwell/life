# Life (Agent-Based) — Details

This document explains the operation, rules, and inner workings of `life`.

## Overview

Each automaton is an independent agent running in its own thread. Agents occupy cells on a 2D grid (the terminal).
They sense neighbors, decide actions locally, and request moves/interactions through a central `Board` which applies
rules atomically and updates the UI.

## Controls

- `s` start/pause, `p` pause
- `r` reseed, `c` clear
- `-`/`+` slower/faster
- `q` quit

## Rendering

- Species letter: `A`–`Z`.
- Weight (0–100) mapped to color buckets with a legend on the bottom bar; 100 draws bold.

## Core Rules

- Movement: one cardinal step per move. Movement probability decreases as weight rises.
- Flee: if a neighbor could eat the actor, the actor attempts to move to the best available cardinal cell maximizing
  distance from predators.
- Eat:
  - Cross‑species: allowed if actor species letter is higher (e.g., `B` can eat `A`).
  - Same‑species: allowed only if the species letter is `A`–`R` and the actor’s weight is < 10. (`S`–`Z` never eat
    their own species.)
  - On success, the actor gains the target’s weight (capped at 100). Target is removed.
- Push: if heavier than target, push target one step directly away (if destination empty). If blocked, apply a recoil
  rule that moves the target back up to `k = (actor.weight + blocker.weight − target.weight)` steps if space allows.
- Spawn: requires the actor to have eaten recently and the partner’s weight to be within ±1 of actor’s weight.
  Newborn appears in a random empty adjacent cell. A spawn backoff TTL prevents rapid reattempts when space is
  unavailable.
- Hunger/weight:
  - If the actor goes multiple cycles without eating, it periodically loses weight (rate scales with load factor).
  - Weight ≤ 0: die and remove.
  - If weight remains 100 for > 100 cycles: die and remove (prevents permanent saturation).

## Concurrency Model

- Each automaton runs in a thread with periodic sleeps based on the configured delay (see tuning).
- All grid mutations happen in `Board` under a mutex; drawing is incremental and thread‑safe.
- Exceptions inside automaton threads are caught; the automaton is removed and an error is logged.

## Tuning

Environment variables:

- `LIFE_THREADS_PER_CORE` — max threads per CPU core (default 20; also bounded by grid size)
- `LIFE_MAX_AUTOMATA` — absolute cap (optional)
- `LIFE_STEP_DELAY_MS` — per‑automaton sleep (default 250; min 5)
- `LOG_LEVEL` — `debug|info|warn|error|none` (default `info`)

Keyboard at runtime: `-` and `+` adjust delay.

## Status Line

- Shows weight legend, thread count (live/cap), delay, and RUNNING/PAUSED.

## Safety and Cleanup

- Signal handling restores the terminal on SIGINT/SIGTERM and supports suspend/resume (Ctrl‑Z / SIGTSTP and SIGCONT).
- On resume, the app reinitializes curses and triggers a full redraw.
