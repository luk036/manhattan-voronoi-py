# AGENTS.md

## Architecture

Dual-language L1 (Manhattan distance) Voronoi diagram library implementing Lee & Wong's divide-and-conquer algorithm.

| Directory | Language | Purpose |
|---|---|---|
| `src/voronoi.js` | JS | Core library (Lee & Wong algorithm) |
| `dist/voronoi.js` | JS | Built library (Babel ES2015 → ES5) |
| `main.js` | JS | Browser demo (SVG rendering) |
| `build/build.js` | JS | Built demo bundle (Browserify + Babel) |
| `index.html` | HTML | Demo page entrypoint |
| `python_ai/voronoi.py` | Python | Direct translation of `src/voronoi.js` |
| `python_ai/demo.py` | Python | Matplotlib equivalent of `main.js` |

## Commands

### JavaScript

```bash
npm install                     # install dev dependencies (gulp, babel, browserify)
npx gulp build                  # build library → dist/voronoi.js
npx gulp build-library          # build library only
npx gulp watch                  # watch + rebuild on change
```

Open `index.html` in a browser to see the demo (loads `build/build.js`).

### Python

```bash
python python_ai/voronoi.py     # smoke test (4 + 16 random sites)
python python_ai/demo.py        # generate python_ai/voronoi_demo.png (64 sites, matplotlib)
```

Python dependencies: `matplotlib` only (no requirements.txt).

## Gotchas

### JS → Python translation pitfalls

These caused bugs during translation and must be preserved:

- **`R.sort()` is in-place mutation.** JS mutates `R` and the mutated order affects the parent `recursiveSplit`'s tree structure. Python `sorted(R)` creates a copy — use `R.sort(key=...)` instead.
- **`cropLArray` sort is descending.** JS uses `angle(B) - angle(A)` (reverse sort). Python equivalent: `key=..., reverse=True`.
- **`checkForOphans` sort picks opposite extreme.** `goUp=True` → largest extreme first; `goUp=False` → smallest extreme first. Python key: `-extremeA if goUp else extremeA`.
- **Division by zero in `findL1Bisector`.** JS `x/0` returns `Infinity` without crashing. Python raises `ZeroDivisionError`. Move slope/intercept computation **after** the `abs(xDistance) == 0` early-return guard.
- **`Array.reverse()` mutates + returns reference.** JS `points = arr.reverse()` mutates the original array. Python `[::-1]` creates a copy. For the polygon-building phase the difference is harmless (only affects rendering direction of merge lines, not geometry), but be aware.
- **`Array.reduce()` with complex accumulator.** Manual `for` loops needed instead — Python's `functools.reduce` doesn't match JS reduce behavior with multi-type accumulators (dict → list → dict).

### JS-specific

- Babel transpiles ES6 → ES5 via `babel-preset-es2015`. The source uses `let`/`const`, arrow functions, template literals, and spread operators.
- `build/build.js` is a Browserify bundle of `main.js` + `src/voronoi.js`. The bundle exposes ES module exports globally.
- No test framework configured. `npm test` will fail with "no test specified".

### Python-specific

- Matplotlib renders with y-axis 0 at bottom. The voronoi algorithm computes with y=0 at top (SVG convention). Demo uses `height - y` to flip.
- No virtualenv/config file — install matplotlib manually.
- `distance()` has an `isinstance` check for dict vs raw point; the JS version always receives raw `[x,y]` arrays.
