## Maximal Visual Madness (mvm)

./mvm -I '/path/to/anim' -o '/path/to/out' -p terr_ -n 25 --kf-hold 140 --ease ease-in-out \
  --stars --stars-layer1 300,1.2 --stars-layer2 120,2.4 --stars-seed 42

./mvm -I '/path/to/your/in' -o '/path/to/out' -p test_ -n 25 --kf-hold 140 --ease ease-in-out \
  --stars --stars-layer1 400,0.6 --stars-layer2 200,1.2 --stars-foreground


./mvm -I '/Users/ulli/Documents/Pandi-Wolve-Love/!clients/Terrassenbad/anim' \
  -o '/Users/ulli/Documents/Pandi-Wolve-Love/!clients/Terrassenbad/out' \
  -p terrsf2_ -n 25 --kf-hold 140 --ease ease-in-out \
  --stars --stars-layer1 400,0.6 --stars-layer2 200,1.2 --stars-foreground \
  -t 8 --frame-threads 4


Audio-reactive and keyframe-interpolation renderer that outputs PNG image sequences ready for Adobe After Effects.

Two workflows are supported:
- Classic audio-driven rendering (generates frames from an audio file using simple DSP features)
- Folder-driven image interpolation (sorts images alphabetically, interpolates between consecutive images, and writes sequentially numbered PNGs)

### Features
- Alphabetical folder-based keyframe interpolation
- Pre-/post-hold of keyframes for readability in titling
- Easing modes (linear, ease-in, ease-out, ease-in-out)
- Optional Gaussian blur passes, including mid-only blur for a retro/demo-scene feel
- Zero-padded 9-digit filenames suitable for direct import as an image sequence in After Effects

### Dependencies
- C++17 compiler
- libpng
- libsndfile
- fmt

On macOS (Homebrew):
```bash
brew install libpng libsndfile fmt
```

### Build
```bash
clang++ -std=c++17 -O2 -o mvm main.cpp -lpng -lsndfile -lfmt
```

### Usage
Run without `--input-dir` to use the classic audio-driven workflow (uses defaults compiled in `main.cpp`).
Use `--input-dir` to enable folder-driven interpolation mode.

```bash
./mvm [options]
```

#### Folder-driven interpolation options
These options are active when `--input-dir` is provided.

- `--input-dir, -I <dir>`: Directory containing keyframe images. Files are read and sorted alphabetically. Supported: .png
- `--out, -o <dir>`: Output directory (default: `/Users/ulli/Documents/mvm/`). PNGs are written here.
- `--prefix, -p <name>`: Output filename prefix (default: `seq_`). Final filenames are `<prefix><9-digit-frame>.png`.
- `--frames, -n <num>`: Number of frames per transition (default: `50`).
- `--start, -S <num>`: Starting frame number for the sequence (default: `1`).
- `--hold, -H <num>`: Post-hold; repeat the destination keyframe for `<num>` extra frames after each transition (default: `0`).
- `--in-hold <num>`: Pre-hold; repeat the source keyframe for `<num>` frames before each transition (default: `0`).
- `--ease <mode>`: Interpolation easing. One of `linear`, `ease-in`, `ease-out`, `ease-in-out` (default: `linear`).
- `--blur <num>`: Apply `<num>` passes of 3x3 Gaussian blur during interpolation (default: `0`).
- `--blur-mid-only`: Only apply blur near the midpoint of the transition.
- `--blur-mid-width <0..0.5>`: Half-width around the midpoint where blur applies (default: `0.2`). Effective only with `--blur`.
- `--threads, -t <num>`: Worker threads used for image operations (blur, rotation, alpha set, scans). Default: `1`.
- `--frame-threads <num>`: Number of frames to render concurrently per transition. Default: `1`.
- `--rotate`: Apply a simple rotation over the duration of each transition.
- `--kf-hold <num>`: Hold duration in frames per keyframe block. When set, each keyframe is shown for `<num>` frames before transitioning; the last keyframe can be held separately.
- `--first-kf-hold <num>`: Override hold for the first keyframe only.
- `--last-kf-hold <num>`: Override hold for the last keyframe only.
- Starfield overlay (parallax):
  - `--stars`: Enable parallax starfield overlay.
  - `--stars-layer1 <count,speed>`: Layer 1 star count and horizontal speed in px/frame (e.g., `200,0.7`).
  - `--stars-layer2 <count,speed>`: Layer 2 star count and horizontal speed in px/frame (e.g., `80,1.4`).
  - `--stars-seed <num>`: Seed for deterministic star positions.
  - `--stars-size <px>`: Base star size scaling factor (default: `3`). Actual star size derives from layer speed and jitter.
  - `--stars-foreground`: Draw stars after image interpolation (ignoring background gating) to ensure visibility.
  - `--stars-twinkle`: Enable per-star alpha twinkle.
  - `--stars-twinkle-period <frames>`: Twinkle cycle length (default: `48`).
  - `--stars-trails <N>`: Draw N trailing steps behind stars (additive).
  - `--stars-trail-fade <0..1>`: Trail attenuation per step (default: `0.75`).
  - `--stars-yoffsets v1,v2,...`: Per-layer vertical offsets (e.g., for VU-like motion).
  - `--stars-deflect`: Deflect stars near title pixels for interaction.
  - `--stars-deflect-mode <attract|repel>`: Direction of deflection (default: `attract`).
  - `--stars-deflect-range <px>`: Sampling radius around text pixels (default: `24`).
  - `--stars-deflect-strength <0..1>`: Strength of deflection (default: `0.5`).
  - `--stars-deflect-occlude`: Don’t draw stars over non-background (default on). Use `--no-stars-deflect-occlude` to disable.
- `--help, -h`: Print help and exit.

Notes:
- Adjacent images must be the same resolution. Mismatched pairs are skipped with a message.
- Output frame numbers are zero-padded to 9 digits for stable import in NLE/VFX tools.
- Frame numbering accounts for pre-hold, transition frames, and post-hold before starting the next transition.

#### Classic (audio-driven) mode
If `--input-dir` is omitted, the program runs the classic audio-driven renderer using the defaults currently set in `main.cpp` (audio filename, output path, fps, resolution, etc.). This path is useful for generating algorithmic keyframes from audio content. Future versions may expose these as CLI flags.

### Examples
- Basic interpolation between images with a readable end hold:
```bash
./mvm -I ./keyframes -o ./out -p title_ -n 48 -H 60
```

- Title cards with pre-roll, eased motion, and subtle mid-only blur:
```bash
./mvm -I ./titles -o ./out -p demo_ -n 48 -S 1 \
  --in-hold 12 --hold 60 \
  --ease ease-in-out \
  --blur 1 --blur-mid-only --blur-mid-width 0.2 \
  -t 8 --frame-threads 4
```

- Linear transitions, no blur, simple end hold:
```bash
./mvm -I /path/to/keys -o /path/to/out -p seq_ -n 40 -H 40
```

### Performance notes
- `--threads` parallelizes per-frame image work (blur, rotation, alpha, and internal scans).
- `--frame-threads` renders multiple frames of a transition at once; keep modest (e.g., 2–4) to balance CPU, memory, and disk I/O.
- Output is deterministic with a fixed merge order; per-thread buffers are reduced consistently.

### Adobe After Effects import
- File → Import → File…
- Select the first PNG in the sequence and enable "Import as: PNG Sequence".
- Set the sequence frame rate to match your project (e.g., 25 fps)
- Each sequence block (pre-hold + transition + post-hold) forms a continuous segment per keyframe pair.

### Tips for titling and demo-scene aesthetics
- Readability: Hold the final keyframe for 1–3 seconds at your timeline fps (e.g., 25–75 frames at 25 fps).
- Motion feel: Use `--ease ease-in-out` and a small `--blur` with `--blur-mid-only` for a subtle analog vibe.
- Consistency: Keep all keyframe source images the same size.


