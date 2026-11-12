# symplane-lab

Find the mirror-symmetry plane of an arbitrary triangle mesh, **certify** the residual on the surface, and **constructively symmetrize** the mesh — all in a pure Node core that the browser re-runs live.

Measured promise: on a rotated, scaled (2.4×), perturbed 1,986-vertex mesh the plane normal is recovered to **0.00396°** (baseline grid search: **7.659°**, i.e. **1,933× coarser**) while detection is **1.43× faster** (247 ms vs 354 ms median). The certificate then quantifies the asymmetry as **max 0.9410% of the bbox diagonal over 5,954 surface samples**; six symmetrization iterations reduce that to **0.1115% (8.44×)** and the 0.5%-tolerance certificate flips from FAIL to PASS.

Upstream: none (original tool; the mesh generator and all algorithms are in this repository).

## What it does

- `detectSymmetry(mesh)` — covariance eigenvectors seeded plane, then a vertex-hash mirror-residual refinement. The target is the reflection plane `n · p = d`. Works in a normalized frame (centroid at origin, bbox diagonal = 1), so the same shape at any scale gives the same verdict.
- `certifySymmetry(mesh, plane, { tolRel })` — reflects **vertices + triangle centroids**, measures the distance from each mirrored sample to the nearest point **on the mesh surface** (uniform triangle grid), and reports `max`, `rms`, `p95`, count over tolerance. A shape with no symmetry abstains (`symmetric:false`) instead of guessing.
- `symmetrize(mesh, plane, iterations)` — every vertex is moved halfway to the closest point of the mirrored surface, iterated. The result is re-certified by the same function (no trust in the repair).
- `surfaceDeviation` / `verify` — same measurement used by the CLI, the browser and the tests.

## How to run (one command each)

```bash
node tools/gen-sample.mjs     # deterministic: writes data/{scene,exact,asym,ref}.json
node bench/verify.mjs         # 8 assertions, exit 0 only if all pass
node bench/bench.mjs          # 5+ timed repetitions -> bench/result.json
python3 -m http.server 4173   # then open http://localhost:4173
```

No `npm install`, no bundler: `index.html` pins `three@0.169.0` via an import map. Drag the three sliders to move the mirror plane; the overlay shows the live max deviation from `src/core.mjs` and the node-baked value for the **same 397 vertex samples**. At slider 0/0/0 the two numbers match bit-for-bit.

## Evidence (real runs, `bench/result.json`)

| What | Result | Baseline / comparison |
| --- | --- | --- |
| Plane normal error | **0.003963°** (rms 0.1984% diag) | 36×7 Fibonacci+offset grid: **7.659°** |
| Exact symmetric mesh | **1.21e-6°**, rms 2.0e-13 | — |
| Asymmetric mesh verdict | **abstain** (`symmetric:false`, rms 1.93% > 1% tol) | would be a false positive if forced |
| Detection time | **247.0 ms** median, cold 262.5 ms, p95 248.8, RSD 1.09% (5 samples) | grid baseline 354.0 ms |
| Surface certificate | **5.87 ms** median for 120 samples (7 reps) | all-faces oracle **91.5 ms** → **15.6×** |
| Certificate agreement | grid max **0.002339798** vs oracle max **0.002339798** | identical on all 120 samples |
| Asymmetry before repair | max **0.9410%** diag, rms 0.0850%, 47/5,954 samples over 0.5% | — |
| After 6 symmetrization steps | max **0.1115%** diag (8.44× lower), certificate **PASS** | repair 2 iterations 135.8 ms |
| Scale invariance | 1e-6, 1e-3, 1, 1e3, 1e6: rms spread **0.0000%**, angle spread **1.2e-10°**, all five verdicts `true` | fixed by frame normalization |
| Oracle independence | `bench/oracle.mjs` scans **all** faces, has its own projection/clamping code and its own epsilon; agreement **8.7e-16** on 161 samples | — |

Run `node bench/verify.mjs` yourself: it prints every number above (8 passed / 0 failed).

## 本轮创新

Most "symmetry check" scripts are a visual aid or a plane-finder that never says *how far* the mesh is from symmetric. This tool turns the mirror plane into a **measurable certificate**: (1) a covariance-eigenvector seed plus a refinement that minimizes the mirrored-vertex residual *in a normalized frame*, (2) a **surface**-based (not vertex-only) deviation measured against the triangle mesh, reported in bbox-diagonal units, and (3) a **constructive** step whose output is re-certified by the same independent measurement. The same `src/core.mjs` runs in Node and in the browser, and the page shows the node-baked number next to the live one for the same sample set — so the interactive view cannot drift from the evidence. Nothing here depends on three.js; the renderer only draws the core's output.

## Known limitations (read these)

- **Samples, not the full surface.** Deviations are measured at vertices and triangle centroids (5,954 points for this mesh). Triangle interiors are not covered by a proof; a badly warped single triangle could pass. This is a sampled certificate, not a continuous Hausdorff bound.
- **Tolerance is relative but fixed.** `tolRel = 0.005` means 0.5% of the bbox diagonal; it is scale-invariant, but it is still a single threshold — a mesh with legitimate 1% design asymmetry will be rejected. The detector's own threshold is `1%` rms of the diagonal, also disclosed.
- **The eigenvector seed assumes the symmetry plane normal is a covariance eigenvector.** That is exact for a reflection-symmetric shape, and the local refinement recovers up to ~0.01 rad of noise; a shape whose asymmetry is comparable to its size (our asymmetric sample: 1.93% rms) is refused rather than forced.
- **Degenerate covariance.** A sphere-like shape has an isotropic covariance matrix, so *any* plane is a candidate; the tool will pick one and the certificate will decide. It does not claim uniqueness of the plane in that case.
- **Performance.** The vertex-hash residual is O(V) per objective evaluation; refinement does 12×6 evaluations on a 1/2 vertex subsample. On meshes above ~50k vertices this becomes slow (no BVH); the certificate grid is O(samples × nearby triangles).
- **Benchmark variance.** The certificate median is ~5.9 ms; RSD is 15.5% because the whole measurement is only a few milliseconds and GC/JIT noise dominates. Absolute spread is 4.6–6.9 ms. The detection RSD is 1.09%.
- **Baseline is deliberately simple**: a 36-normal × 7-offset grid scored on a 1/6 vertex subsample. It is a reference for "scan the parameter space", not a tuned competitor.

## Layout

```
index.html / main.js      renderer + sliders (consumes core, no geometry math)
src/core.mjs              pure ES module: detection, certificate, symmetrization
tools/meshgen.mjs         deterministic mesh generator (seeded)
tools/gen-sample.mjs      writes data/*.json (committed, < 1 MB total)
bench/verify.mjs          8 assertions over real data
bench/bench.mjs          timing + bench/result.json
bench/oracle.mjs         independent all-faces distance oracle
baseline/naive.mjs       grid-search baseline
```

MIT © 2026 llinok111. Data files are generated by `tools/gen-sample.mjs`; rerunning it reproduces byte-identical JSON.
