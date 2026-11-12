// Baseline: brute-force plane search on a fixed Fibonacci-normal grid with a
// fixed offset grid, scored with the same mirror-to-nearest-vertex objective.
// This is the "just scan the parameter space" approach the PCA + local search
// is compared against.
import { reflectPoint } from '../src/core.mjs';

function nearestVertexSq(vertices, px, py, pz) {
  let best = Infinity;
  for (let j = 0; j < vertices.length; j += 3) {
    const dx = vertices[j] - px, dy = vertices[j + 1] - py, dz = vertices[j + 2] - pz;
    const d2 = dx * dx + dy * dy + dz * dz;
    if (d2 < best) best = d2;
  }
  return best;
}

export function fibonacciNormals(count) {
  const out = [];
  const golden = Math.PI * (3 - Math.sqrt(5));
  for (let i = 0; i < count; i++) {
    const y = 1 - (i / (count - 1)) * 2;
    const r = Math.sqrt(Math.max(0, 1 - y * y));
    const phi = golden * i;
    out.push([r * Math.cos(phi), y, r * Math.sin(phi)]);
  }
  return out;
}

export function bruteForcePlane(mesh, opts = {}) {
  const normals = opts.normals ?? 36;
  const offsets = opts.offsets ?? 7;
  const spanRel = opts.spanRel ?? 0.25;
  const v = mesh.vertices;
  const diag = opts.diag;
  const stride = opts.stride ?? 6;
  let cx = 0, cy = 0, cz = 0;
  for (let i = 0; i < v.length; i += 3) { cx += v[i]; cy += v[i + 1]; cz += v[i + 2]; }
  const n = v.length / 3; cx /= n; cy /= n; cz /= n;
  let best = { plane: null, rmsNorm: Infinity };
  for (const nrm of fibonacciNormals(normals)) {
    const base = nrm[0] * cx + nrm[1] * cy + nrm[2] * cz;
    for (let k = 0; k < offsets; k++) {
      const d = base + (k / (offsets - 1) - 0.5) * 2 * spanRel * diag;
      const plane = { n: nrm, d };
      let sum = 0, count = 0;
      for (let i = 0; i < v.length; i += 3 * stride) {
        const q = reflectPoint(plane, [v[i], v[i + 1], v[i + 2]]);
        sum += nearestVertexSq(v, q[0], q[1], q[2]);
        count++;
      }
      const rmsNorm = Math.sqrt(sum / count) / diag;
      if (rmsNorm < best.rmsNorm) best = { plane, rmsNorm };
    }
  }
  return best;
}
