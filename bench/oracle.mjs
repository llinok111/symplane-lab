// Independent oracle: closest point on a triangle by projecting onto the
// triangle plane, then clamping with explicit Voronoi-region tests written
// separately from src/core.mjs, and no spatial index (all faces scanned).
// Own tolerance constant, not shared with the main implementation.
const ORACLE_EPS = 1e-18;

function sub(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
function len2(a) { return dot(a, a); }

export function triangleDistance(p, a, b, c) {
  const ab = sub(b, a), ac = sub(c, a), ap = sub(p, a);
  const n = [ab[1] * ac[2] - ab[2] * ac[1], ab[2] * ac[0] - ab[0] * ac[2], ab[0] * ac[1] - ab[1] * ac[0]];
  const nn = len2(n);
  if (nn > ORACLE_EPS) {
    // Project onto the triangle plane and test whether the projection lies
    // inside using three edge half-spaces.
    const t = dot(sub(p, a), n) / nn;
    const q = [p[0] - t * n[0], p[1] - t * n[1], p[2] - t * n[2]];
    const inside = (u, v, w) => {
      const uv = sub(v, u), uw = sub(w, u), uq = sub(q, u);
      const c1 = uv[0] * uq[1] - uv[1] * uq[0];
      const c2 = uv[0] * uw[1] - uv[1] * uw[0];
      return c1 * c2 >= -ORACLE_EPS;
    };
    if (inside(a, b, c) && inside(b, c, a) && inside(c, a, b)) return Math.sqrt(len2(sub(p, q)));
  }
  // Not inside (or degenerate): minimum over the three edges.
  const edge = (u, v) => {
    const uv = sub(v, u), up = sub(p, u);
    const t = dot(up, uv) / (len2(uv) || ORACLE_EPS);
    const s = Math.max(0, Math.min(1, t));
    return Math.sqrt(len2([p[0] - (u[0] + uv[0] * s), p[1] - (u[1] + uv[1] * s), p[2] - (u[2] + uv[2] * s)]));
  };
  return Math.min(edge(a, b), edge(b, c), edge(c, a));
}

export function bruteSurfaceDistance(mesh, px, py, pz) {
  let best = Infinity;
  const v = mesh.vertices, f = mesh.faces;
  for (let t = 0; t < f.length; t += 3) {
    const a = f[t] * 3, b = f[t + 1] * 3, c = f[t + 2] * 3;
    const d = triangleDistance([px, py, pz],
      [v[a], v[a + 1], v[a + 2]], [v[b], v[b + 1], v[b + 2]], [v[c], v[c + 1], v[c + 2]]);
    if (d < best) best = d;
  }
  return best;
}
