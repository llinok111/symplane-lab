// Pure geometry core for mirror-symmetry detection, surface-deviation
// certificates and constructive symmetrization.
// No three.js, no DOM, no Node APIs, no top-level side effects.
//
// Plane convention: unit normal n = [nx, ny, nz], offset d, with n . p = d.
// Reflection: p' = p - 2 (n . p - d) n.
//
// All thresholds are RELATIVE to the mesh bounding-box diagonal, never
// absolute lengths, so behaviour is invariant across geometric scale.

const GUARD = 1e-300; // denormal guard only; never used as a geometric tolerance

export function bounds(vertices) {
  const n = vertices.length / 3;
  if (n === 0) throw new Error('empty mesh');
  let minX = Infinity, minY = Infinity, minZ = Infinity;
  let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
  for (let i = 0; i < vertices.length; i += 3) {
    const x = vertices[i], y = vertices[i + 1], z = vertices[i + 2];
    if (x < minX) minX = x; if (x > maxX) maxX = x;
    if (y < minY) minY = y; if (y > maxY) maxY = y;
    if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
  }
  const size = [maxX - minX, maxY - minY, maxZ - minZ];
  const center = [(minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2];
  const diag = Math.hypot(size[0], size[1], size[2]);
  return { min: [minX, minY, minZ], max: [maxX, maxY, maxZ], size, center, diag };
}

export function normalizePlane(plane) {
  const [nx, ny, nz] = plane.n;
  const len = Math.hypot(nx, ny, nz);
  if (!(len > GUARD)) throw new Error('degenerate plane normal');
  return { n: [nx / len, ny / len, nz / len], d: plane.d / len };
}

export function planeThrough(normal, point) {
  const n = normalizePlane({ n: normal, d: 0 }).n;
  return { n, d: n[0] * point[0] + n[1] * point[1] + n[2] * point[2] };
}

export function reflectPoint(plane, p) {
  const s = plane.n[0] * p[0] + plane.n[1] * p[1] + plane.n[2] * p[2] - plane.d;
  return [p[0] - 2 * s * plane.n[0], p[1] - 2 * s * plane.n[1], p[2] - 2 * s * plane.n[2]];
}

// ---------------------------------------------------------------- covariance

export function covariance(vertices) {
  const n = vertices.length / 3;
  let cx = 0, cy = 0, cz = 0;
  for (let i = 0; i < vertices.length; i += 3) { cx += vertices[i]; cy += vertices[i + 1]; cz += vertices[i + 2]; }
  cx /= n; cy /= n; cz /= n;
  let xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
  for (let i = 0; i < vertices.length; i += 3) {
    const x = vertices[i] - cx, y = vertices[i + 1] - cy, z = vertices[i + 2] - cz;
    xx += x * x; xy += x * y; xz += x * z; yy += y * y; yz += y * z; zz += z * z;
  }
  return { centroid: [cx, cy, cz], m: [[xx / n, xy / n, xz / n], [xy / n, yy / n, yz / n], [xz / n, yz / n, zz / n]] };
}

// Cyclic Jacobi eigensolver for a real symmetric 3x3 matrix.
export function eigenSym3(m) {
  const a = [[m[0][0], m[0][1], m[0][2]], [m[1][0], m[1][1], m[1][2]], [m[2][0], m[2][1], m[2][2]]];
  const v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let sweep = 0; sweep < 24; sweep++) {
    const off = Math.abs(a[0][1]) + Math.abs(a[0][2]) + Math.abs(a[1][2]);
    const scale = Math.abs(a[0][0]) + Math.abs(a[1][1]) + Math.abs(a[2][2]) + GUARD;
    if (off <= 1e-16 * scale) break;
    for (let p = 0; p < 2; p++) {
      for (let q = p + 1; q < 3; q++) {
        const apq = a[p][q];
        if (Math.abs(apq) <= GUARD) continue;
        const tau = (a[q][q] - a[p][p]) / (2 * apq);
        const t = (tau >= 0 ? 1 : -1) / (Math.abs(tau) + Math.sqrt(tau * tau + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < 3; k++) {
          const akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - s * akq;
          a[k][q] = s * akp + c * akq;
        }
        for (let k = 0; k < 3; k++) {
          const apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - s * aqk;
          a[q][k] = s * apk + c * aqk;
        }
        for (let k = 0; k < 3; k++) {
          const vkp = v[k][p], vkq = v[k][q];
          v[k][p] = c * vkp - s * vkq;
          v[k][q] = s * vkp + c * vkq;
        }
      }
    }
  }
  const pairs = [0, 1, 2].map((k) => ({ value: a[k][k], vector: [v[0][k], v[1][k], v[2][k]] }));
  pairs.sort((p, q) => q.value - p.value);
  return pairs;
}

// ------------------------------------------------------------- spatial grids

function cellKey(ix, iy, iz) { return ix + ',' + iy + ',' + iz; }

export function makeVertexIndex(vertices, cell) {
  const map = new Map();
  const step = cell > 0 ? cell : 1;
  for (let i = 0; i < vertices.length; i += 3) {
    const k = cellKey(Math.floor(vertices[i] / step), Math.floor(vertices[i + 1] / step), Math.floor(vertices[i + 2] / step));
    let bucket = map.get(k);
    if (!bucket) { bucket = []; map.set(k, bucket); }
    bucket.push(i / 3);
  }
  return { vertices, cell: step, map, count: vertices.length / 3 };
}

export function makeTriangleIndex(vertices, faces, cell) {
  const step = cell > 0 ? cell : 1;
  const map = new Map();
  const triCount = faces.length / 3;
  for (let t = 0; t < triCount; t++) {
    const a = faces[t * 3] * 3, b = faces[t * 3 + 1] * 3, c = faces[t * 3 + 2] * 3;
    const x0 = Math.min(vertices[a], vertices[b], vertices[c]);
    const x1 = Math.max(vertices[a], vertices[b], vertices[c]);
    const y0 = Math.min(vertices[a + 1], vertices[b + 1], vertices[c + 1]);
    const y1 = Math.max(vertices[a + 1], vertices[b + 1], vertices[c + 1]);
    const z0 = Math.min(vertices[a + 2], vertices[b + 2], vertices[c + 2]);
    const z1 = Math.max(vertices[a + 2], vertices[b + 2], vertices[c + 2]);
    const ix0 = Math.floor(x0 / step), ix1 = Math.floor(x1 / step);
    const iy0 = Math.floor(y0 / step), iy1 = Math.floor(y1 / step);
    const iz0 = Math.floor(z0 / step), iz1 = Math.floor(z1 / step);
    for (let ix = ix0; ix <= ix1; ix++) for (let iy = iy0; iy <= iy1; iy++) for (let iz = iz0; iz <= iz1; iz++) {
      const k = cellKey(ix, iy, iz);
      let bucket = map.get(k);
      if (!bucket) { bucket = []; map.set(k, bucket); }
      bucket.push(t);
    }
  }
  return { vertices, faces, cell: step, map, triCount };
}

function nearestVertexSq(idx, px, py, pz) {
  const step = idx.cell, map = idx.map, verts = idx.vertices;
  const cx = Math.floor(px / step), cy = Math.floor(py / step), cz = Math.floor(pz / step);
  let best = Infinity;
  for (let r = 0; r <= 12; r++) {
    for (let ix = cx - r; ix <= cx + r; ix++) {
      for (let iy = cy - r; iy <= cy + r; iy++) {
        for (let iz = cz - r; iz <= cz + r; iz++) {
          if (Math.max(Math.abs(ix - cx), Math.abs(iy - cy), Math.abs(iz - cz)) !== r) continue;
          const bucket = map.get(cellKey(ix, iy, iz));
          if (!bucket) continue;
          for (let k = 0; k < bucket.length; k++) {
            const vi = bucket[k] * 3;
            const dx = verts[vi] - px, dy = verts[vi + 1] - py, dz = verts[vi + 2] - pz;
            const d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < best) best = d2;
          }
        }
      }
    }
    if (Number.isFinite(best) && best <= (r * step) * (r * step)) break;
  }
  if (!Number.isFinite(best) || best > (12 * step) * (12 * step)) {
    // Far query: the ring scan cannot certify a bound, so pay the exact O(V) scan.
    best = Infinity;
    for (let vi = 0; vi < verts.length; vi += 3) {
      const dx = verts[vi] - px, dy = verts[vi + 1] - py, dz = verts[vi + 2] - pz;
      const d2 = dx * dx + dy * dy + dz * dz;
      if (d2 < best) best = d2;
    }
  }
  return best;
}

// Squared distance from point p to triangle (a, b, c), following the
// Voronoi-region clamping of Ericson, "Real-Time Collision Detection" 5.1.5.
export function pointTriangleSq(px, py, pz, ax, ay, az, bx, by, bz, cx, cy, cz) {
  const abx = bx - ax, aby = by - ay, abz = bz - az;
  const acx = cx - ax, acy = cy - ay, acz = cz - az;
  const apx = px - ax, apy = py - ay, apz = pz - az;
  const d1 = abx * apx + aby * apy + abz * apz;
  const d2 = acx * apx + acy * apy + acz * apz;
  if (d1 <= 0 && d2 <= 0) return apx * apx + apy * apy + apz * apz;
  const bpx = px - bx, bpy = py - by, bpz = pz - bz;
  const d3 = abx * bpx + aby * bpy + abz * bpz;
  const d4 = acx * bpx + acy * bpy + acz * bpz;
  if (d3 >= 0 && d4 <= d3) return bpx * bpx + bpy * bpy + bpz * bpz;
  const vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) {
    const denom = d1 - d3;
    const v = denom !== 0 ? d1 / denom : 0;
    const qx = ax + abx * v, qy = ay + aby * v, qz = az + abz * v;
    return (px - qx) ** 2 + (py - qy) ** 2 + (pz - qz) ** 2;
  }
  const cpx = px - cx, cpy = py - cy, cpz = pz - cz;
  const d5 = abx * cpx + aby * cpy + abz * cpz;
  const d6 = acx * cpx + acy * cpy + acz * cpz;
  if (d6 >= 0 && d5 <= d6) return cpx * cpx + cpy * cpy + cpz * cpz;
  const vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) {
    const denom = d2 - d6;
    const w = denom !== 0 ? d2 / denom : 0;
    const qx = ax + acx * w, qy = ay + acy * w, qz = az + acz * w;
    return (px - qx) ** 2 + (py - qy) ** 2 + (pz - qz) ** 2;
  }
  const va = d3 * d6 - d5 * d4;
  if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
    const denom = (d4 - d3) + (d5 - d6);
    const w = denom !== 0 ? (d4 - d3) / denom : 0;
    const qx = bx + (cx - bx) * w, qy = by + (cy - by) * w, qz = bz + (cz - bz) * w;
    return (px - qx) ** 2 + (py - qy) ** 2 + (pz - qz) ** 2;
  }
  const denom = va + vb + vc;
  if (!(Math.abs(denom) > GUARD)) {
    // Degenerate triangle: fall back to the closest vertex.
    const da = apx * apx + apy * apy + apz * apz;
    const db = bpx * bpx + bpy * bpy + bpz * bpz;
    const dc = cpx * cpx + cpy * cpy + cpz * cpz;
    return Math.min(da, db, dc);
  }
  const v = vb / denom, w = vc / denom;
  const qx = ax + abx * v + acx * w, qy = ay + aby * v + acy * w, qz = az + abz * v + acz * w;
  return (px - qx) ** 2 + (py - qy) ** 2 + (pz - qz) ** 2;
}

export function closestPointOnMesh(idx, px, py, pz) {
  const step = idx.cell, map = idx.map, verts = idx.vertices, faces = idx.faces;
  const cx = Math.floor(px / step), cy = Math.floor(py / step), cz = Math.floor(pz / step);
  let best = Infinity;
  const scanCell = (ix, iy, iz) => {
    const bucket = map.get(cellKey(ix, iy, iz));
    if (!bucket) return;
    for (let k = 0; k < bucket.length; k++) {
      const t = bucket[k] * 3;
      const a = faces[t] * 3, b = faces[t + 1] * 3, c = faces[t + 2] * 3;
      const d2 = pointTriangleSq(px, py, pz,
        verts[a], verts[a + 1], verts[a + 2],
        verts[b], verts[b + 1], verts[b + 2],
        verts[c], verts[c + 1], verts[c + 2]);
      if (d2 < best) best = d2;
    }
  };
  for (let r = 0; r <= 48; r++) {
    for (let ix = cx - r; ix <= cx + r; ix++) {
      for (let iy = cy - r; iy <= cy + r; iy++) {
        for (let iz = cz - r; iz <= cz + r; iz++) {
          if (Math.max(Math.abs(ix - cx), Math.abs(iy - cy), Math.abs(iz - cz)) !== r) continue;
          scanCell(ix, iy, iz);
        }
      }
    }
    if (Number.isFinite(best) && best <= (r * step) * (r * step)) break;
  }
  if (!Number.isFinite(best)) {
    for (let t = 0; t < faces.length; t += 3) {
      const a = faces[t] * 3, b = faces[t + 1] * 3, c = faces[t + 2] * 3;
      const d2 = pointTriangleSq(px, py, pz,
        verts[a], verts[a + 1], verts[a + 2],
        verts[b], verts[b + 1], verts[b + 2],
        verts[c], verts[c + 1], verts[c + 2]);
      if (d2 < best) best = d2;
    }
  }
  return Math.sqrt(best);
}

// ---------------------------------------------------------------- detection

function mirrorObjective(idx, vertices, plane, diag, stride) {
  // mean squared distance from each mirrored vertex to the nearest vertex
  let sum = 0, count = 0;
  for (let i = 0; i < vertices.length; i += 3 * stride) {
    const p = [vertices[i], vertices[i + 1], vertices[i + 2]];
    const q = reflectPoint(plane, p);
    const d2 = nearestVertexSq(idx, q[0], q[1], q[2]);
    sum += d2; count++;
  }
  return sum / count / (diag * diag);
}

function orthonormalBasis(n) {
  const ax = Math.abs(n[0]), ay = Math.abs(n[1]), az = Math.abs(n[2]);
  const t = ax < ay ? (ax < az ? [1, 0, 0] : [0, 0, 1]) : (ay < az ? [0, 1, 0] : [0, 0, 1]);
  const u = [n[1] * t[2] - n[2] * t[1], n[2] * t[0] - n[0] * t[2], n[0] * t[1] - n[1] * t[0]];
  const ul = Math.hypot(...u);
  const w = [n[1] * u[2] / ul - n[2] * u[1] / ul, n[2] * u[0] / ul - n[0] * u[2] / ul, n[0] * u[1] / ul - n[1] * u[0] / ul];
  return { u: [u[0] / ul, u[1] / ul, u[2] / ul], w };
}

function refinePlane(idx, vertices, plane, diag, stride) {
  const { u, w } = orthonormalBasis(plane.n);
  const params = { a: 0, b: 0, c: 0 };
  const make = (a, b, c) => {
    const n = normalizePlane({ n: [plane.n[0] + a * u[0] + b * w[0], plane.n[1] + a * u[1] + b * w[1], plane.n[2] + a * u[2] + b * w[2]], d: 1 }).n;
    return { n, d: plane.d + c * diag };
  };
  let best = mirrorObjective(idx, vertices, plane, diag, stride);
  let stepA = 0.05, stepC = 0.02 * diag;
  for (let iter = 0; iter < 12; iter++) {
    let improved = false;
    const moves = [[stepA, 0, 0], [-stepA, 0, 0], [0, stepA, 0], [0, -stepA, 0], [0, 0, stepC], [0, 0, -stepC]];
    for (const [da, db, dc] of moves) {
      const cand = make(params.a + da, params.b + db, params.c + dc);
      const val = mirrorObjective(idx, vertices, cand, diag, stride);
      if (val < best - 1e-18) { best = val; params.a += da; params.b += db; params.c += dc; improved = true; }
    }
    if (!improved) { stepA *= 0.5; stepC *= 0.5; if (stepA < 1e-9 && stepC < 1e-9 * Math.max(diag, 1e-12)) break; }
  }
  return { plane: make(params.a, params.b, params.c), rmsNorm: best };
}

export function detectSymmetry(mesh, opts = {}) {
  // Work in a normalized frame (centroid at origin, bbox diagonal = 1) so the
  // search behaves identically at every scale; map the plane back at the end.
  const rawVertices = mesh.vertices;
  const rawBounds = bounds(rawVertices);
  const diag = rawBounds.diag;
  const tolRel = opts.tolRel ?? 0.01;
  const stride = opts.stride ?? Math.max(1, Math.ceil((rawVertices.length / 3) / 3000));
  const vertices = new Float64Array(rawVertices.length);
  for (let i = 0; i < rawVertices.length; i += 3) {
    vertices[i] = (rawVertices[i] - rawBounds.center[0]) / diag;
    vertices[i + 1] = (rawVertices[i + 1] - rawBounds.center[1]) / diag;
    vertices[i + 2] = (rawVertices[i + 2] - rawBounds.center[2]) / diag;
  }
  const { centroid, m } = covariance(vertices);
  const eigen = eigenSym3(m);
  const idx = makeVertexIndex(vertices, 1 / 32);
  const raw = [];
  for (const { value, vector } of eigen) {
    const len = Math.hypot(...vector);
    if (!(len > GUARD)) continue;
    const n = [vector[0] / len, vector[1] / len, vector[2] / len];
    const plane = { n, d: n[0] * centroid[0] + n[1] * centroid[1] + n[2] * centroid[2] };
    const rmsNorm = mirrorObjective(idx, vertices, plane, diag, stride);
    raw.push({ plane, rmsNorm, eigenValue: value, angleDeg: 0 });
  }
  raw.sort((a, b) => a.rmsNorm - b.rmsNorm);
  if (raw.length === 0) throw new Error('could not build a symmetry-plane candidate');
  const start = raw[0];
  const refined = refinePlane(idx, vertices, start.plane, 1, Math.max(stride, 2));
  const dot = Math.abs(start.plane.n[0] * refined.plane.n[0] + start.plane.n[1] * refined.plane.n[1] + start.plane.n[2] * refined.plane.n[2]);
  const candidates = raw.map((c, i) => ({ plane: i === 0 ? refined.plane : c.plane, rmsNorm: i === 0 ? refined.rmsNorm : c.rmsNorm, angleDeg: i === 0 ? Math.acos(Math.min(1, dot)) * 180 / Math.PI : 0, eigenValue: c.eigenValue }));
  const best = candidates[0];
  const rmsNorm = Math.sqrt(Math.max(0, best.rmsNorm));
  // Enforce a consistent normal orientation (largest absolute component positive)
  const n = best.plane.n;
  const flip = (n[0] !== 0 ? n[0] : n[1] !== 0 ? n[1] : n[2]) < 0 ? -1 : 1;
  const planeNorm = { n: [n[0] * flip, n[1] * flip, n[2] * flip], d: best.plane.d * flip };
  const plane = {
    n: planeNorm.n,
    d: planeNorm.d * diag + planeNorm.n[0] * rawBounds.center[0] + planeNorm.n[1] * rawBounds.center[1] + planeNorm.n[2] * rawBounds.center[2],
  };
  return {
    symmetric: rmsNorm <= tolRel,
    plane,
    rmsNorm,
    tolRel,
    diag,
    candidates: candidates.map((c) => ({ rmsNorm: Math.sqrt(Math.max(0, c.rmsNorm)), angleDeg: c.angleDeg, eigenValue: c.eigenValue })),
  };
}

// ------------------------------------------------------------- certificates

export function surfaceSamples(mesh, opts = {}) {
  const includeCentroids = opts.includeCentroids !== false;
  const includeEdges = opts.includeEdges === true;
  const v = mesh.vertices, f = mesh.faces;
  const out = [];
  for (let i = 0; i < v.length; i += 3) out.push(v[i], v[i + 1], v[i + 2]);
  if (includeCentroids) {
    for (let t = 0; t < f.length; t += 3) {
      const a = f[t] * 3, b = f[t + 1] * 3, c = f[t + 2] * 3;
      out.push((v[a] + v[b] + v[c]) / 3, (v[a + 1] + v[b + 1] + v[c + 1]) / 3, (v[a + 2] + v[b + 2] + v[c + 2]) / 3);
    }
  }
  if (includeEdges) {
    for (let t = 0; t < f.length; t += 3) {
      const a = f[t] * 3, b = f[t + 1] * 3;
      out.push((v[a] + v[b]) / 2, (v[a + 1] + v[b + 1]) / 2, (v[a + 2] + v[b + 2]) / 2);
    }
  }
  return Float64Array.from(out);
}

// One-sided (surface-sample -> mirrored surface) deviation statistics in
// absolute and diagonal-normalized units.
export function surfaceDeviation(mesh, plane, opts = {}) {
  const diag = bounds(mesh.vertices).diag;
  const cell = opts.cell ?? diag / 32;
  const idx = makeTriangleIndex(mesh.vertices, mesh.faces, cell);
  const samples = opts.samples ?? surfaceSamples(mesh, opts);
  const dists = [];
  let max = 0, over = 0;
  const tol = opts.tolRel !== undefined ? opts.tolRel * diag : Infinity;
  for (let i = 0; i < samples.length; i += 3) {
    const q = reflectPoint(plane, [samples[i], samples[i + 1], samples[i + 2]]);
    const d = closestPointOnMesh(idx, q[0], q[1], q[2]);
    dists.push(d);
    if (d > max) max = d;
    if (d > tol) over++;
  }
  dists.sort((a, b) => a - b);
  let sum = 0;
  for (const d of dists) sum += d * d;
  const pick = (p) => dists[Math.min(dists.length - 1, Math.floor(p * (dists.length - 1)))];
  return {
    count: dists.length,
    over,
    max, maxNorm: max / diag,
    rms: Math.sqrt(sum / dists.length), rmsNorm: Math.sqrt(sum / dists.length) / diag,
    p95: pick(0.95), median: pick(0.5),
    diag,
  };
}

export function certifySymmetry(mesh, plane, opts = {}) {
  const tolRel = opts.tolRel ?? 0.005;
  const stats = surfaceDeviation(mesh, plane, { ...opts, tolRel });
  return { ok: stats.over === 0, tolRel, tol: tolRel * stats.diag, samples: stats.count, ...stats };
}

export function symmetrize(mesh, plane, iterations = 2) {
  const v = Float64Array.from(mesh.vertices);
  const f = mesh.faces;
  let diag = bounds(v).diag;
  for (let it = 0; it < iterations; it++) {
    const idx = makeTriangleIndex(v, f, diag / 32);
    for (let i = 0; i < v.length; i += 3) {
      const p = [v[i], v[i + 1], v[i + 2]];
      const q = reflectPoint(plane, p);
      // nearest point on the current surface, found by a small local search
      const target = q;
      const cp = closestPointFull(idx, target[0], target[1], target[2]);
      if (Number.isFinite(cp.dist)) {
        v[i] = (p[0] + cp.point[0]) / 2;
        v[i + 1] = (p[1] + cp.point[1]) / 2;
        v[i + 2] = (p[2] + cp.point[2]) / 2;
      }
    }
    diag = bounds(v).diag;
  }
  return { vertices: v, faces: f };
}

// Closest point with coordinates (used by symmetrization).
export function closestPointFull(idx, px, py, pz) {
  const step = idx.cell, map = idx.map, verts = idx.vertices, faces = idx.faces;
  const cx = Math.floor(px / step), cy = Math.floor(py / step), cz = Math.floor(pz / step);
  let best = Infinity, bestPoint = [px, py, pz];
  const triPoint = (t) => {
    const a = faces[t * 3] * 3, b = faces[t * 3 + 1] * 3, c = faces[t * 3 + 2] * 3;
    const ax = verts[a], ay = verts[a + 1], az = verts[a + 2];
    const bx = verts[b], by = verts[b + 1], bz = verts[b + 2];
    const cx3 = verts[c], cy3 = verts[c + 1], cz3 = verts[c + 2];
    const d2 = pointTriangleSq(px, py, pz, ax, ay, az, bx, by, bz, cx3, cy3, cz3);
    if (d2 < best) {
      best = d2;
      bestPoint = triangleClosestPoint(px, py, pz, ax, ay, az, bx, by, bz, cx3, cy3, cz3);
    }
  };
  for (let r = 0; r <= 48; r++) {
    for (let ix = cx - r; ix <= cx + r; ix++) for (let iy = cy - r; iy <= cy + r; iy++) for (let iz = cz - r; iz <= cz + r; iz++) {
      if (Math.max(Math.abs(ix - cx), Math.abs(iy - cy), Math.abs(iz - cz)) !== r) continue;
      const bucket = map.get(cellKey(ix, iy, iz));
      if (!bucket) continue;
      for (let k = 0; k < bucket.length; k++) triPoint(bucket[k]);
    }
    if (Number.isFinite(best) && best <= (r * step) * (r * step)) break;
  }
  if (!Number.isFinite(best)) {
    for (let t = 0; t < faces.length; t += 3) triPoint(t / 3);
  }
  return { dist: Math.sqrt(best), point: bestPoint };
}

// Closest point coordinates with the same region logic as pointTriangleSq.
function triangleClosestPoint(px, py, pz, ax, ay, az, bx, by, bz, cx, cy, cz) {
  const abx = bx - ax, aby = by - ay, abz = bz - az;
  const acx = cx - ax, acy = cy - ay, acz = cz - az;
  const apx = px - ax, apy = py - ay, apz = pz - az;
  const d1 = abx * apx + aby * apy + abz * apz;
  const d2 = acx * apx + acy * apy + acz * apz;
  if (d1 <= 0 && d2 <= 0) return [ax, ay, az];
  const bpx = px - bx, bpy = py - by, bpz = pz - bz;
  const d3 = abx * bpx + aby * bpy + abz * bpz;
  const d4 = acx * bpx + acy * bpy + acz * bpz;
  if (d3 >= 0 && d4 <= d3) return [bx, by, bz];
  const vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) {
    const denom = d1 - d3; const v = denom !== 0 ? d1 / denom : 0;
    return [ax + abx * v, ay + aby * v, az + abz * v];
  }
  const cpx = px - cx, cpy = py - cy, cpz = pz - cz;
  const d5 = abx * cpx + aby * cpy + abz * cpz;
  const d6 = acx * cpx + acy * cpy + acz * cpz;
  if (d6 >= 0 && d5 <= d6) return [cx, cy, cz];
  const vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) {
    const denom = d2 - d6; const w = denom !== 0 ? d2 / denom : 0;
    return [ax + acx * w, ay + acy * w, az + acz * w];
  }
  const va = d3 * d6 - d5 * d4;
  if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
    const denom = (d4 - d3) + (d5 - d6); const w = denom !== 0 ? (d4 - d3) / denom : 0;
    return [bx + (cx - bx) * w, by + (cy - by) * w, bz + (cz - bz) * w];
  }
  const denom = va + vb + vc;
  if (!(Math.abs(denom) > GUARD)) return [ax, ay, az];
  const v = vb / denom, w = vc / denom;
  return [ax + abx * v + acx * w, ay + aby * v + acy * w, az + abz * v + acz * w];
}
