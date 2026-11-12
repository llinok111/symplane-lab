// Deterministic procedural meshes for the symmetry lab. Pure functions only.
// Symmetric shape: radius modulated by |x| (invariant under x -> -x) with
// y^2 / y*z terms that destroy every other mirror plane, so the intended
// symmetry plane is unique. Anisotropic squash keeps that plane but removes
// the spherical degeneracy of the covariance spectrum.

export function makeRand(seed) {
  let a = seed >>> 0;
  return function rand() {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function symmetricRadius(u) {
  const ax = Math.abs(u[0]);
  return 1 + 0.16 * Math.cos(3 * ax) - 0.10 * u[1] * u[1] + 0.07 * u[1] * u[2] + 0.05 * ax * u[1] * u[2];
}

function asymmetricRadius(u) {
  // No mirror plane: odd terms in x and z plus a mixed y*z term break every
  // candidate reflection, which is what the abstention test needs.
  return 1 + 0.22 * u[0] + 0.15 * (u[1] * u[1] - u[2] * u[2])
    + 0.12 * u[0] * u[1] + 0.18 * u[2] + 0.10 * u[1] * u[2];
}

function signedVolume(vertices, faces) {
  let vol = 0;
  for (let t = 0; t < faces.length; t += 3) {
    const a = faces[t] * 3, b = faces[t + 1] * 3, c = faces[t + 2] * 3;
    const ax = vertices[a], ay = vertices[a + 1], az = vertices[a + 2];
    const bx = vertices[b], by = vertices[b + 1], bz = vertices[b + 2];
    const cx = vertices[c], cy = vertices[c + 1], cz = vertices[c + 2];
    vol += ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
  }
  return vol / 6;
}

// UV-sphere topology (pole + rings + pole) with a radial modulation.
export function buildSphereMesh(radiusFn, nu = 64, nv = 32) {
  const verts = [0, radiusFn([0, 1, 0]), 0];
  const ringStart = (iv) => 1 + (iv - 1) * nu;
  for (let iv = 1; iv < nv; iv++) {
    const theta = Math.PI * iv / nv;
    const st = Math.sin(theta), ct = Math.cos(theta);
    for (let iu = 0; iu < nu; iu++) {
      const phi = 2 * Math.PI * iu / nu;
      const u = [st * Math.cos(phi), ct, st * Math.sin(phi)];
      const r = radiusFn(u);
      verts.push(u[0] * r, u[1] * r, u[2] * r);
    }
  }
  const south = verts.length / 3;
  verts.push(0, -radiusFn([0, -1, 0]), 0);

  const faces = [];
  for (let iu = 0; iu < nu; iu++) {
    const a = ringStart(1) + iu, b = ringStart(1) + (iu + 1) % nu;
    faces.push(0, b, a);
  }
  for (let iv = 1; iv < nv - 1; iv++) {
    for (let iu = 0; iu < nu; iu++) {
      const a = ringStart(iv) + iu, b = ringStart(iv) + (iu + 1) % nu;
      const c = ringStart(iv + 1) + iu, d = ringStart(iv + 1) + (iu + 1) % nu;
      faces.push(a, b, c, b, d, c);
    }
  }
  for (let iu = 0; iu < nu; iu++) {
    const a = ringStart(nv - 1) + iu, b = ringStart(nv - 1) + (iu + 1) % nu;
    faces.push(south, a, b);
  }
  if (signedVolume(verts, faces) < 0) {
    for (let i = 0; i < faces.length; i += 3) { const t = faces[i]; faces[i] = faces[i + 1]; faces[i + 1] = t; }
  }
  return { vertices: verts, faces };
}

export function rotatePoint(p, euler) {
  const [rx, ry, rz] = euler;
  let [x, y, z] = p;
  // X
  let c = Math.cos(rx), s = Math.sin(rx);
  let y1 = c * y - s * z, z1 = s * y + c * z;
  y = y1; z = z1;
  // Y
  c = Math.cos(ry); s = Math.sin(ry);
  const x1 = c * x + s * z, z2 = -s * x + c * z;
  x = x1; z = z2;
  // Z
  c = Math.cos(rz); s = Math.sin(rz);
  const x2 = c * x - s * y, y2 = s * x + c * y;
  x = x2; y = y2;
  return [x, y, z];
}

export function transformMesh(mesh, { euler = [0, 0, 0], scale = 1, translate = [0, 0, 0], squash = [1, 1, 1] } = {}) {
  const out = new Array(mesh.vertices.length);
  for (let i = 0; i < mesh.vertices.length; i += 3) {
    const q = rotatePoint([mesh.vertices[i] * squash[0], mesh.vertices[i + 1] * squash[1], mesh.vertices[i + 2] * squash[2]], euler);
    out[i] = q[0] * scale + translate[0];
    out[i + 1] = q[1] * scale + translate[1];
    out[i + 2] = q[2] * scale + translate[2];
  }
  return { vertices: out, faces: mesh.faces.slice() };
}

// Apply a known asymmetry of `amplitude` to a deterministic subset of vertices.
export function perturbMesh(mesh, { count = 61, amplitude = 0.03, seed = 777 } = {}) {
  const rand = makeRand(seed);
  const n = mesh.vertices.length / 3;
  const verts = mesh.vertices.slice();
  const chosen = [];
  // Deterministic distinct indices, away from the poles (no special reason
  // besides keeping the topology/watertightness untouched).
  while (chosen.length < Math.min(count, n - 8)) {
    const i = 4 + Math.floor(rand() * (n - 8));
    if (!chosen.includes(i)) chosen.push(i);
  }
  for (const i of chosen) {
    const theta = rand() * Math.PI * 2, z = rand() * 2 - 1;
    const r = Math.sqrt(1 - z * z);
    const d = [r * Math.cos(theta), z, r * Math.sin(theta)];
    verts[i * 3] += d[0] * amplitude;
    verts[i * 3 + 1] += d[1] * amplitude;
    verts[i * 3 + 2] += d[2] * amplitude;
  }
  return { vertices: verts, faces: mesh.faces.slice(), perturbedVertices: chosen.length };
}

// Ground-truth symmetry plane after transformMesh: local plane x = 0 maps to
// (R * e_x) . p = (R * e_x) . translate.
export function truthPlane({ euler = [0, 0, 0], translate = [0, 0, 0] } = {}) {
  const n = rotatePoint([1, 0, 0], euler);
  return { n, d: n[0] * translate[0] + n[1] * translate[1] + n[2] * translate[2] };
}

export function buildSymmetric({ nu = 64, nv = 32, perturb = 61, amplitude = 0.03, seed = 777 } = {}) {
  const base = buildSphereMesh(symmetricRadius, nu, nv);
  const withSquash = transformMesh(base, { squash: [1, 1.25, 0.8] });
  if (!perturb) return { ...withSquash, perturbedVertices: 0 };
  return perturbMesh(withSquash, { count: perturb, amplitude, seed });
}

export function buildAsymmetric({ nu = 48, nv = 24 } = {}) {
  return transformMesh(buildSphereMesh(asymmetricRadius, nu, nv), { squash: [1.1, 0.9, 1] });
}
