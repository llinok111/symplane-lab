// Timing harness: cold start + >=5 warm repetitions, median/p95/stddev.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { detectSymmetry, certifySymmetry, symmetrize, bounds, surfaceSamples, makeTriangleIndex, closestPointOnMesh, reflectPoint } from '../src/core.mjs';
import { bruteForcePlane } from '../baseline/naive.mjs';
import { bruteSurfaceDistance } from './oracle.mjs';

const here = path.dirname(fileURLToPath(import.meta.url));
const scene = JSON.parse(fs.readFileSync(path.join(here, '..', 'data', 'scene.json'), 'utf8'));
const diag = bounds(scene.vertices).diag;

function stats(ms) {
  const sorted = [...ms].sort((a, b) => a - b);
  const mean = ms.reduce((a, b) => a + b, 0) / ms.length;
  const sd = Math.sqrt(ms.reduce((a, b) => a + (b - mean) ** 2, 0) / ms.length);
  const pick = (p) => sorted[Math.min(sorted.length - 1, Math.floor(p * (sorted.length - 1)))];
  return {
    samples: ms.length,
    values: ms.map((x) => Number(x.toFixed(3))),
    min: Number(sorted[0].toFixed(3)),
    median: Number(pick(0.5).toFixed(3)),
    p95: Number(pick(0.95).toFixed(3)),
    stddev: Number(sd.toFixed(3)),
    rsd: Number((sd / mean).toFixed(4)),
  };
}
function timeOnce(fn) { const t = performance.now(); fn(); return performance.now() - t; }

// Cold start is the very first call in this process.
const coldMs = timeOnce(() => detectSymmetry(scene));
detectSymmetry(scene); // warmup
const detTimes = [];
let det = null;
for (let i = 0; i < 5; i++) detTimes.push(timeOnce(() => { det = detectSymmetry(scene); }));
const detection = { coldMs: Number(coldMs.toFixed(3)), ...stats(detTimes) };

const tBase = performance.now();
const baseline = bruteForcePlane(scene, { diag, normals: 36, offsets: 7, stride: 6 });
const baselineMs = performance.now() - tBase;
const angleOf = (plane) => Math.acos(Math.min(1, Math.abs(
  plane.n[0] * scene.truth.n[0] + plane.n[1] * scene.truth.n[1] + plane.n[2] * scene.truth.n[2]))) * 180 / Math.PI;

// Certificate: grid-accelerated vs all-faces oracle on the same 120 samples.
const samples = surfaceSamples(scene);
// Fair comparison on the SAME 120 surface samples: grid-acceelerated core vs
// the all-faces oracle. Totals over 5954 samples are recorded separately.
const strideSamples = Math.max(1, Math.floor(samples.length / 3 / 120));
const subset = new Float64Array(120 * 3);
for (let n = 0; n < 120; n++) {
  const i = n * strideSamples * 3;
  subset[n * 3] = samples[i]; subset[n * 3 + 1] = samples[i + 1]; subset[n * 3 + 2] = samples[i + 2];
}
const certTimes = [];
let cert = null;
certifySymmetry(scene, det.plane, { tolRel: 0.005, samples: subset }); // warmup
for (let i = 0; i < 7; i++) certTimes.push(timeOnce(() => { cert = certifySymmetry(scene, det.plane, { tolRel: 0.005, samples: subset }); }));
let oracleMax = 0;
const tOracle = performance.now();
for (let i = 0; i < subset.length; i += 3) {
  const q = reflectPoint(det.plane, [subset[i], subset[i + 1], subset[i + 2]]);
  oracleMax = Math.max(oracleMax, bruteSurfaceDistance(scene, q[0], q[1], q[2]));
}
const oracleMs = performance.now() - tOracle;

const repairMs = timeOnce(() => symmetrize(scene, det.plane, 2));

const out = {
  generatedAt: new Date().toISOString(),
  node: process.version,
  mesh: { vertices: scene.vertices.length / 3, faces: scene.faces.length / 3, diag: Number(diag.toFixed(6)) },
  detection: {
    ...detection,
    angleErrorDeg: Number(angleOf(det.plane).toFixed(6)),
    rmsNorm: Number(det.rmsNorm.toFixed(9)),
    baseline: {
      ms: Number(baselineMs.toFixed(3)),
      speedup: Number((baselineMs / detection.median).toFixed(2)),
      angleErrorDeg: Number(angleOf(baseline.plane).toFixed(3)),
      grid: '36 fibonacci normals x 7 offsets, 1/6 vertex subsample',
    },
  },
  certificate: {
    grid: stats(certTimes),
    brute120SamplesMs: Number(oracleMs.toFixed(3)),
    speedup: Number((oracleMs / stats(certTimes).median).toFixed(2)),
    gridMaxNorm: Number(cert.maxNorm.toFixed(9)),
    oracleMaxNorm: Number((oracleMax / diag).toFixed(9)),
    samples: cert.samples,
    tolRel: cert.tolRel,
    ok: cert.ok,
    fullPassSamples: 5954,
  },
  repair: { iterations: 2, ms: Number(repairMs.toFixed(3)) },
};
fs.writeFileSync(path.join(here, 'result.json'), JSON.stringify(out, null, 2) + '\n');
console.log(JSON.stringify(out, null, 2));
