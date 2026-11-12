// Assertions over the real sample data. Exit 0 only if every check passes.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  detectSymmetry, certifySymmetry, symmetrize, surfaceDeviation,
  makeTriangleIndex, closestPointOnMesh, reflectPoint, surfaceSamples, bounds,
} from '../src/core.mjs';
import { buildSymmetric, transformMesh } from '../tools/meshgen.mjs';
import { bruteSurfaceDistance } from './oracle.mjs';

const here = path.dirname(fileURLToPath(import.meta.url));
const data = (f) => JSON.parse(fs.readFileSync(path.join(here, '..', 'data', f), 'utf8'));
const scene = data('scene.json');
const exact = data('exact.json');
const asym = data('asym.json');

let pass = 0, fail = 0;
const ok = (name, cond, detail) => {
  if (cond) { pass++; console.log(`PASS  ${name}  ${detail}`); }
  else { fail++; console.log(`FAIL  ${name}  ${detail}`); }
};

const diag = bounds(scene.vertices).diag;
const det = detectSymmetry(scene);
const dot = Math.abs(det.plane.n[0] * scene.truth.n[0] + det.plane.n[1] * scene.truth.n[1] + det.plane.n[2] * scene.truth.n[2]);
const angle = Math.acos(Math.min(1, dot)) * 180 / Math.PI;
ok('G1/S1 perturbed plane angle', angle < 0.5, `angle=${angle.toFixed(6)} deg, rms=${(det.rmsNorm * 100).toFixed(4)}% diag`);

const exactDet = detectSymmetry(exact);
const exactAngle = Math.acos(Math.min(1,
  Math.abs(exactDet.plane.n[0] * exact.truth.n[0] + exactDet.plane.n[1] * exact.truth.n[1] + exactDet.plane.n[2] * exact.truth.n[2]))) * 180 / Math.PI;
ok('G2/S2 exact mesh plane', exactDet.symmetric && exactAngle < 0.01, `angle=${exactAngle.toExponential(3)} deg, rms=${exactDet.rmsNorm.toExponential(3)}`);

const certBefore = certifySymmetry(scene, det.plane, { tolRel: 0.005 });
ok('G3/S3 asymmetry detected', certBefore.ok === false && certBefore.maxNorm > 0.005,
  `ok=${certBefore.ok}, max=${(certBefore.maxNorm * 100).toFixed(4)}% diag, over=${certBefore.over}/${certBefore.samples}`);

const repaired = symmetrize(scene, det.plane, 6);
const certAfter = certifySymmetry({ vertices: repaired.vertices, faces: repaired.faces }, det.plane, { tolRel: 0.005 });
ok('G4/S4 repair passes certificate', certAfter.ok === true && certAfter.maxNorm < certBefore.maxNorm / 4,
  `max ${(certBefore.maxNorm * 100).toFixed(4)}% -> ${(certAfter.maxNorm * 100).toFixed(4)}% diag (${(certBefore.maxNorm / certAfter.maxNorm).toFixed(2)}x)`);

const asymDet = detectSymmetry(asym);
ok('G5/S5 abstain on asymmetric mesh', asymDet.symmetric === false && asymDet.rmsNorm > 0.01,
  `symmetric=${asymDet.symmetric}, rms=${(asymDet.rmsNorm * 100).toFixed(3)}% diag (threshold 1.0%)`);

// Scale invariance: same shape at 1e-6 .. 1e6 must give the same normalized verdict.
const base = buildSymmetric({ perturb: 61, amplitude: 0.03, seed: 777 });
const scales = [1e-6, 1e-3, 1, 1e3, 1e6];
const rows = [];
for (const s of scales) {
  const m = transformMesh(base, { euler: [0.7, -0.4, 1.1], scale: s, translate: [1.5 * s, -0.7 * s, 2.2 * s] });
  const d = detectSymmetry(m);
  const tr = { n: [1.5 * 0 + Math.cos(0.7) * Math.cos(-0.4), 0, 0], d: 0 };
  // truth normal = R * e_x, computed the same way the generator does
  const t = transformMesh({ vertices: [1, 0, 0], faces: [] }, { euler: [0.7, -0.4, 1.1] });
  const tn = t.vertices;
  const ang = Math.acos(Math.min(1, Math.abs(d.plane.n[0] * tn[0] + d.plane.n[1] * tn[1] + d.plane.n[2] * tn[2]))) * 180 / Math.PI;
  rows.push({ s, rms: d.rmsNorm, ang, sym: d.symmetric, off: d.plane.d / bounds(m.vertices).diag });
}
const rmsSpread = (Math.max(...rows.map((r) => r.rms)) - Math.min(...rows.map((r) => r.rms))) / rows[0].rms;
const angleSpread = Math.max(...rows.map((r) => r.ang)) - Math.min(...rows.map((r) => r.ang));
ok('G6/S6 scale invariance 1e-6..1e6', rows.every((r) => r.sym) && rmsSpread < 0.01 && angleSpread < 0.01,
  `rmsSpread=${(rmsSpread * 100).toFixed(4)}%, angleSpread=${angleSpread.toExponential(2)} deg, verdicts=${rows.map((r) => (r.sym ? 'T' : 'F')).join('')}`);

// Independent oracle: compare distance kernel against bench/oracle.mjs on real samples.
const samples = surfaceSamples(scene);
const idx = makeTriangleIndex(scene.vertices, scene.faces, diag / 32);
let maxDiff = 0, checked = 0;
for (let i = 0; i < samples.length; i += 3 * 37) {
  const q = reflectPoint(det.plane, [samples[i], samples[i + 1], samples[i + 2]]);
  const a = closestPointOnMesh(idx, q[0], q[1], q[2]);
  const b = bruteSurfaceDistance(scene, q[0], q[1], q[2]);
  maxDiff = Math.max(maxDiff, Math.abs(a - b));
  checked++;
}
ok('G7/S7 oracle kernel agreement', checked >= 80 && maxDiff < 1e-9, `checked=${checked} samples, max|grid-oracle|=${maxDiff.toExponential(3)}`);

// Determinism of reported aggregates across two calls.
const again = surfaceDeviation(scene, det.plane);
const ref = JSON.parse(fs.readFileSync(path.join(here, '..', 'data', 'ref.json'), 'utf8'));
ok('G8/S8 baked reference matches live core', Math.abs(again.maxNorm - ref.certificate.maxNorm) === 0 &&
  Math.abs(det.rmsNorm - ref.detectionRmsNorm) === 0,
  `maxNorm live=${again.maxNorm}, baked=${ref.certificate.maxNorm}; rms live=${det.rmsNorm}, baked=${ref.detectionRmsNorm}`);

console.log(`\n${pass} passed / ${fail} failed (exit ${fail ? 1 : 0})`);
process.exit(fail ? 1 : 0);
