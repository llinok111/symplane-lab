// Deterministic sample-data generator. Run: node tools/gen-sample.mjs
// Writes data/scene.json, data/exact.json, data/asym.json and data/ref.json.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { buildSymmetric, buildAsymmetric, transformMesh, truthPlane } from './meshgen.mjs';
import { detectSymmetry, certifySymmetry, surfaceDeviation, symmetrize } from '../src/core.mjs';

const here = path.dirname(fileURLToPath(import.meta.url));
const outDir = path.resolve(here, '..', 'data');

const MAIN_TRANSFORM = { euler: [0.7, -0.4, 1.1], scale: 2.4, translate: [1.5, -0.7, 2.2] };

function pack(name, mesh, transform, extra = {}) {
  const transformed = transformMesh(mesh, transform);
  return {
    name,
    vertices: transformed.vertices.map((x) => Number(x.toFixed(12))),
    faces: transformed.faces,
    truth: truthPlane(transform),
    ...extra,
  };
}

const symBase = buildSymmetric({ perturb: 61, amplitude: 0.03, seed: 777 });
const exactBase = buildSymmetric({ perturb: 0 });
const asymBase = buildAsymmetric();

const scene = pack('sym-rotated-scaled', symBase, MAIN_TRANSFORM, { perturbedVertices: symBase.perturbedVertices });
const exact = pack('sym-exact', exactBase, { euler: [-0.3, 0.9, 0.2], scale: 1, translate: [0, 0, 0] });
const asym = pack('asymmetric', asymBase, { euler: [0.2, 0.1, -0.5], scale: 1, translate: [0, 0, 0] });

function summarize(mesh, truth) {
  const det = detectSymmetry(mesh);
  const cert = certifySymmetry(mesh, det.plane, { tolRel: 0.005 });
  const dot = det.plane.n[0] * truth.n[0] + det.plane.n[1] * truth.n[1] + det.plane.n[2] * truth.n[2];
  const angle = Math.acos(Math.min(1, Math.abs(dot))) * 180 / Math.PI;
  const signedTruthD = (dot >= 0 ? 1 : -1) * truth.d;
  const offErr = Math.abs(det.plane.d - signedTruthD);
  return { detect: det, cert, angleDeg: angle, offsetError: offErr };
}

const mainDet = summarize(scene, scene.truth);
const exactDet = summarize(exact, exact.truth);
const asymDet = detectSymmetry(asym);

const liveStride = 5;
const liveCount = Math.floor(scene.vertices.length / 3 / liveStride);
const liveSamples = new Float64Array(liveCount * 3);
for (let n = 0; n < liveCount; n++) {
  const i = n * liveStride * 3;
  liveSamples[n * 3] = scene.vertices[i];
  liveSamples[n * 3 + 1] = scene.vertices[i + 1];
  liveSamples[n * 3 + 2] = scene.vertices[i + 2];
}
const liveDev = surfaceDeviation(scene, mainDet.detect.plane, { samples: liveSamples, includeCentroids: false });
const repaired = symmetrize(scene, mainDet.detect.plane, 6);
const repairedMesh = { vertices: repaired.vertices, faces: repaired.faces };
const certBefore = certifySymmetry(scene, mainDet.detect.plane, { tolRel: 0.005 });
const certAfter = certifySymmetry(repairedMesh, mainDet.detect.plane, { tolRel: 0.005 });
const devBefore = surfaceDeviation(scene, mainDet.detect.plane);
const devAfter = surfaceDeviation(repairedMesh, mainDet.detect.plane);

fs.mkdirSync(outDir, { recursive: true });
fs.writeFileSync(path.join(outDir, 'scene.json'), JSON.stringify(scene) + '\n');
fs.writeFileSync(path.join(outDir, 'exact.json'), JSON.stringify(exact) + '\n');
fs.writeFileSync(path.join(outDir, 'asym.json'), JSON.stringify(asym) + '\n');
const ref = {
  generator: 'tools/gen-sample.mjs',
  detectedPlane: mainDet.detect.plane,
  detectionRmsNorm: mainDet.detect.rmsNorm,
  detectionAngleErrorDeg: mainDet.angleDeg,
  detectionOffsetErrorRelDiag: mainDet.offsetError / mainDet.detect.diag,
  certificate: { ok: certBefore.ok, maxNorm: certBefore.maxNorm, tolRel: certBefore.tolRel, over: certBefore.over, samples: certBefore.samples },
  repaired: {
    iterations: 6,
    certOk: certAfter.ok,
    maxNormBefore: devBefore.maxNorm,
    maxNormAfter: devAfter.maxNorm,
    rmsNormBefore: devBefore.rmsNorm,
    rmsNormAfter: devAfter.rmsNorm,
  },
  exact: { angleDeg: exactDet.angleDeg, offsetErrorRelDiag: exactDet.offsetError / exactDet.detect.diag, rmsNorm: exactDet.detect.rmsNorm, symmetric: exactDet.detect.symmetric },
  asymmetric: { symmetric: asymDet.symmetric, rmsNorm: asymDet.rmsNorm },
  liveCheck: { stride: liveStride, samples: liveCount, maxNormAtDetected: liveDev.maxNorm },
  vertexCount: scene.vertices.length / 3,
  faceCount: scene.faces.length / 3,
};
fs.writeFileSync(path.join(outDir, 'ref.json'), JSON.stringify(ref, null, 2) + '\n');
console.log(JSON.stringify({
  files: ['scene.json', 'exact.json', 'asym.json', 'ref.json'],
  detAngleDeg: mainDet.angleDeg,
  detRmsNorm: mainDet.detect.rmsNorm,
  certOk: certBefore.ok,
  certMaxNorm: certBefore.maxNorm,
  repairedCertOk: certAfter.ok,
  repairedMaxNorm: devAfter.maxNorm,
  asymSymmetric: asymDet.symmetric,
  asymRms: asymDet.rmsNorm,
}, null, 2));
