// Rendering + interaction only; every number shown is computed by src/core.mjs.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { surfaceDeviation, bounds, normalizePlane } from './src/core.mjs';

const fail = (e) => { document.getElementById('err').textContent = String(e && e.stack || e); };
try {
  const [sceneData, ref] = await Promise.all([
    fetch('./data/scene.json').then((r) => r.json()),
    fetch('./data/ref.json').then((r) => r.json()),
  ]);
  const diag = bounds(sceneData.vertices).diag;
  const center = bounds(sceneData.vertices).center;
  const base = normalizePlane(ref.detectedPlane);

  // Same deterministic vertex subset the generator baked into ref.liveCheck.
  const stride = ref.liveCheck.stride;
  const count = Math.floor(sceneData.vertices.length / 3 / stride);
  const samples = new Float64Array(count * 3);
  for (let n = 0; n < count; n++) {
    const i = n * stride * 3;
    samples[n * 3] = sceneData.vertices[i];
    samples[n * 3 + 1] = sceneData.vertices[i + 1];
    samples[n * 3 + 2] = sceneData.vertices[i + 2];
  }

  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  renderer.setSize(innerWidth, innerHeight);
  document.body.appendChild(renderer.domElement);

  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x0d1117);
  const camera = new THREE.PerspectiveCamera(45, innerWidth / innerHeight, diag / 1000, diag * 100);
  camera.position.set(center[0] + diag * 0.9, center[1] + diag * 0.6, center[2] + diag * 1.1);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.target.set(...center);
  controls.update();

  scene.add(new THREE.HemisphereLight(0xc9d1d9, 0x30363d, 1.1));
  const key = new THREE.DirectionalLight(0xffffff, 1.6);
  key.position.set(center[0] + diag, center[1] + diag, center[2] + diag);
  scene.add(key);

  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(sceneData.vertices, 3));
  geo.setIndex(sceneData.faces);
  geo.computeVertexNormals();
  scene.add(new THREE.Mesh(geo, new THREE.MeshStandardMaterial({
    color: 0x58a6ff, metalness: 0.15, roughness: 0.55, side: THREE.DoubleSide, flatShading: false,
  })));

  // Mirror plane visual: a thin disc aligned with the current plane.
  const planeMesh = new THREE.Mesh(
    new THREE.PlaneGeometry(diag * 1.15, diag * 1.15),
    new THREE.MeshBasicMaterial({ color: 0x3fb950, transparent: true, opacity: 0.10, side: THREE.DoubleSide, depthWrite: false }),
  );
  scene.add(planeMesh);
  const grid = new THREE.GridHelper(diag * 1.15, 12, 0x3fb950, 0x238636);
  grid.rotation.x = Math.PI / 2;
  planeMesh.add(grid);

  const ui = {
    a: document.getElementById('a'), b: document.getElementById('b'), c: document.getElementById('c'),
    live: document.getElementById('live'), baked: document.getElementById('baked'),
    gate: document.getElementById('gate'), match: document.getElementById('match'),
  };

  // Orthonormal basis of the detected plane (same construction as core.mjs).
  const n0 = base.n;
  const t = Math.abs(n0[0]) < Math.abs(n0[1]) ? (Math.abs(n0[0]) < Math.abs(n0[2]) ? [1, 0, 0] : [0, 0, 1]) : (Math.abs(n0[1]) < Math.abs(n0[2]) ? [0, 1, 0] : [0, 0, 1]);
  const uRaw = [n0[1] * t[2] - n0[2] * t[1], n0[2] * t[0] - n0[0] * t[2], n0[0] * t[1] - n0[1] * t[0]];
  const ul = Math.hypot(...uRaw);
  const u = uRaw.map((x) => x / ul);
  const w = [n0[1] * u[2] - n0[2] * u[1], n0[2] * u[0] - n0[0] * u[2], n0[0] * u[1] - n0[1] * u[0]];

  const planeOf = () => {
    const a = Number(ui.a.value), b = Number(ui.b.value), c = Number(ui.c.value);
    const n = normalizePlane({ n: [n0[0] + a * u[0] + b * w[0], n0[1] + a * u[1] + b * w[1], n0[2] + a * u[2] + b * w[2]], d: 1 }).n;
    return { n, d: base.d + c * diag };
  };

  const apply = () => {
    const plane = planeOf();
    const dev = surfaceDeviation({ vertices: sceneData.vertices, faces: sceneData.faces }, plane, { samples, includeCentroids: false });
    ui.live.textContent = (dev.maxNorm * 100).toFixed(4) + '%';
    ui.baked.textContent = (ref.liveCheck.maxNormAtDetected * 100).toFixed(4) + '%';
    const pass = dev.maxNorm <= 0.005;
    ui.gate.textContent = pass ? 'PASS' : 'FAIL';
    ui.gate.className = pass ? 'ok' : 'bad';
    ui.match.textContent = Math.abs(dev.maxNorm - ref.liveCheck.maxNormAtDetected) < 1e-12 ? '(matches node bit-for-bit)' : '';
    planeMesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), new THREE.Vector3(...plane.n).normalize());
    planeMesh.position.copy(new THREE.Vector3(...plane.n).multiplyScalar(plane.d));
    renderer.render(scene, camera);
  };

  for (const el of [ui.a, ui.b, ui.c]) el.addEventListener('input', apply);
  document.getElementById('reset').addEventListener('click', () => { ui.a.value = 0; ui.b.value = 0; ui.c.value = 0; apply(); });
  addEventListener('resize', () => {
    camera.aspect = innerWidth / innerHeight; camera.updateProjectionMatrix(); renderer.setSize(innerWidth, innerHeight); apply();
  });
  document.getElementById('intro').textContent =
    `mesh ${sceneData.vertices.length / 3} verts / ${sceneData.faces.length / 3} faces · detected plane error ${ref.detectionAngleErrorDeg.toFixed(4)}°`;
  apply();
  const tick = () => { controls.update(); renderer.render(scene, camera); requestAnimationFrame(tick); };
  tick();
} catch (e) { fail(e); }
