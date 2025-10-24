import { useEffect, useRef } from 'react';
import * as THREE from 'three';

const colorPalette = [0x35f0ff, 0xff3df0, 0x00ffc6, 0xffb347, 0x9f7aea];

const mapCandidateColor = (candidateId = '') => {
  const index = Math.abs(
    candidateId.split('').reduce((acc, char) => acc + char.charCodeAt(0), 0)
  ) % colorPalette.length;
  return new THREE.Color(colorPalette[index]);
};

const createThreeScene = (container, candidateId, playingRef) => {
  const width = container.clientWidth;
  const height = container.clientHeight;

  const scene = new THREE.Scene();
  scene.background = new THREE.Color('#050912');

  const camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 100);
  camera.position.set(0, 0, 8);

  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setSize(width, height);
  renderer.outputEncoding = THREE.sRGBEncoding;
  container.appendChild(renderer.domElement);

  const ambient = new THREE.AmbientLight(0x66ffff, 0.5);
  scene.add(ambient);

  const point = new THREE.PointLight(0xff3df0, 1.2, 15);
  point.position.set(4, 3, 5);
  scene.add(point);

  const targetGeo = new THREE.TorusKnotGeometry(1.5, 0.35, 128, 24);
  const targetMat = new THREE.MeshStandardMaterial({
    color: 0x142848,
    emissive: 0x06122a,
    metalness: 0.6,
    roughness: 0.25,
  });
  const targetMesh = new THREE.Mesh(targetGeo, targetMat);
  scene.add(targetMesh);

  const ligandColor = mapCandidateColor(candidateId);
  const ligandGeo = new THREE.SphereGeometry(0.5, 48, 32);
  const ligandMat = new THREE.MeshStandardMaterial({
    color: ligandColor,
    emissive: ligandColor.clone().multiplyScalar(0.6),
    roughness: 0.2,
    metalness: 0.7,
  });
  const ligand = new THREE.Mesh(ligandGeo, ligandMat);
  ligand.position.set(2.3, 0, 0);
  scene.add(ligand);

  const orbitLight = new THREE.PointLight(ligandColor.getHex(), 1, 7);
  orbitLight.position.copy(ligand.position);
  scene.add(orbitLight);

  const clock = new THREE.Clock();
  const state = { animationFrame: null };

  const animate = () => {
    const elapsed = clock.getElapsedTime();
    const wobble = playingRef.current ? 0.35 : 0.1;
    targetMesh.rotation.x += 0.0025;
    targetMesh.rotation.y += 0.0035;

    ligand.position.set(
      2.3 + Math.sin(elapsed * 1.8) * wobble,
      Math.cos(elapsed * 1.5) * wobble * 0.6,
      Math.sin(elapsed * 1.3) * wobble * 0.5
    );
    ligand.rotation.y += 0.01;
    orbitLight.position.copy(ligand.position);

    renderer.render(scene, camera);
    state.animationFrame = requestAnimationFrame(animate);
  };

  animate();

  const handleResize = () => {
    const w = container.clientWidth;
    const h = container.clientHeight;
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    renderer.setSize(w, h);
  };

  window.addEventListener('resize', handleResize);

  return {
    dispose() {
      cancelAnimationFrame(state.animationFrame);
      window.removeEventListener('resize', handleResize);
      targetGeo.dispose();
      targetMat.dispose();
      ligandGeo.dispose();
      ligandMat.dispose();
      renderer.dispose();
      container.removeChild(renderer.domElement);
    },
    updateCandidate(nextId) {
      const color = mapCandidateColor(nextId);
      ligand.material.color.copy(color);
      ligand.material.emissive.copy(color.clone().multiplyScalar(0.6));
      orbitLight.color.copy(color);
    },
  };
};

const NGLViewer = ({ flowId, candidateId, playing }) => {
  const containerRef = useRef(null);
  const runtimeRef = useRef(null);
  const playingRef = useRef(playing);

  useEffect(() => {
    playingRef.current = playing;
  }, [playing]);

  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    // In offline mode fall back to Three.js renderer instantly
    runtimeRef.current = createThreeScene(container, candidateId, playingRef);

    return () => {
      runtimeRef.current?.dispose?.();
      runtimeRef.current = null;
    };
  }, []);

  useEffect(() => {
    runtimeRef.current?.updateCandidate?.(candidateId);
  }, [candidateId, flowId]);

  return (
    <div className="relative h-full w-full overflow-hidden rounded-2xl border border-cyan-700/50 bg-gradient-to-br from-[#080d1b] to-[#0d172e]">
      <div ref={containerRef} className="h-full w-full" />
      <div className="pointer-events-none absolute left-3 top-3 rounded-full border border-cyan-500/50 bg-cyan-500/20 px-3 py-1 text-[10px] uppercase tracking-[0.2em] text-cyan-100/80">
        {flowId.replace(/_/g, ' ')}
      </div>
    </div>
  );
};

export default NGLViewer;
