import { useEffect, useRef } from 'react';
import * as THREE from 'three';

const Starfield = () => {
  const containerRef = useRef(null);
  const runtimeRef = useRef(null);

  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    const width = container.clientWidth;
    const height = container.clientHeight;

    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(60, width / height, 1, 1000);
    camera.position.z = 1;

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(window.devicePixelRatio || 1);
    container.appendChild(renderer.domElement);

    const starsGeometry = new THREE.BufferGeometry();
    const starCount = 600;
    const positions = [];
    for (let i = 0; i < starCount; i += 1) {
      const x = THREE.MathUtils.randFloatSpread(600);
      const y = THREE.MathUtils.randFloatSpread(600);
      const z = THREE.MathUtils.randFloatSpread(600);
      positions.push(x, y, z);
    }
    starsGeometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));

    const starsMaterial = new THREE.PointsMaterial({
      color: 0x35f0ff,
      size: 1.2,
      transparent: true,
      opacity: 0.25,
      blending: THREE.AdditiveBlending,
    });

    const stars = new THREE.Points(starsGeometry, starsMaterial);
    scene.add(stars);

    const state = { animationFrame: null };

    const animate = () => {
      stars.rotation.x += 0.0008;
      stars.rotation.y += 0.0005;
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

    runtimeRef.current = {
      dispose() {
        cancelAnimationFrame(state.animationFrame);
        window.removeEventListener('resize', handleResize);
        starsGeometry.dispose();
        starsMaterial.dispose();
        renderer.dispose();
        container.removeChild(renderer.domElement);
      },
    };

    return () => {
      runtimeRef.current?.dispose?.();
      runtimeRef.current = null;
    };
  }, []);

  return <div ref={containerRef} className="pointer-events-none absolute inset-0 -z-10" />;
};

export default Starfield;
