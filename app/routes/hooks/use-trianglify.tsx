import { useState, useEffect, useRef } from "react";

export const useTrianglify = () => {
  const [dimensions, setDimensions] = useState({
    width: 1440,
    height: 900,
    cellSize: 75
  });

  const patternRef = useRef<HTMLDivElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);

  // Update pattern when dimensions change
  useEffect(() => {
    const updatePattern = () => {
      if (typeof window !== 'undefined' && window.trianglify && patternRef.current && containerRef.current) {
        // Get the actual container width for the pattern
        const containerWidth = containerRef.current.clientWidth;

        const pattern = window.trianglify({
          width: containerWidth,
          height: dimensions.height,
          cellSize: dimensions.cellSize
        });

        patternRef.current.innerHTML = '';
        patternRef.current.appendChild(pattern.toSVG());
      }
    };

    updatePattern();
  }, [dimensions]);

  const setWidth = (width: number) => {
    setDimensions({ ...dimensions, width });
  }

  const setHeight = (height: number) => {
    setDimensions({ ...dimensions, height });
  }

  const setCellSize = (cellSize: number) => {
    setDimensions({ ...dimensions, cellSize });
  }

  return {
    patternRef,
    containerRef,
    dimensions,
    setWidth,
    setHeight,
    setCellSize
  }
}