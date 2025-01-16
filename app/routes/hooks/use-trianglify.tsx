import { useState, useEffect, useRef } from "react";


export const useTrianglify = (defaultColorPalette: string[]) => {
  const [dimensions, setDimensions] = useState({
    width: 1440,
    height: 900,
    cellSize: 75, // Adjust for triangle size (smaller = more detailed)
    variance: 0.75, // Adjust for triangle randomness (0-1)
    patternIntensity: 0.5, // Adjust bias for gradient intensity
    xColors: defaultColorPalette,
    yColors: defaultColorPalette,
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
          variance: dimensions.variance,
          cellSize: dimensions.cellSize,
          colorFunction: window.trianglify.colorFunctions.interpolateLinear(dimensions.patternIntensity),
          xColors: dimensions.xColors,
          yColors: dimensions.yColors,
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

  const setShapeVariance = (variance: number) => {
    setDimensions({ ...dimensions, variance });
  }

  const setPatternIntensity = (patternIntensity: number) => {
    setDimensions({ ...dimensions, patternIntensity });
  }

  const setCellSize = (cellSize: number) => {
    setDimensions({ ...dimensions, cellSize });
  }

  const setColorPalette = (colorPalette: string[]) => {
    setDimensions({ ...dimensions, xColors: colorPalette, yColors: colorPalette });
  }

  return {
    patternRef,
    containerRef,
    dimensions,
    setWidth,
    setHeight,
    setCellSize,
    setShapeVariance,
    setPatternIntensity,
    setColorPalette,
  }
}
