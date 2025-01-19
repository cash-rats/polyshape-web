import { useState, useEffect, useRef } from "react";
import { Dimensions } from "~/types";

export const useTrianglify = (defaultColorPalette: string[]) => {
  const [dimensions, setDimensions] = useState<Dimensions>({
    width: 900,
    height: 900,
    cellSize: 75, // Adjust for triangle size (smaller = more detailed)
    variance: 0.75, // Adjust for triangle randomness (0-1)
    patternIntensity: 0.5, // Adjust bias for gradient intensity
    xColors: defaultColorPalette,
    yColors: defaultColorPalette,
  });

  const patternContainerRef = useRef<HTMLDivElement>(null);
  const [triangifyPattern, setTriangifyPattern] = useState<any>(null);

  // Update pattern when dimensions change
  useEffect(() => {
    const updatePattern = () => {
      if (typeof window !== 'undefined' && window.trianglify && patternContainerRef.current) {
        // Get the actual container width for the pattern
        const pattern = window.trianglify({
          variance: dimensions.variance,
          cellSize: dimensions.cellSize,
          colorFunction: window.trianglify.colorFunctions.interpolateLinear(dimensions.patternIntensity),
          xColors: dimensions.xColors,
          yColors: dimensions.yColors,
          width: dimensions.width,
          height: dimensions.height,
        });
        setTriangifyPattern(pattern);
        const patternCanvas = pattern.toCanvas();

        patternCanvas.style.width = 'auto';
        patternCanvas.style.height = 'auto';

        patternCanvas.style.maxWidth = '100%';
        patternCanvas.style.maxHeight = '100%';

        // Set canvas dimensions explicitly
        patternCanvas.style.aspectRatio = `auto ${dimensions.width} / ${dimensions.height}`;
        patternCanvas.style.objectFit = 'contain';
        patternCanvas.style.overflowClipMargin = 'content-box';
        patternCanvas.style.overflow = 'clip';

        patternContainerRef.current.innerHTML = '';
        patternContainerRef.current.appendChild(patternCanvas);
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
    patternContainerRef,
    triangifyPattern,
    dimensions,
    setWidth,
    setHeight,
    setCellSize,
    setShapeVariance,
    setPatternIntensity,
    setColorPalette,
  }
}
