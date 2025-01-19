import { useEffect, useState } from 'react';
import { Dimensions } from "~/types";

export const useTrianglifyPreview = (trianglifyPattern: any, dimensions: Dimensions) => {
  const [previewRef, setPreviewRef] = useState<HTMLDivElement | null>(null);

  useEffect(() => {
    if (previewRef && trianglifyPattern) {
      const previewCanvas = trianglifyPattern.toCanvas();
      previewCanvas.style.width = dimensions.width;
      previewCanvas.style.height = dimensions.height;

      previewCanvas.style.maxWidth = '100%';
      previewCanvas.style.maxHeight = '100%';

      previewCanvas.style.aspectRatio = `${dimensions.width} / ${dimensions.height}`;
      previewCanvas.style.objectFit = 'contain';

      previewRef.innerHTML = '';
      previewRef.appendChild(previewCanvas);
    }
  }, [trianglifyPattern, previewRef, dimensions]);

  return {
    previewRef: setPreviewRef,
  };
};
