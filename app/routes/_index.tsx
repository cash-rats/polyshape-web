import { useTrianglify } from "./hooks/use-trianglify";
import { useColorBrewer } from "./hooks/use-color-brewer";

import Sidebar from "~/components/sidebar";

// TODOs:
//  - [x] Add cell size slider
//  - [x] Add pattern intensity slider
//  - [x] Add pattern triangle variance slider
//  - [x] Fixed color
//  - [x] Selectable color
//  - [ ] pattern canvas should be set to aspect ratio
export default function Index() {
  const { colorBrewsers, defaultColorPalette } = useColorBrewer();

  const {
    patternRef,
    containerRef,
    dimensions,
    setWidth,
    setHeight,
    setCellSize,
    setPatternIntensity,
    setShapeVariance,
    setColorPalette,
  } = useTrianglify(defaultColorPalette);

  return (
    <div className="min-h-screen flex">
      {/* Fixed-width sidebar */}
      <Sidebar
        dimensions={dimensions}
        colorBrewsers={colorBrewsers}
        onChangeWidth={setWidth}
        onChangeHeight={setHeight}
        onChangeCellSize={setCellSize}
        onChangePatternIntensity={setPatternIntensity}
        onChangeShapeVariance={setShapeVariance}
        onChangeColorPalette={setColorPalette}
      />

      {/* Flexible-width main content */}
      <div
        ref={containerRef}
        className="flex-grow bg-background p-6 pt-[40px] px-[60px] pb-[70px]"
      >
        <div
          ref={patternRef}
          className="w-full h-full flex items-center justify-center"
        />
      </div>
    </div>
  );
}
