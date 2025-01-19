import { useTrianglify } from "./hooks/use-trianglify";
import { useColorBrewer } from "./hooks/use-color-brewer";

import Sidebar from "~/components/sidebar";
import ExportDialog from "~/components/export-dialog";

// TODOs:
//  - [x] Add cell size slider
//  - [x] Add pattern intensity slider
//  - [x] Add pattern triangle variance slider
//  - [x] Fixed color
//  - [x] Selectable color
//  - [x] pattern canvas should be set to aspect ratio
export default function Index() {
  const { colorBrewsers, defaultColorPalette } = useColorBrewer();

  const {
    patternContainerRef,
    triangifyPattern,
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
      <div className="flex justify-center items-center flex-grow bg-[#e8e8e8] p-6 pt-[40px] px-[60px] pb-[70px] overflow-y-hidden relative flex-auto z-[1]">
        <div
          ref={patternContainerRef}
          className="flex items-center justify-center h-[765px] w-[765px]"
        />

        {/* Button to export the pattern (centered at bottom) */}
        <div className="absolute bottom-4 left-1/2 transform -translate-x-1/2">
          <ExportDialog triangifyPattern={triangifyPattern} dimensions={dimensions} />
        </div>
      </div>
    </div>
  );
}
