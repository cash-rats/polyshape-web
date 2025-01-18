import { DownloadIcon } from '@radix-ui/react-icons';
import { useTrianglify } from "./hooks/use-trianglify";
import { useColorBrewer } from "./hooks/use-color-brewer";

import Sidebar from "~/components/sidebar";

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
        className="flex-grow bg-[#e8e8e8] p-6 pt-[40px] px-[60px] pb-[70px] overflow-y-hidden relative"
      >
        <div
          ref={patternRef}
          className="w-full h-full flex items-center justify-center"
        />

        {/* Button to export the pattern (centered at bottom) */}
        <div className="absolute bottom-4 left-1/2 transform -translate-x-1/2">
          <button
            type="button"
            className="flex items-center px-4 py-2 bg-white text-[#16C47F] rounded-md hover:bg-gray-100"
            onClick={() => {
              // TODO: Add your export logic here
            }}
          >
            <DownloadIcon className="w-5 h-5 stroke-current text-[#16C47F]" />
            <span className="ml-2">Export</span>
          </button>
        </div>
      </div>
    </div>
  );
}
