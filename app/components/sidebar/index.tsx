import ColorTiles from "../color-tiles";

import { Label } from "~/components/ui/label";
import { Slider } from "~/components/ui/slider";
import { Dimensions } from "~/types";

interface SidebarProps {
  dimensions: Dimensions;
  colorBrewsers: Record<string, chroma.Scale>;
  onChangeWidth: (width: number) => void;
  onChangeHeight: (height: number) => void;
  onChangeCellSize: (cellSize: number) => void;
  onChangePatternIntensity: (patternIntensity: number) => void;
  onChangeShapeVariance: (shapeVariance: number) => void;
  onChangeColorPalette: (colorPalette: string[]) => void;
}

export default function Sidebar({
  dimensions,
  colorBrewsers,
  onChangeWidth,
  onChangeHeight,
  onChangeCellSize,
  onChangePatternIntensity,
  onChangeShapeVariance,
  onChangeColorPalette,
}: SidebarProps) {
  return (
      <div className="w-[400px] flex-shrink-0 border-r bg-background/95 backdrop-blur supports-[backdrop-filter]:bg-background/60 h-screen overflow-y-auto">
        <div className="flex flex-col gap-8 p-6">
          {/* Width Input */}
          <div className="space-y-2">
            <div className="flex items-center justify-start gap-4">
              <div>
                <label htmlFor="width" className="block font-bold">Width</label>
                <input
                  id="width"
                  type="number"
                  value={dimensions.width}
                  onChange={(e) => onChangeWidth(+e.target.value)}
                  className="w-24 rounded border px-2 py-1"
                />
              </div>
              <div>
                <label htmlFor="height" className="block font-bold">Height</label>
                <input
                  id="height"
                  type="number"
                  value={dimensions.height}
                  onChange={(e) => onChangeHeight(+e.target.value)}
                  className="w-24 rounded border px-2 py-1"
                />
              </div>
            </div>
          </div>

          {/* Cell Size Slider */}
          <div className="space-y-4">
            <div className="flex items-center justify-between">
              <Label htmlFor="cellSize">CELL SIZE</Label>
              <span className="w-12 rounded-md border border-transparent px-2 py-0.5 text-right text-sm text-muted-foreground hover:border-border">
                {dimensions.cellSize}px
              </span>
            </div>
            <Slider
              id="cellSize"
              min={5}
              max={100}
              step={5}
              value={[dimensions.cellSize]}
              onValueChange={(value) => onChangeCellSize(value[0])}
              className="[&_[role=slider]]:h-4 [&_[role=slider]]:w-4"
            />
          </div>

          {/* Pattern intensity slider */}
          <div className="space-y-4">
            <div className="flex items-center justify-between">
              <Label htmlFor="patternIntensity">PATTERN INTENSITY</Label>
            </div>
            <Slider
              id="patternIntensity"
              min={0}
              max={1}
              step={0.1}
              value={[dimensions.patternIntensity]}
              onValueChange={(value) => onChangePatternIntensity(value[0])}
              className="[&_[role=slider]]:h-4 [&_[role=slider]]:w-4"
            />
          </div>

          {/* Shape variance slider */}
          <div className="space-y-4">
            <div className="flex items-center justify-between">
              <Label htmlFor="shapeVariance">SHAPE VARIANCE</Label>
            </div>
            <Slider
              id="shapeVariance"
              min={0}
              max={1}
              step={0.1}
              value={[dimensions.variance]}
              onValueChange={(value) => onChangeShapeVariance(value[0])}
              className="[&_[role=slider]]:h-4 [&_[role=slider]]:w-4"
            />
          </div>

          {/* Color Palette Selector*/}
          <div className="space-y-4 flex flex-col items-center">
            <Label htmlFor="colorPalette">COLOR PALETTE</Label>
            {
              Object.keys(colorBrewsers).map((colorDomain) => {
                const colors = colorBrewsers[colorDomain].colors();
                return (<ColorTiles key={colorDomain} colors={colors} onSelect={onChangeColorPalette} />)
              })
            }
          </div>
        </div>
      </div>
  );
}