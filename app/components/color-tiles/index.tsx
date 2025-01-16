interface ColorTilesProps {
  colors: string[];
  onSelect: (colors: string[]) => void
}

export default function ColorTiles({ colors, onSelect }: ColorTilesProps) {
  return (
    <button
      onClick={() => onSelect(colors)}
      className="w-full flex flex-row"
    >
      {colors.map((color) => (
        <div
          key={color}
          className="w-[38px] h-[30px]"
          style={{
            backgroundColor: color,
            width: `${100/colors.length}%`
          }}
        />
      ))}
    </button>
  );
}
