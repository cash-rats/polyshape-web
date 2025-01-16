import { useMemo } from 'react';
import chroma from 'chroma-js';

type BrewerPaletteName = keyof typeof chroma.brewer;

// Produces list of color quantiles based on chroma's predefined color domains, we'll have 9 quantiles per domain.
export const useColorBrewer = () => {
  const userValues = useMemo(() => [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], []);

  const colorBrewsers = Object.keys(chroma.brewer).reduce<Record<string, chroma.Scale>>((acc, colorDomain) => {
    acc[colorDomain] = chroma.scale(colorDomain as BrewerPaletteName).domain(userValues, userValues.length, 'quantiles');
    return acc;
  }, {});

  return {
    colorBrewsers,
  };
}
