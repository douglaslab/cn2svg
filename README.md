# cn2svg - Cadnano SVG exporter

cn2svg a rewrite of the SVG export tool found in the original cadnano ([Svg.as](https://github.com/sdouglas/cadnano/blob/master/edu/harvard/med/cadnano/data/Svg.as)).

## Installation

`pip3 install cn2svg`

## Example Usage

`cn2svg -i cadnanofile.json -o outputdir`

The output should include two separate files: `cadnanofile_slice.svg` and `cadnanofile_path.svg`.

## Dependencies

- [cadnano2.5](https://github.com/cadnano/cadnano2.5)
- [pysvg-py3](https://github.com/alorence/pysvg-py3)

## Notes

Rendering Cadnano designs in SVG format can be useful for making figures and schematics. SVG support was rather limited in Cadnano2 because it relied on the [QtSvgGenerator](https://doc.qt.io/qt-5/qsvggenerator.html) class to provide automatic conversion. The conversion did not retain any grouping based on the Cadnano data structures, and the resulting path strokes were scaled incorrectly, and required extensive manual editing.

## Citing

If you use `cn2svg`, please cite [doi.org/10.1093/nar/gkp436](https://doi.org/10.1093/nar/gkp436).
