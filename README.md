# cn2svg - Cadnano SVG exporter

cn2svg a rewrite of the SVG export tool found in the original cadnano ([Svg.as](https://github.com/sdouglas/cadnano/blob/master/edu/harvard/med/cadnano/data/Svg.as)).

## Installation

`pip3 install cn2svg`

## Example Usage

`cn2svg -i cadnanofile.json -o outputdir`

The output should include two separate files: `cadnanofile_slice.svg` and `cadnanofile_path.svg`.

## Dependencies

- [cadnano2.5](https://github.com/douglaslab/cadnano2.5)
- [pysvg-py3](https://github.com/alorence/pysvg-py3)

## VENV Installation

```
mkvirtualenv myvenv
pip3 install PyQt5==5.10.0 pandas termcolor pysvg-py3
git clone https://github.com/douglaslab/cadnano2.5
cd cadnano2.5/
python3 setup.py install
```

## About

Rendering Cadnano designs in SVG format can be useful for making figures and schematics. However, SVG support was very limited in Cadnano2 because it relied on the [QtSvgGenerator](https://doc.qt.io/qt-5/qsvggenerator.html) class to provide automatic conversion. The conversion did not retain element groups based on the Cadnano data structures, and most path strokes were scaled incorrectly.

As part of developing the new Cadnano Toolkit, we've rewritten the original SVG exporter as a standalone Python script, `cn2svg`. It provides a few new features under the hood, such as hooks for dynamic path styling via unique element id attributes, and element instancing using the `defs` and `use` tags.

## Citing

If you use `cn2svg`, please cite [doi.org/10.1093/nar/gkp436](https://doi.org/10.1093/nar/gkp436).
