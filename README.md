# cn2svg - Cadnano SVG exporter

cn2svg a rewrite of the SVG export tool found in the original cadnano ([Svg.as](https://github.com/sdouglas/cadnano/blob/master/edu/harvard/med/cadnano/data/Svg.as)).

## About

Rendering Cadnano designs in SVG format can be useful for making figures and schematics. However, SVG support was very limited in Cadnano2 because it relied on the [QtSvgGenerator](https://doc.qt.io/qt-5/qsvggenerator.html) class to provide automatic conversion. The conversion did not retain element groups based on the Cadnano data structures, and most path strokes were scaled incorrectly.

As part of developing the new Cadnano Toolkit, we've rewritten the original SVG exporter as a standalone Python script, `cn2svg`. It provides a few new features under the hood, such as hooks for dynamic path styling via unique element id attributes, and element instancing using the `defs` and `use` tags.


## Installation

`pip3 install cn2svg`

## Example Usage

`cn2svg -i cadnanofile.json -o outputdir`

The output should include two separate files: `cadnanofile_slice.svg` and `cadnanofile_path.svg`.

## Sequence support

To overlay DNA sequences onto path oligos, use the following option:

`cn2svg -i cadnanofile.json --seq scaffold.txt`

Sequence alignment and styling is currently formatted for browser viewing.

## Exporting for Adobe Illustrator CS6

By default, cn2svg exports SVGs formatted for web browser viewing. Unfortunately, if we [stylize the document](https://graphicdesign.stackexchange.com/questions/36168/is-there-any-way-to-set-fallback-font-families-in-illustrator-svg) with custom fonts via the `font-family` attribute, we run afoul of [a bug](https://forums.adobe.com/thread/1326594) in Adobe Illustrator CS6, which cannot open the file.

The current workaround is to support the export of a custom SVG files that open in CS6. To invoke this option, use the `--cs6` flag:

`cn2svg -i cadnanofile.json --cs6`

Expected output: `cadnanofile_slice_cs6.svg` and `cadnanofile_path_cs6.svg`.

Path oligo sequences will not match the browser appearance in CS6. You can improve alignment by adjusting the font size and tracking.

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


## Citing

If you use `cn2svg`, please cite [doi.org/10.1093/nar/gkp436](https://doi.org/10.1093/nar/gkp436).


## Development

We use a pyproject.toml-based [build process](https://pip.pypa.io/en/stable/reference/build-system/pyproject-toml/) in pip. This workflow was tested with python 3.12 on macOS in November 2023.

**Setup a dev environment (Mac or Linux)**

* Create a virtualenv: `python3 -m venv ~/virtualenvs/cn2svgdev` 
* Activate virtualenv: `source ~/virtualenvs/cn2svgdev/bin/activate`
* Clone repo: `git clone git@github.com:douglaslab/cn2svg.git`
* Change directory: `cd cn2svg`
* Make desired code edits
* Build and install in [editable mode](https://pip.pypa.io/en/stable/cli/pip_install/#cmdoption-e): `pip install -e .` 
* Test: `cn2svg -i cadnanofile.json`
* Repeat previous 3 steps as needed
