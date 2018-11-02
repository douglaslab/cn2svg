#!/usr/bin/env python3

import argparse
from math import floor
import os
import sys
from typing import Dict, List, Tuple

from pysvg.shape import Circle, Path, Rect, Polygon  # Line, Polyline
from pysvg.structure import G, Svg
from pysvg.text import Text

import cadnano
from cadnano.document import Document


DEFAULT_SLICE_SCALE = 10
DEFAULT_PATH_SCALE = 10


class CadnanoDocument(object):
    """
    Opens and Parses a Cadnano file. Provides accessors for cadnano design
    parameters needed for SVG conversion.
    """
    def __init__(self, cnjsonpath: str) -> None:
        super(CadnanoDocument, self).__init__()
        self.cnjsonpath = cnjsonpath

        app = cadnano.app()
        doc = app.document = Document()
        doc.readFile(cnjsonpath)
        self.part = part = doc.activePart()
        self.part_props = part_props = part.getModelProperties().copy()
        self.vh_order = part_props['virtual_helix_order']
        self.vh_props, self.vh_origins, _ = part.helixProperties()
        self.vh_radius = part.radius()
        self.max_vhelix_length = max(self.vh_props['length'])

        # Determine part coordinate boundaries
        xLL, yLL, xUR, yUR = part.getVirtualHelixOriginLimits()
        print(xLL, yLL, xUR, yUR)
        self.slice_width = xUR-xLL + self.vh_radius * 8
        self.slice_height = yUR-yLL + self.vh_radius * 6
        self.x_offset = self.slice_width/2
        self.y_offset = self.slice_height/2
    # end def

    def getSliceDimensions(self) -> Tuple:
        return self.slice_width, self.slice_height
    # end def

    def getOligoList(self) -> List:
        oligo_list = []
        for oligo in self.part.oligos():
            color = oligo.getColor()
            strands = []
            strand5p = oligo.strand5p()
            for s in strand5p.generator3pStrand():
                strands.append([s.idNum(), s.idx5Prime(), s.idx3Prime(), s.isForward()])
            oligo_list.append([color, strands])
        return oligo_list
    # end def

    def getOligoEndpointsList(self):
        ends5p = []
        ends3p = []
        for oligo in self.part.oligos():
            if not oligo.isCircular():
                color = oligo.getColor()
                strand5p = oligo.strand5p()
                idnum5p = strand5p.idNum()
                idx5p = strand5p.idx5Prime()
                isfwd5p = strand5p.isForward()
                ends5p.append([color, idnum5p, idx5p, isfwd5p])
                strand3p = oligo.strand3p()
                idnum3p = strand3p.idNum()
                idx3p = strand3p.idx3Prime()
                isfwd3p = strand3p.isForward()
                ends3p.append([color, idnum3p, idx3p, isfwd3p])
        return ends5p, ends3p
    # end def

    def insertions(self):
        pass
    # end def

    def skips(self):
        pass
    # end def

# end class


class CadnanoSliceSvg(object):
    """
    Generate a Sliceview SVG for the given Cadnano document cn_doc.
    """
    def __init__(self, cn_doc, output_path, scale=DEFAULT_SLICE_SCALE):
        super(CadnanoSliceSvg, self).__init__()
        self.cn_doc = cn_doc
        self.output_path = output_path
        self._scale = scale
        self._slice_radius_scaled = cn_doc.vh_radius*scale
        self._slice_vh_fontsize = floor(2*self._slice_radius_scaled*0.75)
        self.g_slicevirtualhelices = self.makeSliceVhGroup()
        self.g_slicevirtualhelixlabels = self.makeSliceVhLabelGroup()
        w, h = self.cn_doc.getSliceDimensions()
        scaled_w = w*scale
        scaled_h = h*scale
        self.makeSliceSvg(scaled_w, scaled_h)
    # end def

    def makeSliceVhGroup(self) -> G:
        """
        Creates and returns a 'G' object for Slice Virtual Helices.
        """
        vh_radius = self.cn_doc.vh_radius
        _SLICE_SCALE = self._scale
        g = G()
        g.setAttribute("id", "VirtualHelices")
        for id_num in self.cn_doc.vh_order[::-1]:
            vh_x, vh_y, vh_z = self.cn_doc.vh_origins[id_num]
            x = vh_x + self.cn_doc.x_offset
            y = -vh_y + self.cn_doc.y_offset
            c = Circle(x*_SLICE_SCALE, y*_SLICE_SCALE, vh_radius*_SLICE_SCALE)
            circle_style = 'fill:#f2ca9a; stroke:#cc6600; stroke-width:1;'
            c.set_style(circle_style)
            c.setAttribute("id", "circle_"+self.cn_doc.vh_props['name'][id_num])
            g.addElement(c)
        return g
    # end def

    def makeSliceVhLabelGroup(self) -> G:
        """
        Creates and returns a 'G' object for Slice VirtualHelix Labels.
        """
        _SLICE_SCALE = self._scale
        g = G()
        g.setAttribute("id", "VirtualHelixLabels")
        g.setAttribute("font-size", "%s" % self._slice_vh_fontsize)
        # g.setAttribute("font-family", "'Source Sans Pro', sans-serif")
        g.setAttribute("font-family", "'SourceSansPro-Regular'")
        # g.setAttribute("font-family", "sans-serif")
        g.setAttribute("text-anchor", "middle")
        for id_num in self.cn_doc.vh_order[::-1]:
            vh_x, vh_y, vh_z = self.cn_doc.vh_origins[id_num]
            x = vh_x + self.cn_doc.x_offset
            y = -vh_y + self.cn_doc.y_offset + self.cn_doc.vh_radius/2.
            t = Text('%s' % id_num, x*_SLICE_SCALE, y*_SLICE_SCALE - 1)
            t.setAttribute("id", "label_"+self.cn_doc.vh_props['name'][id_num])
            g.addElement(t)
        return g
    # end def

    def makeSliceSvg(self, width, height) -> Svg:
        slice_svg = Svg(width=width, height=height)
        viewbox = "0 0 %s %s" % (width, height)
        slice_svg.set_viewBox(viewbox)
        slice_svg.set_preserveAspectRatio("xMidYMid meet")
        slice_svg.setAttribute("id", "Cadnano_Slice")  # Main layer name
        slice_svg.addElement(self.g_slicevirtualhelices)  # bottom layer
        slice_svg.addElement(self.g_slicevirtualhelixlabels)  # top layer
        slice_svg.save(self.output_path)
        return slice_svg
    # end def
# end class


class CadnanoPathSvg(object):
    """
    Generate a Pathview SVG for the given Cadnano document cn_doc.
    """
    PATH_X_PADDING = 40
    PATH_Y_PADDING = 40

    def __init__(self, cn_doc, output_path, scale=DEFAULT_PATH_SCALE):
        super(CadnanoPathSvg, self).__init__()
        self.cn_doc = cn_doc
        self.output_path = output_path
        self._scale = scale
        self._path_radius_scaled = cn_doc.vh_radius*scale
        self._path_vh_fontsize = floor(2*self._path_radius_scaled*0.75)
        self._path_vh_margin = self._path_radius_scaled*5
        self._base_width = self.base_height = self._path_radius_scaled

        self.g_pathvirtualhelices = self.makePathVhGroup()
        self.g_pathvirtualhelixlabels = self.makePathVhLabelGroup()
        self.g_pathgridlines = self.makePathGridlinesGroup()
        self.g_patholigos = self.makePathOligosGroup()
        self.g_pathendpoints = self.makePathEndpointsGroup()
        w = self.PATH_X_PADDING*3 + cn_doc.max_vhelix_length*self._base_width
        h = len(cn_doc.vh_order)*self._path_vh_margin + self.PATH_Y_PADDING/2
        self.path_svg = self.makePathSvg(width=w, height=h)

        # self.y_coords = self.mapIdnumsToYcoords()
    # end def

    def mapIdnumsToYcoords(self) -> Dict:
        d = {}
        for i in range(len(self.cn_doc.vh_order)):
            id_num = self.cn_doc.vh_order[i]
            # y0 = self.PATH_Y_PADDING + self._path_vh_margin*i
            # y1 = y0 + self.base_height
            # d[id_num] = [y0, y1]
            y = self.PATH_Y_PADDING + self._path_vh_margin*i
            d[id_num] = y
        return d

    def makePathVhGroup(self) -> G:
        """
        Creates and returns a 'G' object for Path Virtual Helices.
        """
        g = G()
        g.setAttribute("id", "VirtualHelices")
        for i in range(len(self.cn_doc.vh_order)):
            id_num = self.cn_doc.vh_order[i]
            x = self.PATH_X_PADDING
            y = self.PATH_Y_PADDING + self._path_vh_margin*i
            c = Circle(x, y, self.cn_doc.vh_radius*self._scale)
            circle_style = 'fill:#f2ca9a; stroke:#cc6600; stroke-width:1;'
            c.set_style(circle_style)
            c.setAttribute("id", "circle_"+self.cn_doc.vh_props['name'][id_num])
            g.addElement(c)
        return g
    # end def

    def makePathVhLabelGroup(self) -> G:
        """
        Creates and returns a 'G' object for Path VirtualHelix Labels.
        """
        g = G()
        g.setAttribute("id", "VirtualHelixLabels")
        g.setAttribute("font-size", "%s" % self._path_vh_fontsize)
        g.setAttribute("font-family", "'SourceSansPro-Regular'")
        # g.setAttribute("font-family", "'Source Sans Pro', sans-serif")
        # g.setAttribute("font-family", "sans-serif")
        g.setAttribute("text-anchor", "middle")
        for i in range(len(self.cn_doc.vh_order)):
            id_num = self.cn_doc.vh_order[i]
            x = self.PATH_X_PADDING
            y = self.PATH_Y_PADDING + self._path_vh_margin*i + self._path_radius_scaled/2.
            t = Text('%s' % id_num, x, y)
            t.setAttribute("id", "label_"+self.cn_doc.vh_props['name'][id_num])
            g.addElement(t)
        return g
    # end def

    def makePathGridlinesGroup(self) -> G:
        size = self.cn_doc.max_vhelix_length
        print(size)
        _BW = self._base_width
        _BH = self.base_height
        g = G()
        g.setAttribute("id", "GridLines")
        for i in range(len(self.cn_doc.vh_order)):
            id_num = self.cn_doc.vh_order[i]
            x = self.PATH_X_PADDING + self._path_radius_scaled*3
            y = self.PATH_Y_PADDING + self._path_vh_margin*i
            h_lines = " ".join("M %s,%s h %s" % (x, y-_BH+j*_BH, _BW*size) for j in range(3))
            v_lines = " ".join("M %s,%s v %s" % (x+j*_BW, y-_BH, _BH*2) for j in range(size+1))
            p = Path(h_lines + " " + v_lines)
            grid_style = 'fill:none; stroke:#666666; stroke-width:0.25'
            p.set_style(grid_style)
            p.setAttribute("id", "grid_"+self.cn_doc.vh_props['name'][id_num])
            g.addElement(p)
        return g
    # end def

    def makePathOligosGroup(self) -> G:
        id_coords = self.mapIdnumsToYcoords()
        oligo_list = self.cn_doc.getOligoList()
        _BW = self._base_width
        _BH = self.base_height
        _pX = self.PATH_X_PADDING + self._path_radius_scaled*3 + _BW/2
        g = G()
        g.setAttribute("id", "Oligos")
        i = 0
        for color, strands in oligo_list:
            prev_id, prev5, prev3, prevX, prevY = None, None, None, 0, 0
            path_lines = []
            for id_num, idx5, idx3, isfwd in strands:
                x = _pX + idx5*_BW
                if idx5 < idx3:  # top strand
                    y = id_coords[id_num] - _BH/2
                else:
                    y = id_coords[id_num] + _BH/2
                dx = (idx3-idx5)*_BW

                if prev_id is not None and id_num != prev_id:
                    # c dx1 dy1, dx2 dy2, dx dy
                    if isfwd:
                        x1 = x + abs(y-prevY)*0.03
                    else:
                        x1 = x - abs(y-prevY)*0.03
                    y1 = (y+prevY)/2
                    path_lines.append("Q %s %s, %s %s" % (x1, y1, x, y))
                    path_lines.append("h %s" % dx)
                else:
                    path_lines.append("M %s, %s h %s" % (x, y, dx))
                prev_id, prev5, prev3, prevX, prevY = id_num, idx5, idx3, x, y
            p = Path(" ".join(path_lines))
            oligo_style = 'fill:none; stroke:%s; stroke-width:1' % color
            p.set_style(oligo_style)
            p.setAttribute("id", "oligo_%s" % i)
            p.setAttribute("stroke-linejoin", "round")
            i += 1
            g.addElement(p)
        return g
    # end def

    def makePathEndpointsGroup(self) -> G:
        _BW = self._base_width
        _BH = self.base_height
        _pX = self.PATH_X_PADDING + self._path_radius_scaled*3
        id_coords = self.mapIdnumsToYcoords()
        ends5p, ends3p = self.cn_doc.getOligoEndpointsList()
        # print(ends5p, ends3p)
        g = G()
        g.setAttribute("id", "Endpoints")
        i = 0
        for color, idnum5p, idx5p, isfwd5p in ends5p:
            if isfwd5p:
                x = _pX + idx5p*_BW + _BW*.25
                y = id_coords[idnum5p] - _BH + _BH*.05
            else:
                x = _pX + idx5p*_BW
                y = id_coords[idnum5p] + _BH*.05
            r = Rect(x=x, y=y, width=_BW*0.75, height=_BH*.9)
            endpoint_style = 'fill:%s; stroke:none; stroke-width:0.5' % color
            r.set_style(endpoint_style)
            r.setAttribute("id", "end5p_%s" % i)
            i += 1
            g.addElement(r)
        j = 0
        for color, idnum3p, idx3p, isfwd3p in ends3p:
            if isfwd3p:
                x = _pX + idx3p*_BW
                y = id_coords[idnum3p] - _BH
                x1, y1 = x, y
                x2, y2 = x1+_BW*0.75, y1+_BH/2
                x3, y3 = x, y1+_BH
            else:
                x = _pX + idx3p*_BW
                y = id_coords[idnum3p]
                x1, y1 = x+_BW, y
                x2, y2 = x+_BW*0.25, y1+_BH/2
                x3, y3 = x+_BW, y+_BH
            pts = "%s,%s %s,%s %s,%s %s,%s" % (x1, y1, x2, y2, x3, y3, x1, y1)
            p = Polygon(points=pts)
            end3p_style = 'fill:%s; stroke:none; stroke-width:0.5' % color
            p.set_style(end3p_style)
            p.setAttribute("id", "end3p_%s" % j)
            j += 1
            g.addElement(p)
        return g
    # end def

    def makePathSvg(self, width, height) -> Svg:
        viewbox  = "0 0 %s %s" % (width, height)
        path_svg = Svg(width=width, height=height)
        path_svg.set_viewBox(viewbox)
        path_svg.set_preserveAspectRatio("xMinYMid meet")
        path_svg.setAttribute("id", "Cadnano_Path")  # Main layer name
        path_svg.addElement(self.g_pathgridlines)  # bottom layer
        path_svg.addElement(self.g_patholigos)
        path_svg.addElement(self.g_pathendpoints)
        path_svg.addElement(self.g_pathvirtualhelices)
        path_svg.addElement(self.g_pathvirtualhelixlabels)  # top layer
        path_svg.save(self.output_path)
        return path_svg
    # end def
# end class


def main():
    parser  = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', '-i', type=str, required=True, nargs='+', help='Cadnano json file(s)')
    parser.add_argument('--output', '-o', type=str, help='Output directory')
    args    = parser.parse_args()

    if args.input is None:
        parser.print_help()
        sys.exit('Input file not specified')

    json_files          = [filename for filename in args.input if filename.endswith('.json')]
    output_directory    = args.output

    if not json_files:
        parser.print_help()
        sys.exit('Input file(s) is/are not JSON files')

    for design in json_files:
        basename        = os.path.splitext(os.path.basename(design))[0]
        base_path       = os.path.splitext(design)[0]
        cndoc           = CadnanoDocument(design)
        output_slice    = os.path.join(output_directory, '%s_slice.svg' % basename) if output_directory and os.path.exists(output_directory)\
                            else '%s_slice.svg' % base_path
        output_path     = os.path.join(output_directory, '%s_path.svg' % basename) if output_directory and os.path.exists(output_directory) \
                            else '%s_path.svg' % base_path
        CadnanoSliceSvg(cndoc, output_slice)
        CadnanoPathSvg(cndoc, output_path)
# end def


if __name__ == "__main__":
    main()
