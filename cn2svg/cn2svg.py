import argparse
import json
import os
import sys
from collections import defaultdict
from contextlib import redirect_stdout
from math import floor
from typing import Dict, List, Tuple

from pysvg.shape import Circle, Path, Rect, Polygon, Line  # Polyline
from pysvg.structure import Defs, G, Svg, Use
from pysvg.text import Text

import cadnano
from cadnano.document import Document

DEFAULT_ORTHO_SCALE = 10
DEFAULT_PATH_SCALE = 10


class CadnanoDocument(object):
    """
    Opens and Parses a Cadnano file. Provides accessors for cadnano design
    parameters needed for SVG conversion.
    """
    def __init__(self, cnjsonpath: str, sequence:str) -> None:
        super(CadnanoDocument, self).__init__()
        self.cnjsonpath = cnjsonpath

        app = cadnano.app()
        doc = app.document = Document()
        doc.readFile(cnjsonpath)
        self.part = part = doc.activePart()
        self.part_props = part_props = part.getModelProperties().copy()
        self.vh_order = part_props['virtual_helix_order']
        self.vh_props, self.vh_origins = part.helixPropertiesAndOrigins()
        self.vh_radius = part.radius()
        self.max_vhelix_length = max(self.vh_props['length'])
        self.insertions, self.skips = self.getInsertionsAndSkips()

        self.sequence_applied = False
        self.max_oligo_length = 0
        self.applySequenceToSameLengthOligos(sequence)
    # end def

    def getOrthoDimensions(self) -> Tuple:
        return self.ortho_width, self.ortho_height
    # end def

    def applySequenceToSameLengthOligos(self, sequence) -> None:
        seq_length = len(sequence)
        for oligo in self.part.oligos():
            oligo_len = oligo.length()
            if oligo_len > self.max_oligo_length:
                self.max_oligo_length = oligo_len
            if oligo_len == seq_length:
                s = oligo.strand5p()
                fwd = 'fwd' if s.isForward() else 'rev'
                # print('Applying sequence at %s[%s] (%s strand)' % (s.idNum(), s.idx5Prime(), fwd))
                oligo.applySequence(sequence, use_undostack=False)
                self.sequence_applied = True

    def getOligoList(self) -> List:
        oligo_list = []
        for oligo in self.part.oligos():
            color = oligo.getColor()
            strands = []
            strand5p = oligo.strand5p()
            is_staple = None 
            for s in strand5p.generator3pStrand():
                strands.append([s.idNum(), s.idx5Prime(), s.idx3Prime(), s.isForward(), s.sequence()])
                if is_staple is None:
                    is_staple = not s.isForward() if (s.idNum() % 2 == 0) else s.isForward()
            oligo_list.append([is_staple, color, oligo.isCircular(), strands])
        return oligo_list
    # end def

    def getOligoEndpointsList(self) -> Tuple:
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

    def getInsertionsAndSkips(self) -> Tuple:
        part = self.part
        insertion_dict = defaultdict(list)
        skip_dict = defaultdict(list)
        part_insertions = part.insertions()
        for id_num in self.vh_order:
            vh_insertions = part_insertions[id_num]
            for idx, insertion in vh_insertions.items():
                if insertion.isSkip():
                    skip_dict[id_num].append((idx, insertion.length()))
                else:
                    has_fwd, has_rev = part.hasStrandAtIdx(id_num, idx)
                    fwd_col, rev_col = '#cccccc', '#cccccc'  # Defaults
                    if has_fwd:
                        fwd_strand = part.getStrand(True, id_num, idx)
                        fwd_col = fwd_strand.getColor()
                    if has_rev:
                        rev_strand = part.getStrand(False, id_num, idx)
                        rev_col = rev_strand.getColor()
                    insertion_dict[id_num].append((idx, insertion.length(), fwd_col, rev_col))
        return insertion_dict, skip_dict
    # end def
# end class


class CadnanoOrthoSvg(object):
    """
    Generate a Orthoview SVG for the given Cadnano document cn_doc.
    """
    def __init__(self, cn_doc, output_path, scale=DEFAULT_ORTHO_SCALE):
        super(CadnanoOrthoSvg, self).__init__()
        self.cn_doc = cn_doc
        self.output_path = output_path
        self.w = None
        self.h = None
        self._scale = scale
        self._ortho_radius_scaled = cn_doc.vh_radius*scale
        self._ortho_vh_fontsize = floor(2*self._ortho_radius_scaled*0.75)
        self.g_orthovirtualhelices, self.w, self.h, transform = self.makeOrthoVhGroup()
        self.g_orthovirtualhelixlabels = self.makeOrthoVhLabelGroup(transform)
        self.makeOrthoSvg()
    # end def

    def makeOrthoVhGroup(self) -> Tuple:
        """
        Creates and returns a 'G' object for Ortho Virtual Helices.
        """
        _ORTHO_SCALE = self._scale
        vh_radius = self.cn_doc.vh_radius * _ORTHO_SCALE
        g = G()
        g.setAttribute('id', "VirtualHelices")
        minx = 99999
        miny = 99999
        maxx = -99999
        maxy = -99999
        for id_num in self.cn_doc.vh_order[::-1]:
            vh_x, vh_y = self.cn_doc.vh_origins[id_num]
            x = vh_x * _ORTHO_SCALE  # + self.cn_doc.x_offset
            y = -vh_y * _ORTHO_SCALE  # + self.cn_doc.y_offset
            minx = min(minx, x)
            miny = min(miny, y)
            maxx = max(maxx, x)
            maxy = max(maxy, y)
            c = Circle(x, y, vh_radius)
            circle_style = 'fill:#f2ca9a; stroke:#cc6600; stroke-width:1;'
            c.set_style(circle_style)
            c.setAttribute('id', "circle_"+self.cn_doc.vh_props['name'][id_num])
            g.addElement(c)
        minx = minx-vh_radius*2
        miny = miny-vh_radius*2
        width = round(maxx-minx+vh_radius*2, 3)
        height = round(maxy-miny+vh_radius*2, 3)
        transform = 'translate(%.3f %.3f)' % (-minx, -miny)
        g.set_transform(transform)
        return g, width, height, transform
    # end def

    def makeOrthoVhLabelGroup(self, transform) -> G:
        """
        Creates and returns a 'G' object for Ortho VirtualHelix Labels.
        """
        _ORTHO_SCALE = self._scale
        g = G()
        g.setAttribute('id', "VirtualHelixLabels")
        g.setAttribute('font-size', "%s" % self._ortho_vh_fontsize)
        g.setAttribute("font-family", "'Source Sans Pro', sans-serif")
        g.setAttribute("text-anchor", "middle")
        for id_num in self.cn_doc.vh_order[::-1]:
            vh_x, vh_y = self.cn_doc.vh_origins[id_num]
            x = vh_x  # + self.cn_doc.x_offset
            y = -vh_y + self.cn_doc.vh_radius/2.
            t = Text('%s' % id_num, x*_ORTHO_SCALE, y*_ORTHO_SCALE - 1)
            t.setAttribute('id', "label_"+self.cn_doc.vh_props['name'][id_num])
            g.addElement(t)
        g.set_transform(transform)
        return g
    # end def

    def makeOrthoSvg(self) -> Svg:
        ortho_svg = Svg(width=self.w, height=self.h)
        viewbox = "0 0 %s %s" % (self.w, self.h)
        ortho_svg.set_viewBox(viewbox)
        ortho_svg.set_preserveAspectRatio("xMidYMid meet")
        ortho_svg.setAttribute('id', "Cadnano_Ortho")  # Main layer name
        ortho_svg.addElement(self.g_orthovirtualhelices)  # bottom layer
        ortho_svg.addElement(self.g_orthovirtualhelixlabels)  # top layer
        ortho_svg.save(self.output_path)
    # end def
# end class


class CadnanoPathSvg(object):
    """
    Generate a Pathview SVG for the given Cadnano document cn_doc.
    """
    PATH_X_PADDING = 40
    PATH_Y_PADDING = 40

    def __init__(self, cn_doc, output_path, heatmap, scale=DEFAULT_PATH_SCALE):
        super(CadnanoPathSvg, self).__init__()
        self.cn_doc = cn_doc
        self.output_path = output_path
        self._heatmap = heatmap
        self.w = None
        self.h = None
        self._scale = scale
        self._path_radius_scaled = cn_doc.vh_radius*scale
        self._path_vh_fontsize = floor(2*self._path_radius_scaled*0.75)
        self._path_vh_margin = self._path_radius_scaled*(3 if heatmap else 5)
        self._base_width = self._base_height = self._path_radius_scaled

        self.defs = self.makePathDefs()

        self.g_pathvirtualhelices = self.makePathVhGroup()
        self.g_pathvirtualhelixlabels = self.makePathVhLabelGroup()
        self.g_pathgridlines = self.makePathGridlinesGroup()
        self.g_patholigos = self.makePathOligosGroup()
        self.g_pathstaplexos = self.makePathStapleCrossoversGroup()
        self.g_pathendpoints = self.makePathEndpointsGroup()
        self.g_pathinsertions = self.makePathInsertionsGroup()
        self.g_pathskips = self.makePathSkipsGroup()
        if cn_doc.sequence_applied:
            self.g_pathsequences = self.makePathSequencesGroup()

        self.w = round(self.PATH_X_PADDING*2.5 + cn_doc.max_vhelix_length*self._base_width, 3)
        self.h = round(len(cn_doc.vh_order)*self._path_vh_margin + self.PATH_Y_PADDING/2, 3)
        self.path_svg = self.makePathSvg()
    # end def

    def mapIdnumsToYcoords(self) -> Dict:
        d = {}
        for i in range(len(self.cn_doc.vh_order)):
            id_num = self.cn_doc.vh_order[i]
            # y0 = self.PATH_Y_PADDING + self._path_vh_margin*i
            # y1 = y0 + self._base_height
            # d[id_num] = [y0, y1]
            y = self.PATH_Y_PADDING + self._path_vh_margin*i
            d[id_num] = y
        return d

    def makePathVhGroup(self) -> G:
        """
        Creates and returns a 'G' object for Path Virtual Helices.
        """
        g = G()
        g.setAttribute('id', "VirtualHelices")
        for i in range(len(self.cn_doc.vh_order)):
            id_num = self.cn_doc.vh_order[i]
            x = self.PATH_X_PADDING
            y = self.PATH_Y_PADDING + self._path_vh_margin*i
            c = Circle(x, y, self.cn_doc.vh_radius*self._scale)
            circle_style = 'fill:#f2ca9a; stroke:#cc6600; stroke-width:1;'
            c.set_style(circle_style)
            c.setAttribute('id', "circle_"+self.cn_doc.vh_props['name'][id_num])
            g.addElement(c)
        return g
    # end def

    def makePathVhLabelGroup(self) -> G:
        """
        Creates and returns a 'G' object for Path VirtualHelix Labels.
        """
        g = G()
        g.setAttribute('id', "VirtualHelixLabels")
        g.setAttribute('font-size', "%s" % self._path_vh_fontsize)
        g.setAttribute("font-family", "'Source Sans Pro', sans-serif")
        g.setAttribute("text-anchor", "middle")
        for i in range(len(self.cn_doc.vh_order)):
            id_num = self.cn_doc.vh_order[i]
            x = self.PATH_X_PADDING
            y = self.PATH_Y_PADDING + self._path_vh_margin*i + self._path_radius_scaled/2.
            t = Text('%s' % id_num, x, y)
            t.setAttribute('id', "label_"+self.cn_doc.vh_props['name'][id_num])
            g.addElement(t)
        return g
    # end def

    def makePathGridlinesGroup(self) -> G:
        size = self.cn_doc.max_vhelix_length
        _BW = self._base_width
        _BH = self._base_height
        g = G()
        g.setAttribute('id', "GridLines")
        for i in range(len(self.cn_doc.vh_order)):
            id_num = self.cn_doc.vh_order[i]
            x = self.PATH_X_PADDING + self._path_radius_scaled*3
            y = self.PATH_Y_PADDING + self._path_vh_margin*i
            h_lines = " ".join("M %s,%s h %s" % (x, y-_BH+j*_BH, _BW*size) for j in range(3))
            v_lines = " ".join("M %s,%s v %s" % (x+j*_BW, y-_BH, _BH*2) for j in range(int(size)+1))
            p = Path(h_lines + " " + v_lines)
            grid_style = 'fill:none; stroke:#666666; stroke-width:0.25'
            p.set_style(grid_style)
            p.setAttribute('id', "grid_"+self.cn_doc.vh_props['name'][id_num])
            g.addElement(p)
        return g
    # end def

    def makePathOligosGroup(self) -> G:
        heatmap = self._heatmap
        id_coords = self.mapIdnumsToYcoords()
        oligo_list = self.cn_doc.getOligoList()
        _BW = self._base_width
        _BH = self._base_height
        _pX = self.PATH_X_PADDING + self._path_radius_scaled*3 + _BW/2
        g = G()
        g.setAttribute('id', "Oligos")
        i = 0
        for is_staple, color, is_circular, strands in oligo_list:
            prev_id, prev5, prev3, prevX, prevY = None, None, None, 0, 0
            path_lines = []
            if is_circular:
                strands.append(strands[0])
            for id_num, idx5, idx3, isfwd, _ in strands:
                x = _pX + idx5*_BW
                is_top_strand = idx5 < idx3
                if is_top_strand:  # top strand
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
                    if heatmap and is_staple:
                        dy = _BH/2 if idx5 > idx3 else -_BH/2
                        dy1 = _BH/2 if prev5 > prev3 else -_BH/2
                        path_lines.append("M %s, %s" % (x, y))
                    else:
                        path_lines.append("Q %s %s, %s %s" % (x1, y1, x, y))
                    path_lines.append("h %s" % dx)
                else:
                    path_lines.append("M %s, %s h %s" % (x, y, dx))
                prev_id, prev5, prev3, prevX, prevY = id_num, idx5, idx3, x, y
            p = Path(" ".join(path_lines))
            sw = 14 if (heatmap and is_staple) else 2
            oligo_style = 'fill:none; stroke:%s; stroke-width:%s' % (color,sw)
            p.set_style(oligo_style)
            p.setAttribute('id', "oligo_%s" % i)
            p.setAttribute("stroke-linejoin", "round")
            i += 1
            g.addElement(p)
        return g
    # end def

    def makePathStapleCrossoversGroup(self) -> G:
        heatmap = self._heatmap
        id_coords = self.mapIdnumsToYcoords()
        oligo_list = self.cn_doc.getOligoList()
        _BW = self._base_width
        _BH = self._base_height
        _pX = self.PATH_X_PADDING + self._path_radius_scaled*3 + _BW/2
        g = G()
        g.setAttribute('id', "OligoCrossovers")
        i = 0
        for is_staple, color, is_circular, strands in oligo_list:
            prev_id, prev5, prev3, prevX, prevY = None, None, None, 0, 0
            path_lines = []
            if is_circular:
                strands.append(strands[0])
            for id_num, idx5, idx3, isfwd, _ in strands:
                x = _pX + idx5*_BW
                is_top_strand = idx5 < idx3
                if is_top_strand:  # top strand
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
                    # path_lines.append("M %s, %s" % (x1,y1))
                    path_lines.append("Q %s %s, %s %s" % (x1, y1, x, y))
                    path_lines.append("M %s, %s" % (x+dx,y))
                else:
                    path_lines.append("M %s, %s" % (x+dx, y))
                
                prev_id, prev5, prev3, prevX, prevY = id_num, idx5, idx3, x, y
            p = Path(" ".join(path_lines))
            # sw = 2 # if (heatmap and is_staple) else 2
            oligo_style = 'fill:none; stroke:%s; stroke-width:2' % (color)
            p.set_style(oligo_style)
            p.setAttribute('id', "oligo_stapxo_%s" % i)
            p.setAttribute("stroke-linejoin", "round")
            i += 1
            g.addElement(p)
        return g
    # end def


    def makePathEndpointsGroup(self) -> G:
        heatmap = self._heatmap
        _BW = self._base_width
        _BH = self._base_height
        _pX = self.PATH_X_PADDING + self._path_radius_scaled*3
        id_coords = self.mapIdnumsToYcoords()
        ends5p, ends3p = self.cn_doc.getOligoEndpointsList()
        g = G()
        g.setAttribute('id', "Endpoints")
        i = 0
        for color, idnum5p, idx5p, isfwd5p in ends5p:
            is_staple = not isfwd5p if (idnum5p % 2 == 0) else isfwd5p
            # Disable staple endpoints for heatmap
            if heatmap and is_staple:
                continue  

            if isfwd5p:
                x = _pX + idx5p*_BW + _BW*.25
                y = id_coords[idnum5p] - _BH + _BH*.05
            else:
                x = _pX + idx5p*_BW
                y = id_coords[idnum5p] + _BH*.05
            r = Rect(x=x, y=y, width=_BW*0.75, height=_BH*.9)
            endpoint_style = 'fill:%s; stroke:none; stroke-width:0.5' % color
            r.set_style(endpoint_style)
            r.setAttribute('id', "end5p_%s" % i)
            i += 1
            g.addElement(r)
        j = 0
        for color, idnum3p, idx3p, isfwd3p in ends3p:
            is_staple = not isfwd3p if (idnum3p % 2 == 0) else isfwd3p
            # if heatmap and is_staple:
            #     continue

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
            p.setAttribute('id', "end3p_%s" % j)
            j += 1
            g.addElement(p)
        return g
    # end def

    def makePathSequencesGroup(self) -> G:
        """
        Creates and returns a 'G' object for Path oligo sequences.
        """
        _BW = self._base_width
        _BH = self._base_height
        _pX = self.PATH_X_PADDING + self._path_radius_scaled*3
        id_coords = self.mapIdnumsToYcoords()
        oligo_list = self.cn_doc.getOligoList()

        g = G()
        g.setAttribute('id', 'PathOligoSequences')
        g.setAttribute('font-size', '8')
        g.setAttribute('font-family', 'Monaco')

        i = 0
        for is_staple, color, is_circular, strands in oligo_list:
            g_oligo = G()
            g_oligo.setAttribute('id', "oligo_seq_%s" % i)
            i += 1
            for id_num, idx5, idx3, isfwd, seq in strands:
                x = _pX + idx5*_BW + _BW/4 if isfwd else _pX + idx5*_BW - _BW/4
                y = id_coords[id_num] - _BH*1.1 if isfwd else id_coords[id_num] + _BH*2.2
                t = Text('%s' % seq, x, y)
                t.set_style('fill:%s99' % color)
                # strand_width = (idx3-idx5+1)*_BW if isfwd else (idx5-idx3+1)*_BW
                # alignment slightly better with .5 offset
                strand_width = (idx3-idx5+.5)*_BW if isfwd else (idx5-idx3+.5)*_BW
                t.setAttribute('textLength', strand_width)
                if not isfwd:
                    t.set_transform('rotate(180,%s,%s)' % (x+_BW/2, y-_BH/2))

                g_oligo.addElement(t)
            g.addElement(g_oligo)
        return g
    # end def

    def makePathDefs(self) -> Defs:
        """
        Creates Defs object for reusable `Insertion` and `Skip` elements.
        """
        _BW, _BH = self._base_width, self._base_height
        _hBW, _hBH = _BW/2, _BH/2
        _67BW, _67BH = _BW/1.5, _BH/1.5

        p = Path()
        p.appendMoveToPath(_hBW, _hBH)
        p.appendQuadraticCurveToPath(-_67BW, -_BH, 0, -_BH)
        p.appendQuadraticCurveToPath(_67BW, 0, 0, _BH)
        p.set_stroke_linecap("round")
        g_insertion = G()
        g_insertion.setAttribute('id', "insertion")
        g_insertion.addElement(p)

        d = 0.5  # small delta to compensate for round endcap
        line1 = Line(d, d, _BW-d, _BH-d, style="stroke:#cc0000; stroke-width:2")
        line2 = Line(_BW-d, d, d, _BH-d, style="stroke:#cc0000; stroke-width:2")
        line1.set_stroke_linecap("round")
        line2.set_stroke_linecap("round")
        g_skip = G()
        g_skip.setAttribute('id', "skip")
        g_skip.addElement(line1)
        g_skip.addElement(line2)

        path_defs = Defs()
        path_defs.addElement(g_insertion)
        path_defs.addElement(g_skip)
        return path_defs

    def makePathInsertionsGroup(self) -> G:
        part = self.cn_doc.part
        id_coords = self.mapIdnumsToYcoords()
        _BW = self._base_width
        _hBW = _BW/2
        _BH = self._base_height
        _hBH = _BH/2
        _pX = self.PATH_X_PADDING + self._path_radius_scaled*3

        g = G()
        g.setAttribute('id', "PathInsertions")
        g.setAttribute('font-size', "5")
        g.setAttribute("font-family", "'SourceSansPro-Regular'")
        g.setAttribute("text-anchor", "middle")

        insertions = self.cn_doc.insertions

        for id_num in sorted(insertions.keys(), reverse=True):
            fwd_ins_y = id_coords[id_num] - _BH
            rev_ins_y = id_coords[id_num]
            for ins_idx, ins_len, fwd_col, rev_col in insertions[id_num][::-1]:
                ins_g = G()
                ins_g.setAttribute('id', "Ins %s[%s]" % (id_num, ins_idx))
                ins_x = _BW*ins_idx + _pX
                ins_fwd, ins_rev = part.hasStrandAtIdx(id_num, ins_idx)
                if ins_fwd:
                    u = Use()
                    u.setAttribute('xlink:href', '#insertion')
                    u.setAttribute('x', ins_x)
                    u.setAttribute('y', fwd_ins_y)
                    fwd_style = 'fill:none; stroke:%s; stroke-width:2' % fwd_col
                    u.set_style(fwd_style)
                    ins_g.addElement(u)
                    t = Text('%s' % ins_len, ins_x+_hBW, fwd_ins_y)
                    t.set_style('fill:#999999')
                    ins_g.addElement(t)
                if ins_rev:
                    u = Use()
                    u.setAttribute('xlink:href', '#insertion')
                    u.setAttribute('x', ins_x)
                    u.setAttribute('y', rev_ins_y)
                    rev_style = 'fill:none; stroke:%s; stroke-width:2' % rev_col
                    u.set_style(rev_style)
                    u.set_transform('rotate(180,%s,%s)' % (ins_x+_hBW, rev_ins_y+_hBH))
                    ins_g.addElement(u)
                    t = Text('%s' % ins_len, ins_x+_hBW, rev_ins_y+_BH*1.3)
                    t.set_style('fill:#999999')
                    ins_g.addElement(t)
                g.addElement(ins_g)
        # print(g.getXML())
        return g

    def makePathSkipsGroup(self) -> G:
        """
        Creates and returns a 'G' object for Path VirtualHelix Labels.
        Elements are sorted in reverse order so they are listed in
        ascending order in Illustrator.
        """
        _BW = self._base_width
        _BH = self._base_height
        _pX = self.PATH_X_PADDING + self._path_radius_scaled*3
        id_coords = self.mapIdnumsToYcoords()

        part = self.cn_doc.part

        g = G()
        g.setAttribute('id', "PathSkips")
        skips = self.cn_doc.skips

        for id_num in sorted(skips.keys(), reverse=True):
            fwd_skip_y = id_coords[id_num] - _BH
            rev_skip_y = id_coords[id_num]
            for skip_idx, _ in skips[id_num][::-1]:
                skip_g = G()
                skip_g.setAttribute('id', "Skip %s[%s]" % (id_num, skip_idx))
                skip_x = _BW*skip_idx + _pX
                skip_fwd, skip_rev = part.hasStrandAtIdx(id_num, skip_idx)
                if skip_fwd:
                    u = Use()
                    u.setAttribute('xlink:href', '#skip')
                    u.setAttribute('x', skip_x)
                    u.setAttribute('y', fwd_skip_y)
                    skip_g.addElement(u)
                if skip_rev:
                    u = Use()
                    u.setAttribute('xlink:href', '#skip')
                    u.setAttribute('x', skip_x)
                    u.setAttribute('y', rev_skip_y)
                    skip_g.addElement(u)
                g.addElement(skip_g)
        return g
    # end def

    def makePathSvg(self) -> Svg:
        """
        Returns the main `Svg` object for Path view.

        """
        viewbox = "0 0 %s %s" % (self.w, self.h)
        path_svg = Svg(width=self.w, height=self.h)
        path_svg.set_viewBox(viewbox)
        path_svg.set_preserveAspectRatio("xMinYMid meet")
        path_svg.setAttribute('id', "Cadnano_Path")  # Main layer name
        path_svg.addElement(self.defs)
        if not self._heatmap:
            path_svg.addElement(self.g_pathgridlines)  # bottom layer
        path_svg.addElement(self.g_patholigos)
        path_svg.addElement(self.g_pathstaplexos)
        # if not self._heatmap:
        path_svg.addElement(self.g_pathendpoints)
        path_svg.addElement(self.g_pathvirtualhelices)
        path_svg.addElement(self.g_pathvirtualhelixlabels)
        path_svg.addElement(self.g_pathinsertions)
        path_svg.addElement(self.g_pathskips)  # top layer
        if not self._heatmap and self.cn_doc.sequence_applied:
            path_svg.addElement(self.g_pathsequences)
        # else:
        #     print('No sequences were applied. Max oligo length: %s' % self.cn_doc.max_oligo_length, file=sys.stderr)

        path_svg.save(self.output_path)

        return path_svg
    # end def
# end class


class DefaultArgs(argparse.Namespace):
    input      = None  # Cadnano json file
    output     = None  # Output directory
    seq        = None  # Scaffold sequence file
    heatmap    = False # Hide staple crossover quad curves


def parse_args_from_shell(parser):
    parser.add_argument('--input', '-i', type=str, required=True, nargs=1, metavar='FILE',
                        help='Cadnano JSON file')
    parser.add_argument('--output', '-o', type=str, nargs='?', metavar='DIR',
                        help='Output directory')
    parser.add_argument('--seq', '-s', type=str, nargs='?', metavar='SEQUENCE.txt',
                        help='Scaffold sequence file')
    parser.add_argument('--heatmap', '-H', action='store_true',
                        help='Render compact heatmap-friendly style with no staple xovers.')
    return parser.parse_args()


def run(is_notebook_session=False, args=None):
    if args is None:
        # Get arguments from command line or notebook environment
        parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        args = DefaultArgs() if is_notebook_session else parse_args_from_shell(parser)
        if args.input is None:
            parser.print_help()
            sys.exit('Input file not specified')
        design = args.input[0]
        if not design.endswith('.json'):
            parser.print_help()
            sys.exit('Input should be JSON file')
    else:
        # Assume args included in run function call
        design = args.input

    output_directory = args.output

    sequence = ''
    if args.seq is not None:
        valid = 'ACTGactg'
        with open(args.seq) as seqfile:
            sequence = ''.join(seqfile.read().split())
            if all(s in valid for s in sequence) and not is_notebook_session:
                print('Found valid sequence file, %s bases' % len(sequence), file=sys.stderr)
            else:
                sys.exit('Error: %s contains non-[ATCGactg] character(s)' % args.seq)

    basename = os.path.splitext(os.path.basename(design))[0]
    base_path = os.path.splitext(design)[0]
    cndoc = CadnanoDocument(design, sequence)

    # File Path
    if output_directory:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory, exist_ok=True)
        output_ortho = os.path.join(output_directory, '%s_ortho' % basename)
        output_path = os.path.join(output_directory, '%s_path' % basename)
    else:
        output_ortho = '%s_ortho' % base_path
        output_path = '%s_path' % base_path

    output_ortho += '.svg'
    output_path += '.svg'

    ortho_svg = CadnanoOrthoSvg(cndoc, output_ortho)
    path_svg = CadnanoPathSvg(cndoc, output_path, heatmap=args.heatmap)
    dimensions = {'ortho_svg_width': ortho_svg.w,
                  'ortho_svg_height': ortho_svg.h,
                  'path_svg_width': path_svg.w,
                  'path_svg_height': path_svg.h}
    # print(json.dumps(dimensions))

def main():
    try:
        __IPYTHON__
        is_notebook_session = True
    except NameError:
        is_notebook_session = False

    run(is_notebook_session=is_notebook_session)