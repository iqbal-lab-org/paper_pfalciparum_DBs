from abc import ABC
import itertools as it
import numpy as np
from typing import List
from dataclasses import dataclass

SVG_HEADER = (
    '<?xml version="1.0" standalone="no"?>\n'
    '<!DOCTYPE svg PUBLIC " -//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n'
)

DEFAULT_RECT_WIDTH = 15
DEFAULT_RECT_HEIGHT = 20

DEFAULT_OPACITY = 1
DEFAULT_TEXT_FILL = "black"


class BasicSVG(ABC):
    def add_properties(self, **fields):
        for pname, pval in fields.items():
            if pname not in self.VALID_PROPERTIES:
                raise ValueError(f"{pname} not in {self.VALID_PROPERTIES}")
            setattr(self, pname, pval)

    def serialise_properties(self):
        result = ""
        for pname in self.VALID_PROPERTIES:
            pval = getattr(self, pname, None)
            if pval is not None:
                used_pname = pname.replace("_", "-")
                result += f'{used_pname}="{pval}" '
        return result


class Coord:
    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y

    def __repr__(self):
        return f"{self.x},{self.y}"

    def __add__(self, other: "Coord"):
        return Coord(self.x + other.x, self.y + other.y)


class Text(BasicSVG):
    VALID_PROPERTIES = {
        "x",
        "y",
        "dominant-baseline",
        "text-anchor",
        "opacity",
        "fill",
        "font_style",
    }

    def __init__(self, text: str, **properties):
        self.text = text
        self.add_properties(**properties)

    def __repr__(self):
        result = "<text "
        result += self.serialise_properties()
        result += f">{self.text}</text>"
        return result


class Line(BasicSVG):
    VALID_PROPERTIES = {"stroke", "stroke_width"}

    def __init__(self, start: Coord, end: Coord):
        self.start = start
        self.end = end

    def __repr__(self):
        result = f'<line x1="{self.start.x}" y1="{self.start.y}" '
        result += f'x2="{self.end.x}" y2="{self.end.y}" '
        result += self.serialise_properties()
        result = result.rstrip()
        result += "/>"
        return result


@dataclass
class LineSpec:
    index: int
    colour: str
    width: float


LineSpecs = List[LineSpec]


class Grid:
    def __init__(
        self,
        topleft: Coord,
        nrow: int,
        ncol: int,
        grid_spacing: float = 0.0,
        rect_width: float = DEFAULT_RECT_WIDTH,
        rect_height: float = DEFAULT_RECT_HEIGHT,
        vlines: LineSpecs = {},
        hlines: LineSpecs = {},
    ):
        """
        `spacing`: amount of whitespace to add between elements of the grid
        """
        self._size_x = ncol * (rect_width + grid_spacing) - grid_spacing
        self._size_x += sum(map(lambda elem: elem.width, vlines))
        self._size_y = nrow * (rect_height + grid_spacing) - grid_spacing
        self._size_y += sum(map(lambda elem: elem.width, hlines))
        self.rectangles = np.empty((nrow, ncol), dtype=object)
        self.topleft = topleft
        self._nrow = nrow
        self._ncol = ncol
        self.draw_rects(rect_width, rect_height, grid_spacing, vlines, hlines)
        self.draw_lines(vlines, hlines)

    def draw_rects(self, rect_width, rect_height, grid_spacing, vlines, hlines):
        used_vlines = {elem.index: elem for elem in vlines}
        used_hlines = {elem.index: elem for elem in hlines}
        x_offset, y_offset = 0, 0
        for row_idx in range(self._nrow):
            for col_idx in range(self._ncol):
                new_Coord = Coord(x_offset, y_offset)
                self.rectangles[row_idx][col_idx] = Rectangle(
                    new_Coord, rect_width, rect_height
                )
                x_offset += rect_width + grid_spacing
                if col_idx in used_vlines:
                    x_offset += used_vlines[col_idx].width
            y_offset += rect_height + grid_spacing
            if row_idx in used_hlines:
                y_offset += used_hlines[row_idx].width
            x_offset = 0

    def draw_lines(self, vlines, hlines):
        self._lines = list()
        # Note: i want vlines to take precedence over hlines if they overlap,
        # so make hlines first.
        for hline in hlines:
            left_rectangle = self.rectangles[hline.index][0]
            left_coord = Coord(
                left_rectangle._bottomleft.x,
                left_rectangle._bottomleft.y + hline.width / 2,
            )
            right_rectangle = self.rectangles[hline.index][self._ncol - 1]
            right_coord = Coord(
                right_rectangle._bottomright.x,
                right_rectangle._bottomright.y + hline.width / 2,
            )
            new_hline = Line(left_coord, right_coord)
            new_hline.add_properties(stroke=hline.colour, stroke_width=hline.width)
            self._lines.append(new_hline)
        for vline in vlines:
            top_rectangle = self.rectangles[0][vline.index]
            top_coord = Coord(
                top_rectangle._topright.x + vline.width / 2, top_rectangle._topright.y
            )
            bottom_rectangle = self.rectangles[self._nrow - 1][vline.index]
            bottom_coord = Coord(
                bottom_rectangle._bottomright.x + vline.width / 2,
                bottom_rectangle._bottomright.y,
            )
            new_vline = Line(top_coord, bottom_coord)
            new_vline.add_properties(stroke=vline.colour, stroke_width=vline.width)
            self._lines.append(new_vline)

    @property
    def size(self):
        return (self._size_x, self._size_y)

    def __getitem__(self, row_idx: int):
        return self.rectangles[row_idx]

    def __repr__(self):
        result = ""
        for elem in it.chain.from_iterable(self.rectangles):
            result += f"{str(elem)}\n"
        for elem in self._lines:
            result += f"{str(elem)}\n"
        return result


class Rectangle(BasicSVG):
    VALID_PROPERTIES = {"stroke", "stroke_width", "fill", "fill_opacity"}

    def __init__(
        self,
        topleft: Coord,
        width: float = DEFAULT_RECT_WIDTH,
        height: float = DEFAULT_RECT_HEIGHT,
    ):
        self.topleft = topleft
        self.width = width
        self.height = height
        self._topright = self.topleft + Coord(self.width, 0)
        self._bottomleft = self.topleft + Coord(0, self.height)
        self._bottomright = self.topleft + Coord(self.width, self.height)
        self._text = None

    def add_text(self, text: str, opacity=DEFAULT_OPACITY, fill=DEFAULT_TEXT_FILL):
        text_properties = {
            "opacity": opacity,
            "fill": fill,
            "x": str(self.topleft.x + self.width / 2),
            "y": str(self.topleft.y + self.height / 2),
            "text-anchor": "middle",
            "dominant-baseline": "middle",
        }
        self._text = Text(text, **text_properties)

    def add_text_properties(self, **properties):
        self._text.add_properties(**properties)

    def __repr__(self):
        result = f'<rect x="{self.topleft.x}" y="{self.topleft.y}" '
        result += f'width="{self.width}" height="{self.height}" '
        result += self.serialise_properties()
        result = result.rstrip()
        result += "/>"
        if self._text is not None:
            result += str(self._text)
        return result
