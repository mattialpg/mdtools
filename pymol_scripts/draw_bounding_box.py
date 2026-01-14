# -*- coding: utf-8 -*-
from pymol.cgo import BEGIN, LINES, VERTEX, LINEWIDTH, END, COLOR
from pymol import cmd

def draw_bounding_box(selection="(sele)", padding=0.0, linewidth=2.0):
    padding = float(padding)
    linewidth = float(linewidth)

    (minX0, minY0, minZ0), (maxX0, maxY0, maxZ0) = cmd.get_extent(selection)

    # unpadded center & size
    dx0 = maxX0 - minX0
    dy0 = maxY0 - minY0
    dz0 = maxZ0 - minZ0
    Ox = 0.5 * (maxX0 + minX0)
    Oy = 0.5 * (maxY0 + minY0)
    Oz = 0.5 * (maxZ0 + minZ0)

    # apply padding
    minX = minX0 - padding; minY = minY0 - padding; minZ = minZ0 - padding
    maxX = maxX0 + padding; maxY = maxY0 + padding; maxZ = maxZ0 + padding

    dx = maxX - minX; dy = maxY - minY; dz = maxZ - minZ

    print("Box center (%.2f, %.2f, %.2f)" % (Ox, Oy, Oz))
    print("Box dimensions (%.2f, %.2f, %.2f)" % (dx0, dy0, dz0))
    if padding != 0:
        print("Box dimensions + padding (%.2f, %.2f, %.2f)" % (dx, dy, dz))

    box = [
        LINEWIDTH, linewidth,
        BEGIN, LINES,
        COLOR, 1.0, 1.0, 1.0,

        # 4 vertical edges
        VERTEX, minX, minY, minZ, VERTEX, minX, minY, maxZ,
        VERTEX, maxX, minY, minZ, VERTEX, maxX, minY, maxZ,
        VERTEX, maxX, maxY, minZ, VERTEX, maxX, maxY, maxZ,
        VERTEX, minX, maxY, minZ, VERTEX, minX, maxY, maxZ,

        # 4 bottom edges (z = minZ)
        VERTEX, minX, minY, minZ, VERTEX, maxX, minY, minZ,
        VERTEX, maxX, minY, minZ, VERTEX, maxX, maxY, minZ,
        VERTEX, maxX, maxY, minZ, VERTEX, minX, maxY, minZ,
        VERTEX, minX, maxY, minZ, VERTEX, minX, minY, minZ,

        # 4 top edges (z = maxZ)
        VERTEX, minX, minY, maxZ, VERTEX, maxX, minY, maxZ,
        VERTEX, maxX, minY, maxZ, VERTEX, maxX, maxY, maxZ,
        VERTEX, maxX, maxY, maxZ, VERTEX, minX, maxY, maxZ,
        VERTEX, minX, maxY, maxZ, VERTEX, minX, minY, maxZ,

        END
    ]

    cmd.load_cgo(box, "box")
    return "box"

cmd.extend("draw_bounding_box", draw_bounding_box)
