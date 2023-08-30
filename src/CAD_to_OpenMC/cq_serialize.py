"""
MIT License

Copyright (c) 2022 Simon Dibbern

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


import copyreg
from io import BytesIO

import cadquery as cq
import OCP


def _inflate_shape(data: bytes):
    with BytesIO(data) as bio:
        return cq.Shape.importBrep(bio)


def _reduce_shape(shape: cq.Shape):
    with BytesIO() as stream:
        shape.exportBrep(stream)
        return _inflate_shape, (stream.getvalue(),)


def _inflate_transform(*values: float):
    trsf = OCP.gp.gp_Trsf()
    trsf.SetValues(*values)
    return trsf


def _reduce_transform(transform: OCP.gp.gp_Trsf):
    return _inflate_transform, tuple(
        transform.Value(i, j) for i in range(1, 4) for j in range(1, 5)
    )


def register():
    """
    Registers pickle support functions for common CadQuery and OCCT objects.
    """

    for cls in (
        cq.Edge,
        cq.Compound,
        cq.Shell,
        cq.Face,
        cq.Solid,
        cq.Vertex,
        cq.Wire,
    ):
        copyreg.pickle(cls, _reduce_shape)

    copyreg.pickle(cq.Vector, lambda vec: (cq.Vector, vec.toTuple()))
    copyreg.pickle(OCP.gp.gp_Trsf, _reduce_transform)
    copyreg.pickle(
        cq.Location, lambda loc: (cq.Location, (loc.wrapped.Transformation(),))
    )
