__author__ = 'github.com/wardsimon'
__version__ = '0.0.1'

import numpy as np
import pyvista
from .BasePlotter import BasePlotter
from easyCore.Elements.Basic.Lattice import Lattice as coreLattice


class Lattice(BasePlotter):

    def __init__(self, plotter, obj, *args, **kwargs):
        if not isinstance(obj, coreLattice):
            raise AttributeError

        super(Lattice, self).__init__(plotter, obj, *args, **kwargs)

    def plot(self, *args, extent=None, color=(0.45, 0.45, 0.45), centre=(0, 0, 0), **kwargs):

        if extent is None:
            extent = [1, 1, 1]

        lattice: coreLattice = self.obj
        lines = []

        # Add major lattice points
        for z in np.arange(extent[2]):
            for y in np.arange(extent[1]):
                for x in np.arange(extent[0]):
                    lines.append(pyvista.Line(pointa=lattice.get_cartesian_coords([x, y, z]), pointb=lattice.get_cartesian_coords([x + 1, y, z])))
                    lines.append(pyvista.Line(pointa=lattice.get_cartesian_coords([x, y, z]), pointb=lattice.get_cartesian_coords([x, y + 1, z])))
                    lines.append(pyvista.Line(pointa=lattice.get_cartesian_coords([x, y, z]), pointb=lattice.get_cartesian_coords([x, y, z + 1])))

                    lines.append(pyvista.Line(pointa=lattice.get_cartesian_coords([x + 1, y + 1, z + 1]),
                                           pointb=lattice.get_cartesian_coords([x + 1, y, z + 1])))
                    lines.append(pyvista.Line(pointa=lattice.get_cartesian_coords([x + 1, y + 1, z + 1]),
                                           pointb=lattice.get_cartesian_coords([x + 1, y + 1, z])))
                    lines.append(pyvista.Line(pointa=lattice.get_cartesian_coords([x + 1, y + 1, z + 1]),
                                           pointb=lattice.get_cartesian_coords([x, y + 1, z + 1])))

        # Draw x lines
        for x in np.arange(extent[0]):
            y = extent[1]
            z = 0
            lines.append(pyvista.Line(pointa=lattice.get_cartesian_coords([x, y, z]), pointb=lattice.get_cartesian_coords([x + 1, y, z])))
            z = extent[2]
            y = 0
            lines.append(pyvista.Line(pointa=lattice.get_cartesian_coords([x, y, z]), pointb=lattice.get_cartesian_coords([x + 1, y, z])))

        # Draw y lines
        for y in np.arange(extent[1]):
            x = 0
            z = extent[2]
            lines.append(pyvista.Line(pointa=lattice.get_cartesian_coords([x, y, z]), pointb=lattice.get_cartesian_coords([x, y + 1, z])))
            x = extent[0]
            z = 0
            lines.append(pyvista.Line(pointa=lattice.get_cartesian_coords([x, y, z]), pointb=lattice.get_cartesian_coords([x, y + 1, z])))

        # Draw y lines
        for z in np.arange(extent[2]):
            x = 0
            y = extent[1]
            lines.append(pyvista.Line(pointa=lattice.get_cartesian_coords([x, y, z]), pointb=lattice.get_cartesian_coords([x, y, z + 1])))
            x = extent[0]
            y = 0
            lines.append(pyvista.Line(pointa=lattice.get_cartesian_coords([x, y, z]), pointb=lattice.get_cartesian_coords([x, y, z + 1])))


        centre_point = lattice.get_cartesian_coords(extent)/2
        for line in lines:
            line.transform(-centre_point)
            self.actors.append(self.plotter.add_mesh(line, color=color))
        self.plotter.reset_camera()
