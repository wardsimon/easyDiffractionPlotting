# __author__ = 'github.com/wardsimon'
# __version__ = '0.0.1'
#
# from typing import Dict
#
# import numpy as np
# import pyvista
# import vtk
# from .BasePlotter import BasePlotter
# from easyCore.Elements.Basic.Site import PeriodicAtoms
#
#
# def all_orbits(atoms, extent, atom_tolerance, centre) -> Dict[str, np.ndarray]:
#     """
#     Generate all atomic positions from the atom array and symmetry operations over an extent.
#
#     :return:  dictionary with keys of atom labels, containing numpy arrays of unique points in the extent
#     (0, 0, 0) -> obj.extent
#     :rtype: Dict[str, np.ndarray]
#     """
#
#     offsets = np.array(np.meshgrid(range(0, extent[0] + 1),
#                                    range(0, extent[1] + 1),
#                                    range(0, extent[2] + 1))).T.reshape(-1, 3)
#
#     orbits = atoms.get_orbits()
#     for orbit_key in orbits.keys():
#         orbit = orbits[orbit_key]
#         site_positions = np.apply_along_axis(np.add, 1, offsets, orbit).reshape((-1, 3)) - centre
#         orbits[orbit_key] = \
#             site_positions[np.all(site_positions >= - atom_tolerance, axis=1) &
#                            np.all(site_positions <= extent + atom_tolerance, axis=1),
#             :] + centre
#     return orbits
#
#
# class Atoms(BasePlotter):
#
#     __slots__ = ['extent', 'atom_tolerance', 'radius', 'color']
#
#     def __init__(self, plotter, obj, extent=(1, 1, 1), atom_tolerance=0.01, radius=0.5, color=(1, 1, 1)):
#         super(Atoms, self).__init__(plotter=plotter, obj=obj)
#         if isinstance(radius, float):
#             radius = {atom.lable.raw_value: radius for atom in obj}
#         if isinstance(radius, list):
#             radius = {atom.lable.raw_value: rad for atom, rad in zip(obj, radius)}
#         self.radius = radius
#
#
#
#     def plot(self, *args, **kwargs):
#
#
#         orbits = all_orbits(self.obj, extent, atom_tolerance)
#
#         for key, atom in orbits.items():
#
#
#
#
#
