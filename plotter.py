__author__ = 'github.com/wardsimon'
__version__ = '0.0.1'

import numpy as np
import vtk
import pyvista

from easyCore.Symmetry.Bonding import generate_bonds

import random


def get_random_color(pastel_factor=0.5):
    return [(x + pastel_factor) / (1.0 + pastel_factor) for x in [random.uniform(0, 1.0) for i in [1, 2, 3]]]


def color_distance(c1, c2):
    return sum([abs(x[0] - x[1]) for x in zip(c1, c2)])


def generate_new_color(existing_colors, pastel_factor=0.5):
    max_distance = None
    best_color = None
    for i in range(0, 100):
        color = get_random_color(pastel_factor=pastel_factor)
        if not existing_colors:
            return color
        best_distance = min([color_distance(color, c) for c in existing_colors])
        if not max_distance or best_distance > max_distance:
            max_distance = best_distance
            best_color = color
    return best_color


class plotter:
    def __init__(self, phaseObj, extent=[1, 1, 1], max_distance=3):

        self.canvas = None
        self.componentsShow = {'lattice': [],
                               'atoms':   [],
                               'bonds':   []
                               }

        self._max_distance = max_distance
        self._phaseObj = phaseObj
        self._extent = extent

        self.update()

    @property
    def phaseObj(self):
        return self._phaseObj

    @phaseObj.setter
    def phaseObj(self, value):
        self._phaseObj = value
        self.update()

    @property
    def max_distance(self):
        return self._max_distance

    @max_distance.setter
    def max_distance(self, value):
        self._max_distance = np.array(value)
        self.generate_bonds()

    @property
    def extent(self):
        return self._extent

    @extent.setter
    def extent(self, value):
        self._extent = value
        self.update()

    def plot(self):
        if self.canvas is not None:
            actors = [*self.componentsShow['lattice'],
                      *self.componentsShow['bonds'],
                      *self.componentsShow['atoms']]
            for actor in actors:
                self.canvas.add_actor(actor, render=False)
            self.canvas.reset_camera()

    def update(self):
        self.generate_lattice()
        self.generate_bonds()
        self.generate_atoms()

    def generate_bonds(self):

        bonds = generate_bonds(self.phaseObj, max_distance=self.max_distance)
        all_atoms = self.phaseObj.get_orbits(magnetic_only=False)
        all_atoms_r = np.vstack([np.array(all_atoms[key]) for key in all_atoms.keys()])

        # generate all cell translations
        cTr1, cTr2, cTr3 = np.mgrid[0:(self.extent[0] + 1), 0:(self.extent[1] + 1), 0:(self.extent[2] + 1)]
        # cell origin translations: Na x Nb x Nc x 1 x 1 x3
        pos = np.stack([cTr1, cTr2, cTr3], axis=0).reshape((3, -1, 1), order='F')
        R1 = (pos + all_atoms_r[bonds.atom1, :].T.reshape((3, -1, 1)).transpose((0, 2, 1))).reshape((3, -1), order='F')
        R2 = (pos + (all_atoms_r[bonds.atom2, :].T + bonds.dl).reshape((3, -1, 1)).transpose((0, 2, 1))).reshape(
            (3, -1), order='F')
        IDX = np.tile(bonds.idx, (pos.shape[1]))

        lattice = self.phaseObj.cell
        pos = []
        radius = 0.05
        colors = []
        for _ in range(int(bonds.nSym)):
            colors.append(generate_new_color(colors, pastel_factor=0.9))

        def rotation_matrix_from_vectors(vec1, vec2):
            """ Find the rotation matrix that aligns vec1 to vec2
            :param vec1: A 3d "source" vector
            :param vec2: A 3d "destination" vector
            :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
            """
            a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
            v = np.cross(a, b)
            c = np.dot(a, b)
            s = np.linalg.norm(v)
            kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
            if np.all(kmat == 0):
                return np.eye(3)
            rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
            return rotation_matrix

        for r1, r2, idx in zip(R1.T, R2.T, IDX):
            if np.any(r1 < 0) | np.any(r2 < 0) | \
                    np.any(r1 > self.extent) | np.any(r2 > self.extent):
                # In this silly example we only think about extent = [1, 1, 1]
                continue

            sp = lattice.get_cartesian_coords(r1)
            ep = lattice.get_cartesian_coords(r2)
            cp = (ep - sp) / 2
            polyDataSource = vtk.vtkCylinderSource()
            polyDataSource.SetRadius(radius)
            length = np.linalg.norm(ep - sp)
            polyDataSource.SetHeight(length)

            transform = vtk.vtkTransform()
            rot = rotation_matrix_from_vectors([0, 1, 0], cp)
            np.linalg.norm(sp - ep)
            wxyz = np.identity(4)
            wxyz[0:3, 0:3] = rot
            wxyz[0:3, 3] = sp + cp

            m = vtk.vtkMatrix4x4()
            m.DeepCopy(wxyz.ravel())
            transform.SetMatrix(m)
            transformPD = vtk.vtkTransformPolyDataFilter()
            transformPD.SetTransform(transform)
            transformPD.SetInputConnection(polyDataSource.GetOutputPort())

            transform2 = vtk.vtkTransform()
            wxyz2 = np.identity(4)
            wxyz2[0:3, 3] = - lattice.get_cartesian_coords(0.5*np.array(self.extent))
            m2 = vtk.vtkMatrix4x4()
            m2.DeepCopy(wxyz2.ravel())
            transform2.SetMatrix(m2)
            transformPD2 = vtk.vtkTransformPolyDataFilter()
            transformPD2.SetTransform(transform2)
            transformPD2.SetInputConnection(transformPD.GetOutputPort())

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(transformPD2.GetOutputPort())
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(*colors[idx])
            pos.append(actor)
        self.componentsShow['bonds'] = pos

    def generate_lattice(self):

        lattice = self.phaseObj.cell
        actors = []
        radius = 0.035

        def rotation_matrix_from_vectors(vec1, vec2):
            """ Find the rotation matrix that aligns vec1 to vec2
            :param vec1: A 3d "source" vector
            :param vec2: A 3d "destination" vector
            :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
            """
            a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
            v = np.cross(a, b)
            c = np.dot(a, b)
            s = np.linalg.norm(v)
            kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
            if np.all(kmat == 0):
                return np.eye(3)
            rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
            return rotation_matrix

        def create_cylinder(edge):
            s_in, e_in = edge
            sp = lattice.get_cartesian_coords(s_in)
            ep = lattice.get_cartesian_coords(e_in)
            cp = (ep - sp) / 2
            polyDataSource = vtk.vtkCylinderSource()
            polyDataSource.SetRadius(radius)
            length = np.linalg.norm(ep - sp)
            polyDataSource.SetHeight(length)

            transform = vtk.vtkTransform()
            rot = rotation_matrix_from_vectors([0, 1, 0], cp)
            np.linalg.norm(sp - ep)
            wxyz = np.identity(4)
            wxyz[0:3, 0:3] = rot
            wxyz[0:3, 3] = sp + cp
            # wxyz[2, 3] = wxyz[2, 3] + 10
            # transform.Concatenate(*wxyx.reshape(1, -1).tolist())
            m = vtk.vtkMatrix4x4()
            m.DeepCopy(wxyz.ravel())
            transform.SetMatrix(m)
            transformPD = vtk.vtkTransformPolyDataFilter()
            transformPD.SetTransform(transform)
            transformPD.SetInputConnection(polyDataSource.GetOutputPort())

            transform2 = vtk.vtkTransform()
            wxyz2 = np.identity(4)
            wxyz2[0:3, 3] = - lattice.get_cartesian_coords(0.5*np.array(self.extent))
            m2 = vtk.vtkMatrix4x4()
            m2.DeepCopy(wxyz2.ravel())
            transform2.SetMatrix(m2)
            transformPD2 = vtk.vtkTransformPolyDataFilter()
            transformPD2.SetTransform(transform2)
            transformPD2.SetInputConnection(transformPD.GetOutputPort())

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(transformPD2.GetOutputPort())
            return mapper

        extent = self.extent
        lines = []

        # Add major lattice points
        for z in np.arange(extent[2]):
            for y in np.arange(extent[1]):
                for x in np.arange(extent[0]):
                    lines.append(create_cylinder(([x, y, z], [x + 1, y, z])))
                    lines.append(create_cylinder(([x, y, z], [x, y + 1, z])))
                    lines.append(create_cylinder(([x, y, z], [x, y, z + 1])))

                    lines.append(create_cylinder(([x + 1, y + 1, z + 1], [x + 1, y, z + 1])))
                    lines.append(create_cylinder(([x + 1, y + 1, z + 1], [x + 1, y + 1, z])))
                    lines.append(create_cylinder(([x + 1, y + 1, z + 1], [x, y + 1, z + 1])))

        # Draw x lines
        for x in np.arange(extent[0]):
            y = extent[1]
            z = 0
            lines.append(create_cylinder(([x, y, z], [x + 1, y, z])))
            z = extent[2]
            y = 0
            lines.append(create_cylinder(([x, y, z], [x + 1, y, z])))

        # Draw y lines
        for y in np.arange(extent[1]):
            x = 0
            z = extent[2]
            lines.append(create_cylinder(([x, y, z], [x, y + 1, z])))
            x = extent[0]
            z = 0
            lines.append(create_cylinder(([x, y, z], [x, y + 1, z])))

        # Draw y lines
        for z in np.arange(extent[2]):
            x = 0
            y = extent[1]
            lines.append(create_cylinder(([x, y, z], [x, y, z + 1])))
            x = extent[0]
            y = 0
            lines.append(create_cylinder(([x, y, z], [x, y, z + 1])))

        for mapp in lines:
            actor = vtk.vtkActor()
            actor.SetMapper(mapp)
            actors.append(actor)
        
        self.componentsShow['lattice'] = actors

    def generate_atoms(self):
        atoms = self.phaseObj.all_sites(extent=np.array(self.extent))
        lattice = self.phaseObj.cell
        actors = []
        colors = ["#ff7f50", "#4682b4", "#6b8e23", "#d2691e", "#5f9ea0", "#8fbc8f", "#6495ed"]
        radii = [(0.1 + 0.004 * i) * min(lattice.lengths) for i in range(len(colors))]
        res = 50

        for atom_index, atom_name in enumerate(atoms.keys()):
            radius = radii[atom_index]
            color = hex_to_rgb(colors[atom_index])
            for atom in atoms[atom_name]:
                polyDataSource = vtk.vtkSphereSource()
                polyDataSource.SetRadius(radius)
                polyDataSource.SetPhiResolution(res)
                polyDataSource.SetThetaResolution(res)
                polyDataSource.SetCenter(*lattice.get_cartesian_coords(atom))

                transform = vtk.vtkTransform()
                wxyz = np.identity(4)
                wxyz[0:3, 3] = - lattice.get_cartesian_coords(0.5*np.array(self.extent))
                m = vtk.vtkMatrix4x4()
                m.DeepCopy(wxyz.ravel())
                transform.SetMatrix(m)
                transformPD = vtk.vtkTransformPolyDataFilter()
                transformPD.SetTransform(transform)
                transformPD.SetInputConnection(polyDataSource.GetOutputPort())

                mapper = vtk.vtkPolyDataMapper()
                mapper.SetInputConnection(transformPD.GetOutputPort())
                actor = vtk.vtkActor()
                actor.SetMapper(mapper)
                actor.GetProperty().SetColor(color)
                actor.GetProperty().SetAmbient(0.35)
                actor.GetProperty().SetDiffuse(0.55)
                actor.GetProperty().SetSpecular(0.5)
                actor.GetProperty().SetSpecularPower(7.0)
                actors.append(actor)
        self.componentsShow['atoms'] = actors

    def clearScene(self):
        actors = [*self.componentsShow['lattice'],
                  *self.componentsShow['bonds'],
                  *self.componentsShow['atoms']]
        if self.canvas is not None:
            self.canvas.remove_actor(actors, render=False)
            self.canvas.update()


class QMLPlotter(plotter):
    def __init__(self, phaseObj, **kwargs):
        super(QMLPlotter, self).__init__(phaseObj)
        self.canvas = pyvista.Plotter(**kwargs)

    def plot(self):
        super(winPlotter, self).plot()
        self.canvas.update()


class winPlotter(plotter):
    def __init__(self, phaseObj, **kwargs):
        super(winPlotter, self).__init__(phaseObj)
        self.canvas = pyvista.Plotter(**kwargs)

    def plot(self):
        super(winPlotter, self).plot()
        self.canvas.show()


class nbPlotter(winPlotter):
    def __init__(self, phaseObj):
        super(nbPlotter, self).__init__(phaseObj, notebook=True)


class nbInteractPlotter(plotter):

    def __init__(self, phaseObj):
        super(nbInteractPlotter, self).__init__(phaseObj)
        self.canvas = pyvista.PlotterITK()

    def plot(self):
        super(nbInteractPlotter, self).plot()
        self.canvas.show()


# class qtWindow(Qt.QMainWindow, plotter):
#     def __init__(self, phaseObj, parent=None, show=True):
#         Qt.QMainWindow.__init__(self, parent)
#         plotter.__init__(self, phaseObj)
#
#         # create the frame
#         self.frame = Qt.QFrame()
#         vlayout = Qt.QVBoxLayout()
#
#         # Button Layout
#         btn_layout = Qt.QHBoxLayout()
#         btn_layout.addItem(Qt.QSpacerItem(0, 0, Qt.QSizePolicy.Expanding, Qt.QSizePolicy.Minimum))
#         btn1 = Qt.QPushButton("Button 1")
#         # btn1.clicked.connect(self.viewA)
#         btn2 = Qt.QPushButton("Button 2")
#         btn_layout.addWidget(btn1)
#         btn_layout.addWidget(btn2)
#
#         # add the vtki interactor object
#         self.vtk_widget = pyvista.QtInteractor(self.frame)
#         vlayout.addWidget(self.vtk_widget)
#         vlayout.addLayout(btn_layout)
#
#         self.frame.setLayout(vlayout)
#         self.setCentralWidget(self.frame)
#
#         # simple menu to demo functions
#         mainMenu = self.menuBar()
#         fileMenu = mainMenu.addMenu('File')
#         exitButton = Qt.QAction('Exit', self)
#         exitButton.setShortcut('Ctrl+Q')
#         exitButton.triggered.connect(self.close)
#         fileMenu.addAction(exitButton)
#
#         self.plotter = plotter(spinWobj, self.vtk_widget)
#
#         # allow adding a sphere
#         meshMenu = mainMenu.addMenu('Add Elements')
#         self.add_lattice_action = Qt.QAction('Lattice', self)
#         self.add_lattice_action.triggered.connect(self.add_lattice)
#         meshMenu.addAction(self.add_lattice_action)
#         self.add_atoms_action = Qt.QAction('Atoms', self)
#         self.add_atoms_action.triggered.connect(self.add_atoms)
#         meshMenu.addAction(self.add_atoms_action)
#         self.add_bonds_action = Qt.QAction('Bonds', self)
#         self.add_bonds_action.triggered.connect(self.add_bonds)
#         meshMenu.addAction(self.add_bonds_action)
#
#         RMmeshMenu = mainMenu.addMenu('Remove Elements')
#         self.rm_lattice_action = Qt.QAction('Lattice', self)
#         self.rm_lattice_action.triggered.connect(self.rm_lattice)
#         RMmeshMenu.addAction(self.rm_lattice_action)
#         self.rm_atoms_action = Qt.QAction('Atoms', self)
#         self.rm_atoms_action.triggered.connect(self.rm_atoms)
#         RMmeshMenu.addAction(self.rm_atoms_action)
#         self.rm_bonds_action = Qt.QAction('Bonds', self)
#         self.rm_bonds_action.triggered.connect(self.rm_bonds)
#         RMmeshMenu.addAction(self.rm_bonds_action)
#
#         self.vtk_widget.background_color = 'white'
#         if show:
#             self.show()

def hex_to_rgb(hex):
    hex = hex.lstrip('#')
    rgb = [int(hex[i:i + 2], 16) for i in (0, 2, 4)]
    rgb = [i / 255 for i in rgb]
    return rgb
