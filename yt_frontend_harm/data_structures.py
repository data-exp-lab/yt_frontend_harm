"""Skeleton data structures."""

# -----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import os
import numpy as np
import weakref

from yt.geometry.unstructured_mesh_handler import UnstructuredIndex
from yt.data_objects.index_subobjects.unstructured_mesh import SemiStructuredMesh
from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.utilities.file_handler import HDF5FileHandler
from itertools import chain, product
from .fields import IHarmFieldInfo


class IHarmGrid(AMRGridPatch):

    _id_offset = 0

    def __init__(self, id, index, level):
        AMRGridPatch.__init__(self, id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        LE, RE = self.index.grid_left_edge[id, :], self.index.grid_right_edge[id, :]
        self.dds = self.ds.arr((RE - LE) / self.ActiveDimensions, "code_length")
        if self.ds.dimensionality < 2:
            self.dds[1] = 1.0
        if self.ds.dimensionality < 3:
            self.dds[2] = 1.0
        self.field_data["dx"], self.field_data["dy"], self.field_data["dz"] = self.dds

    def __repr__(self):
        return "IHarmGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


class IHarmMKSMesh(SemiStructuredMesh):

    _index_offset = 0

    def __init__(self, mesh_id, filename, conn, coords, index, blocks, dims):
        SemiStructuredMesh.__init__(self, mesh_id, filename, conn, coords, index)
        self.mesh_blocks = blocks
        self.mesh_dims = dims


class IHarmMKSHierarchy(UnstructuredIndex):

    grid = IHarmGrid
    _dataset_type = "iharm"

    def __init__(self, ds, dataset_type=_dataset_type):
        self.dataset = weakref.proxy(ds)
        self.directory = os.path.dirname(self.dataset.filename)
        self.dataset_type = dataset_type
        self.index_filename = self.dataset.filename
        self._handle = ds._handle
        UnstructuredIndex.__init__(self, ds, dataset_type)

    def _initialize_mesh(self):
        X1, X2, X3 = self.dataset._get_boundaries()
        X1 = np.linspace(X1[0], X1[-1], self.dataset._x1_info[2] + 1)
        X2 = np.linspace(X2[0], X2[-1], self.dataset._x2_info[2] + 1)
        X3 = np.linspace(X3[0], X3[-1], self.dataset._x3_info[2] + 1)
        X1m, X2m, X3m = self.dataset._map_coordinates(X1, X2, X3)
        self.meshes = []
        cis = np.fromiter(
            chain.from_iterable(product([0, 1], [0, 1], [0, 1])),
            dtype=np.int64,
            count=8 * 3,
        )
        cis.shape = (8, 3)
        nx1 = X1m.size
        nx2 = X2m.size
        nx3 = X3m.size
        coords = np.zeros((nx1, nx2, nx3, 3), dtype="float64", order="C")
        coords[:, :, :, 0] = X1m[:, None, None]
        coords[:, :, :, 1] = X2m[None, :, None]
        coords[:, :, :, 2] = X3m[None, None, :]
        coords.shape = (nx1 * nx2 * nx3, 3)
        cycle = np.rollaxis(np.indices((nx1 - 1, nx2 - 1, nx3 - 1)), 0, 4)
        cycle.shape = ((nx1 - 1) * (nx2 - 1) * (nx3 - 1), 3)
        off = cis + cycle[:, np.newaxis]
        conn = ((off[:, :, 0] * nx2) + off[:, :, 1]) * nx2 + off[:, :, 2]
        mesh = IHarmMKSMesh(
            0,
            self.index_filename,
            conn,
            coords,
            self,
            np.array(1),
            np.array([nx1 - 1, nx2 - 1, nx3 - 1]),
        )
        self.meshes.append(mesh)

    def _detect_output_fields(self):
        self.field_list = [("iharm", s) for s in self.dataset._field_map]

    def _count_grids(self):
        """
        This sets the number of different grids we expect. Since the iharm
        format is not AMR, we can safely assume one here.
        """
        self.num_grids = 1

    def _parse_index(self):
        """
        TODO maybe the best thing to do here is to split up the grid
             according to some sort of optimization. For now, we stick
             to only one patch.
        """
        ngrids = self.num_grids
        self.grid_left_edge = np.zeros((ngrids, 3), dtype="float64")
        self.grid_right_edge = np.zeros((ngrids, 3), dtype="float64")
        self.grid_dimensions = np.zeros((ngrids, 3), dtype="int64")

        N1 = self._handle["header"]["n1"][()]
        N2 = self._handle["header"]["n2"][()]
        N3 = self._handle["header"]["n3"][()]
        for i in range(ngrids):
            X1, X2, X3 = self.dataset._get_boundaries()
            X1m, X2m, X3m = self.dataset._map_coordinates(X1, X2, X3)
            self.grid_left_edge[i] = np.array([X1m[0], X2m[0], X3m[0]], dtype="float64")
            self.grid_right_edge[i] = np.array(
                [X1m[0], X2m[1], X3m[1]], dtype="float64"
            )
            self.grid_dimensions[i] = np.array([N1, N2, N3], dtype="int64")

        self.grid_left_edge = self.ds.arr(self.grid_left_edge, "code_length")
        self.grid_right_edge = self.ds.arr(self.grid_right_edge, "code_length")

        self.grids = np.empty(ngrids, dtype="object")
        for i in range(ngrids):
            self.grids[i] = self.grid(i, self, 0)

        self.grid_particle_count = np.zeros([ngrids, 1], dtype="int64")
        self.max_level = 0

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()


class IHarmDataset(Dataset):

    _index_class = IHarmMKSHierarchy
    _field_info_class = IHarmFieldInfo
    _dataset_type = "iharm"
    _metric_type = None
    _metric_params = None
    _grid_filename = None

    _x1_info = [0.0, 1.0, 128]
    _x2_info = [0.0, 1.0, 128]
    _x3_info = [0.0, 1.0, 128]

    def __init__(
        self,
        filename,
        dataset_type=_dataset_type,
        storage_filename=None,
        units_override=None,
    ):

        self._handle = HDF5FileHandler(filename)
        self.fluid_types += ("iharm",)
        Dataset.__init__(self, filename, dataset_type, units_override=units_override)
        self.storage_filename = storage_filename
        self.filename = filename

    def _set_code_unit_attributes(self):
        """
        Load known units from parameter file (if required as for GRRMHD)
        or otherwise set based on defaults
        """

        # TODO actually load for GRRMHD
        self.length_unit = self.quan(1.0, "code_length")
        self.mass_unit = self.quan(1.0, "code_mass")
        self.time_unit = self.quan(1.0, "code_time")
        # TODO deal with magnetic fields

        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        # self.length_unit = self.quan(1.0, "cm")
        # self.mass_unit = self.quan(1.0, "g")
        # self.time_unit = self.quan(1.0, "s")
        # self.time_unit = self.quan(1.0, "s")
        #
        # These can also be set:
        # self.velocity_unit = self.quan(1.0, "cm/s")
        # self.magnetic_unit = self.quan(1.0, "gauss")

    def _get_boundaries(self):
        """
        Returns array of three 2D tuples corresponding to extent of
        boundary for this file.
        """
        dds = np.array(
            [self._geom_params[f"dx{i + 1}"] for i in range(3)], dtype="float64"
        )
        LE = np.array(
            [self._geom_params[f"startx{i + 1}"] for i in range(3)], dtype="float64"
        )
        RE = LE + self.domain_dimensions * dds
        return np.vstack((LE, RE)).T

    def _map_coordinates(self, X1, X2, X3, metric=None):
        """
        Convert X1,X2,X3 bounds to cartesian/spherical bounds
        according to metric. Assumes self._metric_params has
        been set.
        """
        if metric is None:
            metric = self._metric_type

        if metric == "MMKS":
            R = np.exp(X1)
            H = np.pi * X2 + (1.0 - self._metric_params["hslope"]) / 2.0 * np.sin(
                2.0 * np.pi * X2
            )
            P = X3
            return (R, H, P)

    def _parse_parameter_file(self):

        header = self._handle["header"]

        # general init
        self.unique_identifier = self.parameter_filename.__hash__()
        self.current_time = self._handle["t"][()]
        self._periodicity = (False, False, False)
        self.num_ghost_zones = 0
        # self._grid_filename = header["gridfile_name"][()]

        # get geometry
        self._metric_type = header["metric"][()].decode("ascii", "ignore")
        if self._metric_type.endswith("MKS"):
            self.geometry = "spherical"
            self._metric_params = {
                os.path.basename(param.name): param[()]
                for param in header["geom"][self._metric_type.lower()].values()
            }
        else:
            raise NotImplementedError

        self._geom_params = {}
        for key in ("dx1", "dx2", "dx3", "startx1", "startx2", "startx3"):
            self._geom_params[key] = header["geom"][key][()]
        self.domain_dimensions = np.array(
            [header[key][()] for key in ("n1", "n2", "n3")], dtype="int64"
        )

        # get and set grid dimensions
        X1, X2, X3 = self._get_boundaries()
        self._x1_info = [X1[0], X1[1], self.domain_dimensions[0]]
        self._x2_info = [X2[0], X2[1], self.domain_dimensions[1]]
        self._x3_info = [X3[0], X3[1], self.domain_dimensions[2]]
        (X1m, X2m, X3m) = self._map_coordinates(X1, X2, X3)
        self.domain_left_edge = np.array([X1m[0], X2m[0], X3m[0]], dtype="float64")
        self.domain_right_edge = np.array([X1m[1], X2m[1], X3m[1]], dtype="float64")
        self.domain_width = self.domain_right_edge - self.domain_left_edge

        # set fields directly in the field
        i = 0
        self._field_map = {}
        for name in header["prim_names"]:
            self._field_map[name.decode("ascii", "ignore")] = ("prims", i)

        # TODO add jcon, other fields of interest?

        self.refine_by = 1  # should this be zero?

        self.dimensionality = 3  # TODO figure out how many of the N's are 1

        self.field_ordering = "fortran"

        # not a cosmological simulation, so default the following to zero
        self.cosmological_simulation = 0
        self.current_redshift = 0.0
        self.omega_lambda = 0.0
        self.omega_matter = 0.0
        self.hubble_constant = 0.0

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            if args[0].endswith(".ih5"):
                return True
        except:
            pass
        return False
