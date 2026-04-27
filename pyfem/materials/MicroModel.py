# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Tuple

import h5py
import numpy as np

from pyfem.io.InputReader import InputRead
from pyfem.io.HDF5Writer import HDF5Writer
from pyfem.io.MeshWriter import MeshWriter
from pyfem.materials.BaseMaterial import BaseMaterial
from pyfem.solvers.Solver import Solver
from pyfem.util.dataStructures import Properties


class MicroModel(BaseMaterial):

  def __init__(self, props):

    BaseMaterial.__init__(self, props)

    if not hasattr(self, "rveFile"):
      raise ValueError("MicroModel requires an 'rveFile' property.")

    self.rank = getattr(props, "rank", 2)
    self.perturbation = getattr(self, "perturbation", 1.0e-7)
    self.sampleID = getattr(self, "sampleID", -1)
    self.exportFormat = getattr(self, "exportFormat", "vtu")
    self.exportPrefix = getattr(self, "exportPrefix", None)
    self.exportSamples = self._normalizeExportSamples(getattr(self, "exportSamples", []))

    if self.rank != 2:
      raise NotImplementedError("MicroModel currently supports 2D small-strain elements only.")

    self.outLabels = ["S11", "S22", "S12"]
    self.outData = np.zeros(3)

    self._committed_state = self._readMicroProblem(self.rveFile)
    self._trial_state = deepcopy(self._committed_state)

  def _readMicroProblem(self, rveFile: str) -> dict:

    fname = str(Path(rveFile).expanduser())
    props, globdat = InputRead(fname)

    if not hasattr(props, "models") or "rve" not in props.models:
      raise ValueError("MicroModel expects the micro problem to define models = ['rve'].")

    if not hasattr(props, "rve") or getattr(props.rve, "type", None) != "RVE":
      raise ValueError("MicroModel expects an RVE model in the micro problem.")

    props.solver = self._buildSolverProps(getattr(props, "solver", None))
    globdat.active = True

    return {"props": props, "globdat": globdat}

  def _buildSolverProps(self, solverProps: Properties | None) -> Properties:

    props = Properties()
    props.type = "NonlinearSolver"
    props.iterMax = getattr(solverProps, "iterMax", 20)
    props.tol = getattr(solverProps, "tol", 1.0e-6)
    props.maxCycle = 1
    props.loadTable = [1.0]

    return props

  def _normalizeExportSamples(self, exportSamples) -> list[int]:

    if exportSamples in [None, False]:
      return []

    if exportSamples == "all":
      return [-1]

    if isinstance(exportSamples, (int, np.integer)):
      return [int(exportSamples)]

    return [int(iSam) for iSam in exportSamples]

  def _normalizeFormats(self, fileFormat = None) -> list[str]:

    formats = fileFormat if fileFormat is not None else self.exportFormat

    if formats in [None, False]:
      return []

    if isinstance(formats, str):
      formats = [formats]

    normalized = []

    for fmt in formats:
      key = str(fmt).lower()

      if key in ["vtk", "vtu"]:
        normalized.append("vtu")
      elif key in ["h5", "hdf5"]:
        normalized.append("h5")
      else:
        raise ValueError(f"Unsupported micro export format '{fmt}'. Use 'vtu' or 'h5'.")

    return normalized

  def _shouldExport(self) -> bool:

    if len(self.exportSamples) == 0 or self.exportFormat in [None, False]:
      return False

    return -1 in self.exportSamples or self.sampleID in self.exportSamples

  def _getExportPrefix(self, prefix = None) -> str:

    if prefix is not None:
      exportPrefix = Path(prefix)
    elif self.exportPrefix is not None:
      exportPrefix = Path(self.exportPrefix)
    else:
      cycle = getattr(self.solverStat, "cycle", 0)
      rvePath = Path(self.rveFile)
      exportPrefix = rvePath.with_suffix("").parent / f"{rvePath.stem}_gp{self.sampleID}_macro_t{cycle}"

    exportPrefix.parent.mkdir(parents=True, exist_ok=True)

    return str(exportPrefix)

  def _makeWriterProps(self, writerName: str, writerType: str) -> Properties:

    writerProps = Properties({"type": writerType})
    props = Properties({writerName: writerProps})
    props.currentModule = writerName

    return props

  def storeMicroModel(self, fileFormat = None, prefix: str | None = None, state: dict | None = None) -> None:

    formats = self._normalizeFormats(fileFormat)

    if len(formats) == 0:
      return

    if state is None:
      state = self._committed_state

    globdat = state["globdat"]
    exportPrefix = self._getExportPrefix(prefix)

    oldPrefix = globdat.prefix
    oldCycle = globdat.solverStatus.cycle

    globdat.prefix = exportPrefix
    globdat.solverStatus.cycle = max(1, getattr(self.solverStat, "cycle", oldCycle))

    try:
      if "vtu" in formats:
        meshProps = self._makeWriterProps("microVtu", "MeshWriter")
        writer = MeshWriter(meshProps, globdat)
        writer.prefix = exportPrefix
        writer.writeCycle(globdat.state, meshProps, globdat)
        writer.writePvd()

      if "h5" in formats:
        h5Props = self._makeWriterProps("microH5", "HDF5Writer")
        writer = HDF5Writer(h5Props, globdat)
        writer.prefix = exportPrefix
        writer.singleFile = False

        with h5py.File(exportPrefix + ".h5", "w") as h5file:
          writer.writeCycle(h5file, globdat)
    finally:
      globdat.prefix = oldPrefix
      globdat.solverStatus.cycle = oldCycle

  def _solveMicroProblem(self, macroDStrain: np.ndarray) -> Tuple[np.ndarray, dict]:

    state = deepcopy(self._committed_state)
    props = state["props"]
    globdat = state["globdat"]

    unitStrain = np.asarray(macroDStrain, dtype=float)

    props.rve.unitStrain = unitStrain.tolist()

    for model in getattr(globdat.models, "modelss", []):
      if hasattr(model, "unitStrain"):
        model.unitStrain = unitStrain.copy()

    globdat.active = True
    globdat.dofs.createConstrainer()
    globdat.solverStatus.cycle = 0
    globdat.solverStatus.time = 0.0
    globdat.solverStatus.time0 = 0.0
    globdat.solverStatus.iiter = 0
    globdat.solverStatus.lam = 0.0

    solver = Solver(props, globdat)
    solver.run(props, globdat)

    return np.array(globdat.equivStress, copy=True), state

  def _computeTangent(self, baseStress: np.ndarray, macroDStrain: np.ndarray) -> np.ndarray:

    tang = np.zeros((3, 3))

    for i in range(3):
      perturbed = np.array(macroDStrain, copy=True)
      perturbed[i] += self.perturbation

      sigma, _ = self._solveMicroProblem(perturbed)
      tang[:, i] = (sigma - baseStress) / self.perturbation

    return tang

  def getStress(self, deformation):

    macroDStrain = np.asarray(deformation.dstrain, dtype=float)

    if len(macroDStrain) != 3:
      raise ValueError("MicroModel expects a 2D small-strain kinematics vector of length 3.")

    stress, state = self._solveMicroProblem(macroDStrain)

    if hasattr(state["globdat"], "equivTangent"):
      tang = np.array(state["globdat"].equivTangent, copy=True)
    else:
      tang = self._computeTangent(stress, macroDStrain)

    self._trial_state = state
    self.outData[:] = stress

    return stress, tang

  def commitHistory(self) -> None:

    BaseMaterial.commitHistory(self)
    self._committed_state = deepcopy(self._trial_state)

    if self._shouldExport():
      self.storeMicroModel(state=self._committed_state)
