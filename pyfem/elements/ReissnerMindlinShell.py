# SPDX-License-Identifier: MIT
# Copyright (c) 2011-2026 Joris J.C. Remmers

"""Reissner-Mindlin shell element implementation."""

from numpy import array, cross, eye, sum, zeros
from numpy.linalg import norm

from pyfem.elements.Composite import Laminate
from pyfem.util.matrixUtils import skew
from pyfem.util.shapeFunctions import getElemShapeData, getShapeQuad4
from pyfem.util.utilFunctions import one_minus_cos_over_x2, sin_over_x

from .Element import Element


class PostProcessPoint:
    """Container describing a through-thickness output location."""


class CurrentBasisData:
    """Container for the current shell basis vectors."""


class KinematicOperatorData:
    """Container for the kinematic operators used in linearization."""


class ReferenceBasisData:
    """Container for the reference shell basis vectors and mappings."""


class ReissnerMindlinShell(Element):
    """Four-node Reissner-Mindlin shell element with layered laminate support."""

    dofTypes = ["u", "v", "w", "rx", "ry", "rz"]

    def __init__(self, elnodes, props):
        """Initialize the shell element and its laminate properties."""
        super().__init__(elnodes, props)

        self.material = Laminate(props)
        self.family = "SHELL"

        self.tangentEps = 1.0e-7
        if hasattr(props, "tangentPerturbation"):
            self.tangentEps = props.tangentPerturbation

        self.zetaSample = array([-0.5773502691896257, 0.5773502691896257])
        self.zetaWeight = array([1.0, 1.0])

        self.drillingScale = 1.0e-6
        if hasattr(props, "drillingScale"):
            self.drillingScale = props.drillingScale

        self.inertia = self.material.getMassInertia()

        self.initPostProcessing()

    def getTangentStiffness(self, elemdat):
        """Assemble the element internal force vector and tangent stiffness."""
        elemdat.fint[:] = 0.0
        elemdat.stiff[:, :] = 0.0

        s_data = self.getShapeData(elemdat.coords)
        node_directors = self.getReferenceDirectors(elemdat.coords)

        n_nel = elemdat.coords.shape[0]
        n_dof = len(self.dofTypes) * n_nel

        for shape_data in s_data:
            ref = self.getReferenceBasis(elemdat.coords, node_directors, shape_data)

            for i_lay, (layer, zeta_data) in enumerate(self.iterateLayers()):
                c_mat = self.getLayerMatrix(layer)

                for zeta, weight in zeta_data:
                    cur = self.getCurrentBasis(
                        elemdat.coords,
                        node_directors,
                        ref,
                        shape_data,
                        elemdat.state,
                        zeta,
                    )
                    op = self.getKinematicOperators(
                        shape_data,
                        ref,
                        cur,
                        elemdat.state,
                        node_directors,
                        zeta,
                    )
                    strain = self.getStrainFromBasis(cur)
                    stress = c_mat @ strain

                    elemdat.fint += (
                        op.B.transpose() @ stress
                        * shape_data.weight
                        * weight
                    )
                    elemdat.stiff += (
                        op.B.transpose() @ (c_mat @ op.B)
                        + self.getGeometricStiffness(op, stress, n_dof)
                    ) * shape_data.weight * weight

                    if hasattr(self, "globdat"):
                        self.storeLayerOutput(stress, i_lay, zeta, shape_data.weight * weight)

        self.addDrillingContribution(
            elemdat.stiff,
            elemdat.fint,
            elemdat.coords,
            elemdat.state,
        )
        elemdat.stiff = 0.5 * (elemdat.stiff + elemdat.stiff.transpose())

    def getInternalForce(self, elemdat):
        """Compute the element internal force vector."""
        elemdat.fint = self.getLocalInternalForce(
            elemdat.coords,
            elemdat.state,
            True,
        )

    def getMassMatrix(self, elemdat):
        """Assemble the consistent and lumped element mass matrices."""
        s_data = self.getShapeData(elemdat.coords)

        n_nel = elemdat.coords.shape[0]
        n_dof = len(self.dofTypes) * n_nel

        mass = zeros(shape=(n_dof, n_dof))

        for shape_data in s_data:
            for i_nod in range(n_nel):
                for j_nod in range(n_nel):
                    shp = shape_data.h[i_nod] * shape_data.h[j_nod] * shape_data.weight

                    for k in range(3):
                        mass[6 * i_nod + k, 6 * j_nod + k] += self.inertia[0] * shp

                    for k in range(3, 6):
                        mass[6 * i_nod + k, 6 * j_nod + k] += self.inertia[2] * shp

        elemdat.mass = mass
        elemdat.lumped = sum(mass)

    def getLocalInternalForce(self, coords, state, storeOutput):
        """Return the local internal force vector for the given state."""
        s_data = self.getShapeData(coords)
        node_directors = self.getReferenceDirectors(coords)

        n_nel = coords.shape[0]
        fint = zeros(len(self.dofTypes) * n_nel)

        for shape_data in s_data:
            ref = self.getReferenceBasis(coords, node_directors, shape_data)

            for i_lay, (layer, zeta_data) in enumerate(self.iterateLayers()):
                c_mat = self.getLayerMatrix(layer)

                for zeta, weight in zeta_data:
                    cur = self.getCurrentBasis(
                        coords,
                        node_directors,
                        ref,
                        shape_data,
                        state,
                        zeta,
                    )
                    op = self.getKinematicOperators(
                        shape_data,
                        ref,
                        cur,
                        state,
                        node_directors,
                        zeta,
                    )
                    strain = self.getStrainFromBasis(cur)
                    stress = c_mat @ strain

                    fint += op.B.transpose() @ stress * shape_data.weight * weight

                    if storeOutput and hasattr(self, "globdat"):
                        self.storeLayerOutput(
                            stress,
                            i_lay,
                            zeta,
                            shape_data.weight * weight,
                        )

        self.addDrillingContribution(None, fint, coords, state)

        return fint

    def getStrainVector(self, coords, nodeDirectors, ref, shapeData, state, zeta):
        """Compute the generalized strain vector at a thickness coordinate."""
        cur = self.getCurrentBasis(
            coords,
            nodeDirectors,
            ref,
            shapeData,
            state,
            zeta,
        )
        return self.getStrainFromBasis(cur)

    def getStrainFromBasis(self, cur):
        """Build the generalized shell strain vector from basis vectors."""
        strain = zeros(5)

        strain[0] = 0.5 * ((cur.a1 @ cur.a1) - 1.0)
        strain[1] = 0.5 * ((cur.a2 @ cur.a2) - 1.0)
        strain[2] = cur.a1 @ cur.a2
        strain[3] = cur.a1 @ cur.d
        strain[4] = cur.a2 @ cur.d

        return strain

    def getCurrentBasis(self, coords, nodeDirectors, ref, shapeData, state, zeta):
        """Construct the current covariant basis vectors and director."""
        cur = CurrentBasisData()

        g1 = zeros(3)
        g2 = zeros(3)
        dd1 = zeros(3)
        dd2 = zeros(3)
        dvec = zeros(3)

        for i_nod in range(len(shapeData.h)):
            base = 6 * i_nod

            x_vec = coords[i_nod, :] + state[base : base + 3]
            director = self.getDirector(
                state[base + 3 : base + 6],
                nodeDirectors[i_nod],
            )

            g1 += shapeData.dhdxi[i_nod, 0] * x_vec
            g2 += shapeData.dhdxi[i_nod, 1] * x_vec
            dd1 += shapeData.dhdxi[i_nod, 0] * director
            dd2 += shapeData.dhdxi[i_nod, 1] * director
            dvec += shapeData.h[i_nod] * director

        g1 += zeta * dd1
        g2 += zeta * dd2

        cur.a1 = ref.c1[0] * g1 + ref.c1[1] * g2
        cur.a2 = ref.c2[0] * g1 + ref.c2[1] * g2
        cur.d = dvec

        return cur

    def getKinematicOperators(
        self,
        shapeData,
        ref,
        cur,
        state,
        nodeDirectors,
        zeta,
    ):
        """Build the linearized kinematic operators for the shell element."""
        n_nel = len(shapeData.h)
        n_dof = 6 * n_nel

        op = KinematicOperatorData()
        op.B = zeros(shape=(5, n_dof))
        op.A1 = []
        op.A2 = []
        op.D = []

        for i_nod in range(n_nel):
            base = 6 * i_nod

            dmat = self.getDirectorJacobian(
                state[base + 3 : base + 6],
                nodeDirectors[i_nod],
            )

            alpha1 = (
                ref.c1[0] * shapeData.dhdxi[i_nod, 0]
                + ref.c1[1] * shapeData.dhdxi[i_nod, 1]
            )
            alpha2 = (
                ref.c2[0] * shapeData.dhdxi[i_nod, 0]
                + ref.c2[1] * shapeData.dhdxi[i_nod, 1]
            )

            a1_mat = zeros(shape=(3, 6))
            a2_mat = zeros(shape=(3, 6))
            d_mat = zeros(shape=(3, 6))

            a1_mat[:, 0:3] = alpha1 * eye(3)
            a2_mat[:, 0:3] = alpha2 * eye(3)

            a1_mat[:, 3:6] = zeta * alpha1 * dmat
            a2_mat[:, 3:6] = zeta * alpha2 * dmat
            d_mat[:, 3:6] = shapeData.h[i_nod] * dmat

            op.A1.append(a1_mat)
            op.A2.append(a2_mat)
            op.D.append(d_mat)

            op.B[0, base : base + 6] = cur.a1 @ a1_mat
            op.B[1, base : base + 6] = cur.a2 @ a2_mat
            op.B[2, base : base + 6] = (cur.a2 @ a1_mat) + (cur.a1 @ a2_mat)
            op.B[3, base : base + 6] = (cur.d @ a1_mat) + (cur.a1 @ d_mat)
            op.B[4, base : base + 6] = (cur.d @ a2_mat) + (cur.a2 @ d_mat)

        return op

    def getGeometricStiffness(self, op, stress, nDof):
        """Assemble the stress-dependent geometric stiffness contribution."""
        stiff = zeros(shape=(nDof, nDof))

        for i_nod in range(len(op.A1)):
            i_base = 6 * i_nod
            for j_nod in range(len(op.A1)):
                j_base = 6 * j_nod

                block = stress[0] * (op.A1[i_nod].transpose() @ op.A1[j_nod])
                block += stress[1] * (op.A2[i_nod].transpose() @ op.A2[j_nod])
                block += stress[2] * (
                    (op.A1[i_nod].transpose() @ op.A2[j_nod])
                    + (op.A2[i_nod].transpose() @ op.A1[j_nod])
                )
                block += stress[3] * (
                    (op.A1[i_nod].transpose() @ op.D[j_nod])
                    + (op.D[i_nod].transpose() @ op.A1[j_nod])
                )
                block += stress[4] * (
                    (op.A2[i_nod].transpose() @ op.D[j_nod])
                    + (op.D[i_nod].transpose() @ op.A2[j_nod])
                )

                stiff[i_base : i_base + 6, j_base : j_base + 6] += block

        return stiff

    def addDrillingContribution(self, stiff, fint, coords, state):
        """Add the artificial drilling stiffness and corresponding forces."""
        area = self.getArea(coords)
        c_mem = max(1.0, self.getMembraneStiffnessScale())
        k_drill = self.drillingScale * area * c_mem

        for i_nod in range(coords.shape[0]):
            fint[6 * i_nod + 5] += k_drill * state[6 * i_nod + 5]
            if stiff is not None:
                stiff[6 * i_nod + 5, 6 * i_nod + 5] += k_drill

    def getMembraneStiffnessScale(self):
        """Return a scalar measure of the membrane stiffness level."""
        amat = self.material.getA()
        return abs(amat[0, 0]) + abs(amat[1, 1]) + 2.0 * abs(amat[2, 2])

    def iterateLayers(self):
        """Yield layer data together with through-thickness integration points."""
        for i_lay, layer in enumerate(self.material.layers):
            z0 = self.material.h[i_lay]
            z1 = self.material.h[i_lay + 1]
            zdat = []

            for xi, weight in zip(self.zetaSample, self.zetaWeight):
                zeta = 0.5 * (z0 + z1) + 0.5 * (z1 - z0) * xi
                zdat.append((zeta, 0.5 * (z1 - z0) * weight))

            yield layer, zdat

    def getLayerMatrix(self, layer):
        """Return the constitutive matrix for a laminate layer."""
        name = layer.mat
        theta = layer.theta

        c_mat = zeros(shape=(5, 5))
        c_mat[:3, :3] = self.material.materials[name].getQbar(theta)
        c_mat[3:, 3:] = (
            self.material.shearCorr
            * self.material.materials[name].getQshearbar(theta)
        )

        return c_mat

    def getDirector(self, rot, d0):
        """Rotate a reference director using the nodal rotation vector."""
        theta = norm(rot)

        if theta < 1.0e-14:
            return d0.copy()

        axis = rot / theta
        ax = skew(axis)

        rotation = eye(3)
        rotation += sin_over_x(theta) * theta * ax
        rotation += one_minus_cos_over_x2(theta) * theta * theta * (ax @ ax)

        return rotation @ d0

    def getDirectorJacobian(self, rot, d0):
        """Compute the director derivative with respect to the rotation vector."""
        jac = zeros(shape=(3, 3))

        for i_dir in range(3):
            h_val = self.getPerturbation(rot[i_dir])
            drot = zeros(3)
            drot[i_dir] = h_val

            d_plus = self.getDirector(rot + drot, d0)
            d_min = self.getDirector(rot - drot, d0)

            jac[:, i_dir] = (d_plus - d_min) / (2.0 * h_val)

        return jac

    def getShapeData(self, coords):
        """Return quadrature and shape-function data for the element."""
        if coords.shape[0] != 4:
            raise NotImplementedError(
                "ReissnerMindlinShell currently only supports Quad4 elements."
            )

        return getElemShapeData(coords, elemType="Quad4")

    def getPerturbation(self, value):
        """Return a finite-difference perturbation scaled to the local magnitude."""
        scale = max(1.0, abs(value), self.material.thick)
        return self.tangentEps * scale

    def getReferenceDirectors(self, coords):
        """Construct reference directors at the shell nodes."""
        node_xi = [
            array([-1.0, -1.0]),
            array([1.0, -1.0]),
            array([1.0, 1.0]),
            array([-1.0, 1.0]),
        ]

        dirs = zeros(shape=(coords.shape[0], 3))

        for i_nod, xi_vec in enumerate(node_xi):
            sdat = getShapeQuad4(xi_vec)
            g1_vec = coords.transpose() @ sdat.dhdxi[:, 0]
            g2_vec = coords.transpose() @ sdat.dhdxi[:, 1]
            dirs[i_nod, :] = self.unit(cross(g1_vec, g2_vec))

        return dirs

    def getReferenceBasis(self, coords, nodeDirectors, shapeData):
        """Construct the reference basis and local in-plane mapping tensors."""
        ref = ReferenceBasisData()

        g1_vec = coords.transpose() @ shapeData.dhdxi[:, 0]
        g2_vec = coords.transpose() @ shapeData.dhdxi[:, 1]

        g_cov = zeros(shape=(2, 2))
        g_cov[0, 0] = g1_vec @ g1_vec
        g_cov[0, 1] = g1_vec @ g2_vec
        g_cov[1, 0] = g_cov[0, 1]
        g_cov[1, 1] = g2_vec @ g2_vec

        det_g = g_cov[0, 0] * g_cov[1, 1] - g_cov[0, 1] * g_cov[1, 0]
        if abs(det_g) < 1.0e-14:
            raise RuntimeError("Degenerated shell metric.")

        g_inv = zeros(shape=(2, 2))
        g_inv[0, 0] = g_cov[1, 1] / det_g
        g_inv[0, 1] = -g_cov[0, 1] / det_g
        g_inv[1, 0] = -g_cov[1, 0] / det_g
        g_inv[1, 1] = g_cov[0, 0] / det_g

        g_sup1 = g_inv[0, 0] * g1_vec + g_inv[0, 1] * g2_vec
        g_sup2 = g_inv[1, 0] * g1_vec + g_inv[1, 1] * g2_vec

        ref.e3 = self.unit(nodeDirectors.transpose() @ shapeData.h)

        tmp = g1_vec - (g1_vec @ ref.e3) * ref.e3
        if norm(tmp) < 1.0e-14:
            tmp = g2_vec - (g2_vec @ ref.e3) * ref.e3

        ref.e1 = self.unit(tmp)
        ref.e2 = self.unit(cross(ref.e3, ref.e1))
        ref.e1 = self.unit(cross(ref.e2, ref.e3))

        ref.c1 = array([ref.e1 @ g_sup1, ref.e1 @ g_sup2])
        ref.c2 = array([ref.e2 @ g_sup1, ref.e2 @ g_sup2])

        return ref

    def getArea(self, coords):
        """Compute the shell area from two triangular sub-surfaces."""
        return 0.5 * norm(cross(coords[1] - coords[0], coords[3] - coords[0])) + 0.5 * norm(
            cross(coords[2] - coords[1], coords[3] - coords[1])
        )

    def unit(self, a_vec):
        """Return a unit vector in the direction of the input vector."""
        a_len = norm(a_vec)
        if a_len < 1.0e-14:
            raise RuntimeError("Zero-length vector encountered in shell basis construction.")

        return a_vec / a_len

    def storeLayerOutput(self, stress, iLay, zeta, weight):
        """Store selected through-thickness stress outputs for post-processing."""
        for post_process_point in self.postProcess:
            if abs(zeta - post_process_point.z) < 0.51 * self.material.layers[iLay].thick:
                self.appendNodalOutput(post_process_point.labels, stress[:3], weight)

        self.appendNodalOutput(["q13", "q23"], stress[3:], weight)

    def initPostProcessing(self):
        """Initialize through-thickness post-processing output locations."""
        layer_count = self.material.layerCount()

        self.postProcess = []

        point = PostProcessPoint()
        point.z = -0.5 * self.material.thick
        point.labels = ["s11bot", "s22bot", "s12bot"]
        self.postProcess.append(point)

        if layer_count > 1:
            i_mid = int(0.5 * layer_count)
            point = PostProcessPoint()
            point.z = 0.5 * (self.material.h[i_mid] + self.material.h[i_mid + 1])
            point.labels = ["s11mid", "s22mid", "s12mid"]
            self.postProcess.append(point)

        point = PostProcessPoint()
        point.z = 0.5 * self.material.thick
        point.labels = ["s11top", "s22top", "s12top"]
        self.postProcess.append(point)
