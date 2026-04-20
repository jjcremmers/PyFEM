# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from .Element import Element

import numpy as np
from scipy.linalg import norm
from math import atan2, sin, cos, tan
from pyfem.elements.BeamUtil import glob2loc,loc2glob,BeamRotation


#===============================================================================
#
#===============================================================================

class ShapeFunc():

    def __init__( self ):
        self.n       = np.zeros(2)
        self.dn      = np.zeros(2)
        self.reduced = False
 
class Kinematic:
    pass    

class Beam3D ( Element ):

    #dofs per element
    dofTypes = [ 'u' , 'v' , 'w' , 'rx' , 'ry' , 'rz' ]

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    def __init__ ( self, elnodes , props ):
    
        Element.__init__( self, elnodes , props )

        self.length = None
        self.trans = np.zeros(shape=(3,3))
        self.referenceCoords = None
        self.referenceInitialized = False
        self.orientation = getattr( self, "orientation", None )
        self.sectionStiffness = self.buildSectionStiffness()
        self.sectionStiffnessReduced = self.sectionStiffness[:3, :3].copy()
        self.sectionStiffnessCurvature = self.sectionStiffness[3:, 3:].copy()
        (
            self.massPerLength,
            self.rotaryInertia1,
            self.rotaryInertia2,
            self.polarRotaryInertia,
        ) = self.buildMassProperties()
    
        gsr = ( -5.77350269189626e-1 , 5.77350269189626e-1 , 0.0 )
        gsa = ( 1.0 , 1.0 , 2.0 )
    
        self.shapes = [ShapeFunc() for i in range(3)]
    
        for iInt,shp in enumerate(self.shapes):
    
            shp.s = gsr[iInt]
        
            shp.n[0]   = 0.5*( 1.0 - shp.s)
            shp.n[1]   = 0.5*( 1.0 + shp.s)
        
            if iInt == 2:
                shp.reduced = True
                          
        self.beamRot = BeamRotation(3)
                                  
        self.family = "BEAM"

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    def __type__ ( self ):
        return self.__class__.__name__

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    def initializeReferenceGeometry( self, coords ):

        t = coords[1, :] - coords[0, :]

        self.length = norm( t )

        if self.length < 1.0e-14:
            raise RuntimeError( "Beam3D requires a non-zero element length." )

        t = t / self.length

        nx = self.getLocalXAxis( t )
        ny = np.cross( t, nx )
        ny /= norm( ny )

        self.trans[0, :] = nx
        self.trans[1, :] = ny
        self.trans[2, :] = t

        for iInt,shp in enumerate(self.shapes):
            shp.dn[0]  = -1.0 / self.length
            shp.dn[1]  =  1.0 / self.length
            shp.weight =  0.5 * self.length * ( 2.0 if iInt == 2 else 1.0 )

        self.referenceCoords = coords.copy()
        self.referenceInitialized = True

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    def ensureReferenceGeometry( self, coords ):

        if not self.referenceInitialized:
            self.initializeReferenceGeometry( coords )
            return

        if norm( coords - self.referenceCoords ) > 1.0e-12:
            raise RuntimeError(
                "Beam3D reference geometry changed after initialization."
            )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    def getLocalXAxis( self, t ):

        ref_vec = self.getOrientationVector( t )

        nx = ref_vec - ( ref_vec @ t ) * t
        nx_norm = norm( nx )

        if nx_norm < 1.0e-14:
            raise RuntimeError( "Beam3D failed to construct a local x-axis." )

        return nx / nx_norm

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    def getOrientationVector( self, t ):

        if self.orientation is not None:
            if len( self.orientation ) != 3:
                raise RuntimeError( "Beam3D orientation must contain three components." )

            ref_vec = np.zeros( 3 )
            ref_vec[:] = self.orientation[:]
            ref_norm = norm( ref_vec )

            if ref_norm < 1.0e-14:
                raise RuntimeError( "Beam3D orientation vector must be non-zero." )

            ref_vec /= ref_norm

            if abs( ref_vec @ t ) > 1.0 - 1.0e-8:
                raise RuntimeError(
                    "Beam3D orientation vector must not be parallel to the beam axis."
                )

            return ref_vec

        seed = np.zeros( 3 )
        seed[2] = 1.0

        if abs( t @ seed ) > 0.99:
            seed[:] = 0.0
            seed[0] = 1.0

        return seed

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    def buildSectionStiffness( self ):

        if hasattr( self, "sectionStiffness" ):
            data = self.sectionStiffness

            section = np.zeros( shape=( 6, 6 ) )

            if len( data ) == 36:
                for i in range( 6 ):
                    for j in range( 6 ):
                        section[i, j] = data[6 * i + j]
            elif len( data ) == 6:
                for i in range( 6 ):
                    if len( data[i] ) != 6:
                        raise RuntimeError(
                            "Beam3D sectionStiffness must be a 6 x 6 matrix."
                        )

                    for j in range( 6 ):
                        section[i, j] = data[i][j]
            else:
                raise RuntimeError(
                    "Beam3D sectionStiffness must contain 36 values or define "
                    "a 6 x 6 matrix."
                )

            coupling = np.zeros( shape=( 3, 3 ) )
            coupling[:] = section[:3, 3:]

            if norm( coupling ) > 1.0e-12 or norm( section[3:, :3] ) > 1.0e-12:
                raise RuntimeError(
                    "Beam3D sectionStiffness coupling between strain and curvature "
                    "blocks is not supported by the current reduced/full integration."
                )

            return section

        self.EA = getattr( self, "EA", self.E * self.A )
        self.GAy = getattr( self, "GAy", None )
        self.GAz = getattr( self, "GAz", None )
        self.EIy = getattr( self, "EIy", None )
        self.EIz = getattr( self, "EIz", None )
        self.GJ = getattr( self, "GJ", self.G * self.J )

        if self.GAy is None:
            if hasattr( self, "ky" ):
                self.GAy = self.ky * self.G * self.A
            else:
                self.GAy = 5.0 / 6.0 * self.G * self.A

        if self.GAz is None:
            if hasattr( self, "kz" ):
                self.GAz = self.kz * self.G * self.A
            else:
                self.GAz = 5.0 / 6.0 * self.G * self.A

        if self.EIy is None:
            self.EIy = self.E * self.Iy

        if self.EIz is None:
            self.EIz = self.E * self.Ix

        section = np.zeros( shape=( 6, 6 ) )

        section[0, 0] = self.GAy
        section[1, 1] = self.GAz
        section[2, 2] = self.EA
        section[3, 3] = self.EIy
        section[4, 4] = self.EIz
        section[5, 5] = self.GJ

        return section

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    def buildMassProperties( self ):

        mass_per_length = getattr( self, "massPerLength", None )
        rotary_inertia_1 = getattr( self, "rotaryInertia1", None )
        rotary_inertia_2 = getattr( self, "rotaryInertia2", None )
        polar_rotary_inertia = getattr( self, "polarRotaryInertia", None )

        if mass_per_length is None:
            mass_per_length = getattr( self, "rhoA", self.rho * self.A )

        if rotary_inertia_1 is None:
            rotary_inertia_1 = getattr( self, "rhoIx", self.rho * self.Ix )

        if rotary_inertia_2 is None:
            rotary_inertia_2 = getattr( self, "rhoIy", self.rho * self.Iy )

        if polar_rotary_inertia is None:
            polar_rotary_inertia = getattr( self, "rhoJ", self.rho * self.J )

        if mass_per_length <= 0.0:
            raise RuntimeError( "Beam3D massPerLength must be positive." )

        if rotary_inertia_1 < 0.0 or rotary_inertia_2 < 0.0 or polar_rotary_inertia < 0.0:
            raise RuntimeError( "Beam3D rotary inertia terms must be non-negative." )

        return mass_per_length, rotary_inertia_1, rotary_inertia_2, polar_rotary_inertia

#------------------------------------------------------------------------------
#  getTangentStiffness
#------------------------------------------------------------------------------

    def getTangentStiffness ( self, elemdat ):
        
        self.ensureReferenceGeometry( elemdat.coords )

        kin = self.glob2loc( elemdat )
        
        fint  = np.zeros(12)
        stiff = np.zeros(shape=(12,12))
    
        for iInt,shp in enumerate(self.shapes):
    
            self.updateKinematics( iInt , shp , kin )

            xi1 , xi2 = self.getXi( shp , kin.dphi )
        
            nmspat = self.getNMspat ( kin , shp.reduced )
        
            fint[:6] += (xi1 @ nmspat) * shp.weight
            fint[6:] += (xi2 @ nmspat) * shp.weight        
        
            cspat = self.getCspat ( kin.lambda_ , shp.reduced )
        
            stiff[:6,:6] += (xi1 @ (cspat @ xi1.T)) * shp.weight 
            stiff[:6,6:] += (xi1 @ (cspat @ xi2.T)) * shp.weight 
            stiff[6:,:6] += (xi2 @ (cspat @ xi1.T)) * shp.weight 
            stiff[6:,6:] += (xi2 @ (cspat @ xi2.T)) * shp.weight
        
            psi1 , psi2 = self.getPsi( shp )
        
            b = self.getGeomStiff( nmspat , kin.dphi )
        
            stiff[:6,:6] += (psi1 @ (b @ psi1.T)) * shp.weight 
            stiff[:6,6:] += (psi1 @ (b @ psi2.T)) * shp.weight 
            stiff[6:,:6] += (psi2 @ (b @ psi1.T)) * shp.weight 
            stiff[6:,6:] += (psi2 @ (b @ psi2.T)) * shp.weight    
        
        loc2glob( elemdat.fint  , self.trans , fint  )
        loc2glob( elemdat.stiff , self.trans , stiff )   
            
#------------------------------------------------------------------------------
#  getInternalForce
#------------------------------------------------------------------------------

    def getInternalForce ( self, elemdat ):
        self.ensureReferenceGeometry( elemdat.coords )

        kin = self.glob2loc( elemdat )
            
        fint  = np.zeros(12)
    
        for iInt,shp in enumerate(self.shapes):
    
            self.updateKinematics( iInt , shp , kin )

            xi1 , xi2 = self.getXi( shp , kin.dphi )
        
            nmspat = self.getNMspat ( kin , shp.reduced )
        
            fint[:6] += (xi1 @ nmspat) * shp.weight
            fint[6:] += (xi2 @ nmspat) * shp.weight        
                
        loc2glob( elemdat.fint  , self.trans , fint  )
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def commit ( self, elemdat ):

        self.beamRot.commit()        
         
#------------------------------------------------------------------------------
#  getCspat
#------------------------------------------------------------------------------
    
    def getCspat( self , lambda_ , reduced ):

        pi = np.zeros(shape=(6,6))

        pi[:3,:3] = lambda_
        pi[3:,3:] = lambda_
    
        cmat = np.zeros(shape=(6,6))
  
        if reduced:
            cmat[:3, :3] = self.sectionStiffnessReduced
        else:
            cmat[3:, 3:] = self.sectionStiffnessCurvature
  
        return pi @ (cmat @ pi.T)  

#-------------------------------------------------------------------------------
#  getNMspat
#-------------------------------------------------------------------------------


    def getNMspat( self , kin , reduced ):

        pi = np.zeros(shape=(6,6))

        pi[:3,:3] = kin.lambda_
        pi[3:,3:] = kin.lambda_
       
        nmmat = np.zeros(6)
     
        if reduced:
            nmmat[:3] = self.sectionStiffnessReduced @ kin.gamma
        else:
            nmmat[3:] = self.sectionStiffnessCurvature @ kin.kappa
 
  
        return pi @ nmmat

#------------------------------------------------------------------------------
#  getGeomStiff
#------------------------------------------------------------------------------

    def getGeomStiff( self, nmspat, dphi ):

        b = np.zeros(shape=(9,9))
    
        ndotp = nmspat[:3] @ dphi

        b[7,0] =  nmspat[2]
        b[8,0] = -nmspat[1]
        b[6,1] = -nmspat[2]
        b[8,1] =  nmspat[0]
        b[6,2] =  nmspat[1]
        b[7,2] = -nmspat[0]

        b[1,6] = -nmspat[2]
        b[2,6] =  nmspat[1]
        b[0,7] =  nmspat[2]
        b[2,7] = -nmspat[0]
        b[0,8] = -nmspat[1]
        b[1,8] =  nmspat[0]

        b[4,6] = -nmspat[5]
        b[5,6] =  nmspat[4]
        b[3,7] =  nmspat[5]
        b[5,7] = -nmspat[3]
        b[3,8] = -nmspat[4]
        b[4,8] =  nmspat[3]

        b[6,6] = nmspat[0] * dphi[0] - ndotp
        b[7,6] = nmspat[1] * dphi[0]
        b[8,6] = nmspat[2] * dphi[0]
        b[6,7] = nmspat[0] * dphi[1]
        b[7,7] = nmspat[1] * dphi[1] - ndotp
        b[8,7] = nmspat[2] * dphi[1]
        b[6,8] = nmspat[0] * dphi[2]
        b[7,8] = nmspat[1] * dphi[2]
        b[8,8] = nmspat[2] * dphi[2] - ndotp

        return b
 
#-------------------------------------------------------------------------------
#  getXi
#-------------------------------------------------------------------------------

    def getXi( self , shp , dphi ):
    
        xi1 = np.zeros(shape=(6,6))
        xi2 = np.zeros(shape=(6,6))
        
        for i in range(6):
            xi1[i, i] = shp.dn[0]

        xi1[4, 0] = -shp.n[0] * dphi[2]
        xi1[5, 0] =  shp.n[0] * dphi[1]
        xi1[3, 1] =  shp.n[0] * dphi[2]
        xi1[5, 1] = -shp.n[0] * dphi[0]
        xi1[3, 2] = -shp.n[0] * dphi[1]
        xi1[4, 2] =  shp.n[0] * dphi[0]
   
        for i in range(6):
            xi2[i, i] = shp.dn[1]

        xi2[4, 0] = -shp.n[1] * dphi[2]
        xi2[5, 0] =  shp.n[1] * dphi[1]
        xi2[3, 1] =  shp.n[1] * dphi[2]
        xi2[5, 1] = -shp.n[1] * dphi[0]
        xi2[3, 2] = -shp.n[1] * dphi[1]
        xi2[4, 2] =  shp.n[1] * dphi[0]
        
        return xi1,xi2
        
#-------------------------------------------------------------------------------
#  getPsi
#-------------------------------------------------------------------------------

    def getPsi( self , shp ):
    
        psi1 = np.zeros(shape=(6,9))
        psi2 = np.zeros(shape=(6,9))
        
   
        for i in range(6):
            psi1[i, i] = shp.dn[0]

        psi1[3, 6] = shp.n[0]
        psi1[4, 7] = shp.n[0]
        psi1[5, 8] = shp.n[0]
       
        for i in range(6):
            psi2[i, i] = shp.dn[1]

        psi2[3, 6] = shp.n[1]
        psi2[4, 7] = shp.n[1]
        psi2[5, 8] = shp.n[1]
        
        return psi1,psi2

#-------------------------------------------------------------------------------
#  updateKinematics
#-------------------------------------------------------------------------------

    def updateKinematics( self, iInt , shp , kin ):
    
         kin.theta  = shp.n[0]  * kin.ddis[3:6] + shp.n[1]  * kin.ddis[9:]
         kin.dtheta = shp.dn[0] * kin.ddis[3:6] + shp.dn[1] * kin.ddis[9:]
    
         kin.dphi   = shp.dn[0] * kin.disp[:3] + shp.dn[1] * kin.disp[6:9]
         kin.dphi[2] += 1.0

         self.beamRot.update(iInt, kin )

         kin.gamma    = kin.lambda_.T @ kin.dphi
         kin.gamma[2] -= 1.0                  
         
#-------------------------------------------------------------------------------
#  glob2loc
#-------------------------------------------------------------------------------
        
    def glob2loc( self , elemdat ):
    
        kin = Kinematic()
            
        kin.disp = np.zeros(12)
        kin.ddis = np.zeros(12)
        
        glob2loc( kin.disp , self.trans , elemdat.state  )
        glob2loc( kin.ddis , self.trans , elemdat.Dstate )
        
        return kin
        
#-------------------------------------------------------------------------------
#  getMassMatrix
#-------------------------------------------------------------------------------

    def getMassMatrix( self , elemdat ):
        self.ensureReferenceGeometry( elemdat.coords )
        
        l  = self.length
        l2 = self.length*self.length
        mu = self.massPerLength

        mass = np.zeros( shape=(12,12) )

        mass[0,0]   = 1./3.
        mass[1,1]   = 13./35. + 6.*self.rotaryInertia1/(5.*mu*l2)
        mass[2,2]   = 13./35. + 6.*self.rotaryInertia2/(5.*mu*l2)
        mass[3,3]   = self.polarRotaryInertia/(3.0*mu)
        mass[4,4]   = l2/105.0 + 2.0*self.rotaryInertia2/(15.0*mu)
        mass[5,5]   = l2/105.0 + 2.0*self.rotaryInertia1/(15.0*mu)
        mass[6,6]   = mass[0,0]
        mass[7,7]   = mass[1,1]
        mass[8,8]   = mass[2,2]
        mass[9,9]   = mass[3,3]
        mass[10,10] = mass[4,4]
        mass[11,11] = mass[5,5]

        mass[6,0]   = 1./6.
        mass[0,6]   = mass[6,0]
        
        mass[5,1]   = 11.0*l/210.0 + self.rotaryInertia1/(10.0*mu*l)
        mass[1,5]   = mass[5,1]       
        mass[7,1]   = 9.0/70.0 - 6.0*self.rotaryInertia1/(5.0*mu*l2)
        mass[1,7]   = mass[7,1]
        mass[11,1]  = -13.0*l/420.0 + self.rotaryInertia1/(10.0*mu*l)
        mass[1,11]  = mass[11,1]        
        
        mass[4,2]   = -11.0*l/210.0 - self.rotaryInertia2/(10.0*mu*l)
        mass[2,4]   = mass[4,2]     
        mass[8,2]   = 9.0/70.0 - 6.0*self.rotaryInertia2/(5.0*mu*l2)
        mass[2,8]   = mass[8,2]          
        mass[10,2]  = 13.0*l/420.0 - self.rotaryInertia2/(10.0*mu*l)
        mass[2,10]  = mass[10,2]
        
        mass[9,3]   = self.polarRotaryInertia/(6.0*mu)
        mass[3,9]   = mass[9,3]

        mass[8,4]   = -13.0*l/420.0 + self.rotaryInertia2/(10.0*mu*l)
        mass[4,8]   = mass[8,4]
        mass[10,4]  = -l2/140.0 - self.rotaryInertia2/(30.0*mu)
        mass[4,10]  = mass[10,4]

        mass[7,5]   = 13.0*l/420.0 - self.rotaryInertia1/(10.0*mu*l)
        mass[5,7]   = mass[7,5]
        mass[11,5]  = -l2/140.0 - self.rotaryInertia1/(30.0*mu)
        mass[5,11]  = mass[11,5]

        mass[11,7]  = -11.0*l/210.0 - self.rotaryInertia1/(10.0*mu*l)
        mass[7,11]  = mass[11,7]       

        mass[10,8]  = 11.0*l/210.0 + self.rotaryInertia2/(10.0*mu*l)
        mass[8,10]  = mass[10,8]

        mass *= mu*l

        loc2glob( elemdat.mass , self.trans , mass )

        elemdat.lumped = sum(elemdat.mass)
        
