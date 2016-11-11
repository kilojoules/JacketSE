#-------------------------------------------------------------------------------
# Name:        VarTrees.py
# Purpose:   Stores common VarTrees that are used by multiple modules (files).
#            This avoids circular calls.
#
# Author:      rdamiani
#
# Created:     17/02/2014
# Copyright:   (c) rdamiani 2014
# Licence:     Apache C 2014
#-------------------------------------------------------------------------------
from openmdao.main.api import VariableTree,  Component, Assembly, set_as_top
from openmdao.main.datatypes.api import Int, Float, Array, VarTree, Bool, Dict,Instance
from commonse.Tube import Tube
from commonse.RigidMember import RigidMember
import numpy as np
#______________________________________________________________________________#

def RNAprops(self):
        """Basic Inertial and Geometric Properties of RNA"""
        self.add_param("RNApropos:mass", 0., units='kg', pass_by_object=True) # #RNA mass [kg]
        self.add_param("RNApropos:I", np.zeros(6), units='kg*m**2', pass_by_object=True) # 'RNA [IXX,IYY,IZZ,IXY,IXZ,IYZ] @tower top flange')
        self.add_param("RNApropos:CMoff", np.zeros(3), units='m', pass_by_object=True) # 'RNA CM [x,y,z] offset from Tower Top Flange'
        self.add_param("RNApropos:Thoff", np.zeros(3), units='m', pass_by_object=True) # 'Rotor Hub Center [x,y,z] offset from Tower Top Flange
        self.add_param("RNApropos:rna_weightM", True, pass_by_object=True) # flag to consider or not the RNA weight effect on Moment
        self.add_param("RNApropos:yawangle", 0., units='m', pass_by_object=True) # YAW angle CCW RH rule, to account for possible nacelle weight contribution
       # mass = Float(   units='kg',    desc='RNA mass')    #RNA mass [kg]
       # I = Array(np.zeros(6), dtype=np.float, units='kg*m**2',desc='RNA [IXX,IYY,IZZ,IXY,IXZ,IYZ] @tower top flange')
       # CMoff=Array(np.zeros(3), dtype=np.float,units='m', desc='RNA CM [x,y,z] offset from Tower Top Flange')       # RNA CMx,y,zoff [m]
       # Thoff=Array(np.zeros(3), dtype=np.float,units='m',desc='Rotor Hub Center [x,y,z] offset from Tower Top Flange')       # Thrust point of application [m]
       # rna_weightM = Bool(True, units=None, desc='flag to consider or not the RNA weight effect on Moment')
        # yawangle=Float(0.,   units='deg', desc='YAW angle CCW RH rule, to account for possible nacelle weight contribution.')       # RNA Yaw angle [deg]

def TPmassprops(self):
        """Basic Inertial and Geometric Properties of TP additional mass for towerSE"""
        self.add_param('TPmassprops:mass', 0., units='kg', pass_by_object=True) # TP lumped mass
        self.add_param("TPmassprops:I", np.zeros(6), units='', pass_by_object=True) # TP lumped mass [IXX,IYY,IZZ,IXY,IXZ,IYZ] @deck height (=tower-base flange)
        self.add_param("TPmassprops:CMoff", np.zeros(3), units='kg*m**2', pass_by_object=True) # TP lumped mass CM [x/DTP,y/DTP,z/TPlength] offset from deck height
 #   mass=Float(   units='kg',    desc='TP lumped mass')    #TP mass [kg]
 #   I    =Array(np.zeros(6), dtype=np.float, units='kg*m**2',desc='TP lumped mass [IXX,IYY,IZZ,IXY,IXZ,IYZ] @deck height (=tower-base flange)')
 #   CMoff=Array(np.zeros(3), dtype=np.float,units=None, desc='TP lumped mass CM [x/DTP,y/DTP,z/TPlength] offset from deck height')
 #   #Thoff=Array(np.zeros(3), dtype=np.float,units='m',desc='Rotor Hub Center [x,y,z] offset from Tower Top Flange')       # Thrust point of application [m]
 #   #rna_weightM = Bool(True, units=None, desc='flag to consider or not the RNA weight effect on Moment')
 #   #yawangle=Float(0.,   units='deg', desc='YAW angle CCW RH rule, to account for possible nacelle weight contribution.')       # RNA Yaw angle [deg]


def JcktGeoOutputs(self)::
        """Node Coordinates and Member Connectivity."""
        self.add_param("JcktGeoOutputs:nNodesJckt", 0, pass_by_object=True) # Total Number of nodes in the Jacket, no Tower
        self.add_param("JcktGeoOutputs:nodes", np.array([]), units='m', pass_by_object=True) # Node''s coordinates in the Substructure
        self.add_param("JcktGeoOutputs:radii", np.array([]), units='m', pass_by_object=True) # Node''s Radii
        self.add_param("JcktGeoOutputs:Reacts", np.array([]), pass_by_object=True) # Node IDs for reactions + Fixity values (1/0)'
        self.add_param("JcktGeoOutputs:nmems", 0, pass_by_object=True) # Total Number of Members in the Jacket, no Tower
        self.add_param("JcktGeoOutputs:mems", np.array([]), units='', pass_by_object=True) # Member Connectivity Node i Node j for every member (i.e.,element)
        self.add_param("JcktGeoOutputs:props", np.array([]), units='', pass_by_object=True) # Jacket Member's xsectional and material properties
        self.add_param("JcktGeoOutputs:XnsfM", np.array([]), pass_by_object=True) # Jacket Member''s Global to Local DIrection Cosine Matrices [3,3,nelems]
        self.JcktGeoOutputs_TubeObjs = Tube() # Object of Class Tube for Jacket Elements. Tube Objects, one per element
        self.add_param("JcktGeoOutputs:jnt_masses", np.array([]), units='', pass_by_object=True) # Jacket Concentrated Masses at Joints: (joint no, mass, IXX, IYY, IZZ)
        self.add_param("JcktGeoOutputs:jnt_masses_yaw", np.array([]), units='', pass_by_object=True) # Jacket Concentrated Masses at Joints: (joint no, mass, IXX, IYY, IZZ) accounting for RNA yaw
    #nNodesJckt = Int(  units=None,                 desc='Total Number of nodes in the Jacket, no Tower')
    #nodes =      Array(units='m',  dtype=np.float, desc='Node''s coordinates in the Substructure')
    #radii =      Array(units='m',  dtype=np.float, desc='Node''s Radii')
    #Reacts   =   Array(units=None, dtype=int,      desc='Node IDs for reactions + Fixity values (1/0)')  #Fixed for the time being with all fixity (6 dofs per node fixed)
    #nmems  =     Int(  units=None,                 desc='Total Number of Members in the Jacket, no Tower')
    #mems   =     Array(units=None, dtype=int,      desc='Member Connectivity Node i Node j for every member (i.e.,element)')
    #props  =     Array(            dtype=np.float, desc='Jacket Member''s xsectional and material properties')
    #XnsfM  =     Array(            dtype=np.float, desc='Jacket Member''s Global to Local DIrection Cosine Matrices [3,3,nelems]')
    #TubeObjs=    Instance(Klass=Tube,                  desc='Object of Class Tube for Jacket Elements. Tube Objects, one per element')
    #jnt_masses=  Array(            dtype=np.float, desc='Jacket Concentrated Masses at Joints: (joint no, mass, IXX, IYY, IZZ)')
    #jnt_masses_yaw=  Array(            dtype=np.float, desc='Jacket Concentrated Masses at Joints: (joint no, mass, IXX, IYY, IZZ) accounting for RNA yaw')

def TwrGeoOutputs(self):
    """Basic Geometric Outputs needed to build Tower of Jacket."""
    self.add_param('TwrGeoOutputs', , units=, pass_by_object=True)
    # TwrObj   =  Instance(Klass=Tube,                             desc='Object of Class Tube for Tower portion of Jacket')
    self.TwrGeoOutputs_TwrObj = Tube()      # Object of Class Tube for Tower portion of Jacket
    # Twr2RNAObj =Instance(Klass=RigidMember,                      desc='Object of Class RigidMember for Rigid Tower portion')
    self.TwrGeoOutputs_Twr2RNAObj = RigidMember()
    # joints   = Array(np.array([]),    units='m', dtype=np.float, desc='Pile Joint (with legs) Coordinates (3,nlegs=nfaces)')
    self.add_param('TwrGeoOutputs:joints', np.array([]), units='m', pass_by_object=True) # Pile Joint (with legs) Coordinates (3,nlegs=nfaces)
    # nodes    = Array(np.empty((3,10)),units='m', dtype=np.float, desc='Tower ALL Nodes'' Coordinates (3,nNodes)')
    self.add_param('TwrGeoOutputs:nodes', np.empty((3,10)), units='m', pass_by_object=True)  # Tower ALL Nodes'' Coordinates (3,nNodes)
    # nNodes   = Int(0,           units=None,                      desc='Number of Tower Nodes INCLUDING joint at TP')
    self.add_param('TwrGeoOutputs:nNodes', 0, pass_by_object=True)
    # nElems   = Int(0,           units=None,                      desc='Number of of elements in the flexible portion of tower')
    self.add_param('TwrGeoOutputs:nElems', 0, pass_by_object=True) # Number of of elements in the flexible portion of tower
    # mass =     Float(           units='kg',                      desc='Tower Mass')
    self.add_param('TwrGeoOutputs:mass', 0., units='kg',pass_by_object=True) # Tower Mass
    # TopMass  = Array(np.zeros([10]),            dtype=np.float, desc='Tower Top mass, Ixx, Iyy, Izz,Ixy,Ixz,Iyz,CMxoff,CMyoff,CMzoff from RNA properties in input')
    self.add_param('TwrGeoOutputs:TopMass', np.zeros([10]), pass_by_object=True) # Tower Top mass, Ixx, Iyy, Izz,Ixy,Ixz,Iyz,CMxoff,CMyoff,CMzoff from RNA properties in input
    # TopMass_yaw  = Array(np.zeros([10]),        dtype=np.float, desc='Tower Top mass, Ixx, Iyy, Izz,Ixy,Ixz,Iyz,CMxoff,CMyoff,CMzoff from RNA properties in input including yaw angle w.r.t. Global XYZ')
    self.add_param('TwrGeoOutputs:TopMass_yaw', np.zeros([10]), pass_by_object=True) # Tower Top mass, Ixx, Iyy, Izz,Ixy,Ixz,Iyz,CMxoff ,CMyoff,CMzoff from RNA properties in input including yaw angle w.r.t. Global XYZ
    # TwrlumpedMass=Array(np.zeros([1,11]),dtype=np.float, desc='Concentrated masses along tower: first column z''s from base of tower; 2nd through 11th column: mass and Ixx,Iyy,Izz,Ixy,Ixz,Iyz,CMxoff,CMyoff,CMzoff values')
    self.add_param('TwrGeoOutputs:TwrlumpedMass', np.zeros([1,11]), pass_by_object=True) # Concentrated masses along tower: first column z''s from base of tower; 2nd through 11th column: mass and Ixx,Iyy,Izz,Ixy,Ixz,Iyz,CMxoff,CMyoff,CMzoff values
    # HH       = Float(            units='m', desc='Hub-Height')
    self.add_param('TwrGeoOutputs:HH', 0., units='m', pass_by_object=True) # Hub-Height
    #Htwr     = Float(            units='m', desc='Tower Length')
    self.add_param('TwrGeoOutputs:Htwr', 0., units='m', pass_by_object=True) # Tower Length
    # Htwr2     = Float(            units='m', desc='Tower at constant cross section Length')
    self.add_param('TwrGeoOutputs:Htwr2', 0., units='m', pass_by_object=True) # Tower at constant cross section Length
    # Thoff_yaw= Array(np.zeros([3]),units='m',  dtype=np.float, desc='Tower-top to Hub-center vector in yawed coordinate system (TO BE MOVED SOMEWHERE ELSE at one point)')
    self.add_param('TwrGeoOutputs:Thoff_yaw', np.zeros([3]), units='m', pass_by_object=True) # Tower-top to Hub-center vector in yawed coordinate system (TO BE MOVED SOMEWHERE ELSE at one point)
    rna_yawedcm  = Array(np.zeros([3]),        dtype=np.float, desc='Tower Top mass CMxoff,CMyoff,CMzoff in yawed coordinate system(TO BE MOVED SOMEWHERE ELSE at one point)')
    self.add_param('TwrGeoOutputs:rna_yawedcm', np.zeros([3]), pass_by_object=True) # Tower Top mass CMxoff,CMyoff,CMzoff in yawed coordinate system(TO BE MOVED SOMEWHERE ELSE at one point)

def Frame3DDaux(self):
    """General Frame3DD parameters"""
    # sh_fg  = Int(1,  units=None, desc='Shear Deformation Flag: 1=Yes, 0=No.')
    self.add_param('TwrGeoOutputs:sh_fg', 1, pass_by_object=True) # Shear Deformation Flag: 1=Yes, 0=No.
    # geo_fg = Int(0,  units=None, desc='Geometric Stiffness Effects (for buckling analysis) Flag: 1=Yes, 0=No.')
    self.add_param('TwrGeoOutputs:geo_fg', 0, pass_by_object=True)  # Geometric Stiffness Effects (for buckling analysis) Flag: 1=Yes, 0=No.
    # exagg  = Int(1,  units=None, desc='Shear Deformation Flag: 1=Yes, 0=No.')
    self.add_param('TwrGeoOutputs:exagg', 1, pass_by_object=True)  # Shear Deformation Flag: 1=Yes, 0=No.
    # deltaz =Float(10., units='m',  desc='member z-axis increment for internal forces calc')
    self.add_param('TwrGeoOutputs:deltaz', 10., units='m', pass_by_object=True)  # member z-axis increment for internal forces calc
    # gvector=Array(np.array([0.,0,-9.8065]), units='m/s**2', desc='Inertial Frame Acceleration. For Gravity acceleration gz must be <0.')
    self.add_param('TwrGeoOutputs:gvector', np.array([0.,0,-9.8065]), units='m/s**2',pass_by_object=True)  # Inertial Frame Acceleration. For Gravity acceleration gz must be <0
    # nModes= Int(6,  units=None, desc='Number of desired dynamic modes (nModes)')
    self.add_param('TwrGeoOutputs:nModes', 6, pass_by_object=True)  # Number of desired dynamic modes (nModes)
    # nModesAn=Int(6, units=None, desc='Number of desired dynamic modes to animate(nModesAn)')
    self.add_param('TwrGeoOutputs:nModesAn', 6, pass_by_object=True)  # Number of desired dynamic modes to animate(nModesAn)
    # Mmethod=Int(1,  units=None, desc='Dynamic Eigenvalue Method: 1= Subspace-Nacobi iteration, 2= Stodola (matrix iteration) method')
    self.add_param('TwrGeoOutputs:Mmethod', 1, pass_by_object=True)  # Dynamic Eigenvalue Method: 1= Subspace-Nacobi iteration, 2= Stodola (matrix iteration) method'
    # lump=   Int(0,  units=None, desc='Flag: 0= consistent mass matrix, 1= lumped mass matrix0')
    self.add_param('TwrGeoOutputs:lump', 0, pass_by_object=True)  # Flag: 0= consistent mass matrix, 1= lumped mass matrix0
    # tol= Float(1e-9,units=None, desc='Frequency convergence tolerance')
    self.add_param('TwrGeoOutputs:tol', 1e-9, pass_by_object=True)  # Frequency convergence tolerance
    # shift= Float(0.,  units=None,  desc='frequency shift-factor for rigid body modes, make 0 for pos.def. [K]')
    self.add_param('TwrGeoOutputs:shift', 0., pass_by_object=True)  # frequency shift-factor for rigid body modes, make 0 for pos.def. [K]
    # exagg_modal=Float(10.,units=None, desc='Exaggeration factor for modal mesh deformations plotting')
    self.add_param('TwrGeoOutputs:exagg_modal', 10., pass_by_object=True)  # Exaggeration factor for modal mesh deformations plotting
    # pan=Float(2,units=None,desc='animation pan rate 0=no panning')
    self.add_param('TwrGeoOutputs:pan', 10., pass_by_object=True)  # animation pan rate 0=no panning
    Redux=Int(2,units=None,desc='matrix condensation method ...  0=none, 1=static, 2=Guyan, 3=dynamic')
    self.add_param('TwrGeoOutputs:Redux', 2, pass_by_object=True) # matrix condensation method ...  0=none, 1=static, 2=Guyan, 3=dynamic
def main():
    pass

if __name__ == '__main__':
    main()
