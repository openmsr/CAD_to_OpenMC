from initialize_materials import create_materials
import openmc
import matplotlib.colors as mcolors
from matplotlib import colormaps as cm
import sys
import math
import pathlib as pl
import numpy as np
import itertools as it
in2cm = 2.54

def create_materials(operating_temp):
    #fuel-like materials
    safety_block = openmc.Material(name='safety_block', temperature=operating_temp)
    safety_block.add_element('Mo',1.1399e-3,'ao')
    safety_block.add_nuclide('U233',4.6384e-6,'ao')
    safety_block.add_nuclide('U234',4.7068e-4,'ao')
    safety_block.add_nuclide('U235',4.2848e-2,'ao')
    safety_block.add_nuclide('U236',3.1140e-4,'ao')
    safety_block.add_nuclide('U238',2.3253e-3,'ao')
    safety_block.set_density('g/cc',7288/401.298)
    safety_block.volume=401.298

    ISP = openmc.Material(name='ISP', temperature=operating_temp)
    ISP.add_element('Mo',1.7248e-3,'ao')
    ISP.add_nuclide('U233',4.6326e-6,'ao')
    ISP.add_nuclide('U234',4.7007e-4,'ao')
    ISP.add_nuclide('U235',4.2791e-2,'ao')
    ISP.add_nuclide('U236',3.1101e-4,'ao')
    ISP.add_nuclide('U238',2.3249e-3,'ao')
    ISP.set_density('g/cc',(3511.0+3312.0)/386.869)
    ISP.volume=386.869
    subassembly_plate_inner = ISP.clone()
    subassembly_plate_inner.name = 'subassembly_plate_inner'

    fuel_ring = openmc.Material(name='fuel_ring', temperature=operating_temp)
    fuel_ring.add_element('Mo',1.7030e-3,'ao')
    fuel_ring.add_nuclide('U233',4.6343e-6,'ao')
    fuel_ring.add_nuclide('U234',4.7016e-4,'ao')
    fuel_ring.add_nuclide('U235',4.2801e-2,'ao')
    fuel_ring.add_nuclide('U236',3.1112e-4,'ao')
    fuel_ring.add_nuclide('U238',2.3328e-3,'ao')
    fuel_ring_mass=8378.0+8134.0+8046.0+6909.0+8166.0+9283.0
    fuel_ring.set_density('g/cc',fuel_ring_mass/2782.931)
    fuel_ring.volume=2782.931

    control_rod = openmc.Material(name='control_rod', temperature=operating_temp)
    control_rod.add_element('Mo',1.5550e-3,'ao')
    control_rod.add_nuclide('U233',4.7315e-6,'ao')
    control_rod.add_nuclide('U234',4.8008e-4,'ao')
    control_rod.add_nuclide('U235',4.3703e-2,'ao')
    control_rod.add_nuclide('U236',3.1765e-4,'ao')
    control_rod.add_nuclide('U238',2.3767e-3,'ao')
    control_rod.set_density('g/cc',824.0/44.88)
    control_rod.volume=44.88

    burst_rod = openmc.Material(name='burst_rod', temperature=operating_temp)
    burst_rod.add_element('Mo',1.5909e-3,'ao')
    burst_rod.add_nuclide('U233',4.8405e-6,'ao')
    burst_rod.add_nuclide('U234',4.9114e-4,'ao')
    burst_rod.add_nuclide('U235',4.4710e-2,'ao')
    burst_rod.add_nuclide('U236',3.2497e-4,'ao')
    burst_rod.add_nuclide('U238',2.4314e-3,'ao')
    burst_rod.set_density('g/cc',826.0/43.975)
    burst_rod.volume=43.975

    #structure materials
    ss303 = openmc.Material(name='ss303', temperature=operating_temp)
    ss303.add_element('C',0.075,'wo')
    ss303.add_element('Si',1.00,'wo')
    ss303.add_element('P',0.1,'wo')
    ss303.add_element('S',0.3,'wo')
    ss303.add_element('Cr',18.00,'wo')
    ss303.add_element('Mn',1.00,'wo')
    ss303.add_element('Fe',70.225,'wo')
    ss303.add_element('Ni',9.0,'wo')
    ss303.add_element('Mo',0.30,'wo')
    ss303.set_density('g/cc',8.0)
    spindle = ss303.clone()
    spindle.name='spindle'
    safety_block_base = ss303.clone()
    safety_block_base.name='safety_block_base'
    clamp_support = ss303.clone()
    clamp_support.name='clamp_support'

    vascomax300 = openmc.Material(name='vascomax300', temperature=operating_temp)
    vascomax300.add_element('C',0.02,'wo')
    vascomax300.add_element('Al',0.1,'wo')
    vascomax300.add_element('Si',0.05,'wo')
    vascomax300.add_element('P',0.005,'wo')
    vascomax300.add_element('S',0.005,'wo')
    vascomax300.add_element('Ti',0.730,'wo')
    vascomax300.add_element('Mn',0.05,'wo')
    vascomax300.add_element('Fe',66.94,'wo')
    vascomax300.add_element('Co',8.8,'wo')
    vascomax300.add_element('Ni',18.5,'wo')
    vascomax300.add_element('Mo',4.8,'wo')
    vascomax300.set_density('g/cc',8.0)
    clamp = vascomax300.clone()
    clamp.name='clamp'

    sae4340 = openmc.Material(name='sae4340', temperature=operating_temp)
    sae4340.add_element('C',0.405,'wo')
    sae4340.add_element('Si',0.225,'wo')
    sae4340.add_element('P',0.018,'wo')
    sae4340.add_element('S',0.02,'wo')
    sae4340.add_element('Cr',0.8,'wo')
    sae4340.add_element('Mn',0.725,'wo')
    sae4340.add_element('Fe',95.757,'wo')
    sae4340.add_element('Ni',1.8,'wo')
    sae4340.add_element('Mo',0.250,'wo')
    sae4340.set_density('g/cc',7.85)
    subassembly_cover_plate = sae4340.clone()
    subassembly_cover_plate.name='subassembly_cover_plate'
    support_pad_ring = sae4340.clone()
    support_pad_ring.name='support_pad_ring'
    bearing_ring = sae4340.clone()
    bearing_ring.name='bearing_ring'

    aluminium = openmc.Material(name='aluminium', temperature=operating_temp)
    aluminium.add_element('Mg',1.0,'wo')
    aluminium.add_element('Al',97.23,'wo')
    aluminium.add_element('Si',0.6,'wo')
    aluminium.add_element('Ti',0.075,'wo')
    aluminium.add_element('Cr',0.5*(0.04+0.35),'wo')
    aluminium.add_element('Mn',0.075,'wo')
    aluminium.add_element('Fe',0.35,'wo')
    aluminium.add_element('Cu',0.5*(0.15+0.4),'wo')
    aluminium.add_element('Zn',0.125,'wo')
    aluminium.set_density('g/cc',2.7)
    mounting_plate = aluminium.clone()
    mounting_plate.name='mounting_plate'

    helium=openmc.Material(name='helium', temperature=operating_temp)
    helium.add_nuclide('He3',1.0,'ao')
    helium.set_density('g/cc',1.23456e-5)

    mats = [fuel_ring,safety_block,control_rod,burst_rod,ISP,subassembly_plate_inner,ss303,
        spindle, safety_block_base,clamp_support, clamp, vascomax300, subassembly_cover_plate,support_pad_ring, bearing_ring,sae4340,mounting_plate,aluminium, helium]
    return mats

class GIV_reactor:

    #           CASE#    BRz          CR1z      CR2z
    def __init__(self, core_h5m, burst_rod_h5m,control_rod_h5m, CAD_rods=True,
            fuel=None,operating_temp=273.15,batches=20,inactive=5,particles=1000,
            padding=0, burst_rod_z=0, ctrl_rod_z=0, output_dir='.', case=None, verbose=False):

        for ff in [core_h5m, burst_rod_h5m, control_rod_h5m]:
            if (not pl.Path(ff).exists()):
                print(f'That file \"{ff}\" does not appear to exist')
                raise IOError
        self.fuel=fuel

        self.core_h5m = core_h5m
        self.burst_rod_h5m = burst_rod_h5m
        self.control_rod_h5m = control_rod_h5m

        self.padding=padding

        try:
            iter(ctrl_rod_z)
            self.ctrl_rod_z=shim_rod_z
        except TypeError as e:
            self.CRz=[ctrl_rod_z]*2
        self.BRz=burst_rod_z
        self.verbose=verbose

        self.case=case
        #possibly override by use of one the predfined cases
        if case is not None and case in range(1,6):
            self.Brz, *self.CRz = self.experiment_cases()[case]

        self.operating_temp=operating_temp
        self.batches=batches
        self.particles=particles
        self.inactive=inactive
        self.vol_trig=None

        self.model=None
        self.geometry=None
        self.materials=None
        self.settings=None
        self.tallies=None
        self.plots=None

        self.output_dir=output_dir

    def experiment_cases(self):
        a = 17.3757 #mounting plate to origin distance
        c1_fi = 7.712*in2cm #c-rod 1 fully inserted, top to mounting plate
        c2_fi = 7.239*in2cm #c-rod 2 fully inserted, top to mounting plate
        br_fi = 7.272*in2cm #b-rod fully inserted, top top mounting plate
        br_fw = 4.302*in2cm #b-rod fully withdrawn, top top mounting plate
        #cr positions are given in table 5 as inches withdrawn
        cases={
            1: ( br_fi-a, c1_fi-a-4.001*in2cm, c2_fi-a-0.449*in2cm ),
            2: ( br_fi-a, c1_fi-a-1.998*in2cm, c2_fi-a-1.666*in2cm ),
            3: ( br_fi-a, c1_fi-a-0.493*in2cm, c2_fi-a-3.794*in2cm ),
            4: ( br_fw-a, c1_fi-a-0.469*in2cm, c2_fi-a-0.447*in2cm ),
            5: ( br_fi-a, c1_fi-a-0.319*in2cm, c2_fi-a-0.656*in2cm ),
        }
        if self.verbose:
            print('# case#    BR z    CR1_z    CR2_z')
            for i in cases.keys():
                print(i,cases[i])

        return cases

    def bld_materials(self):
        mats=create_materials(self.operating_temp)
        self.materials=openmc.Materials(mats)


    def control_rod_csg_univ(self):
        cr_mat = [m for m in self.materials if m.name=='control_rod']
        helium = [m for m in self.materials if m.name=='helium']
        co = openmc.ZCylinder(r0=2.1844/2.0)
        ci = openmc.ZCylinder(r0=0.9525/2.0)
        tpo = openmc.ZPlane(z0=0)
        tpi = openmc.ZPlane(z0=-1.905)
        bpo = openmc.ZPlane(z0=-12.70)
        bpi = openmc.ZPlane(z0=-10.795)
        c = openmc.Cell( fill=cr_mat, region =( (+ci & -co & -top & +bpo) | (-ci & -tpi & +bpi) ) )
        c_prime = openmc.Cell( fill=helium, region=~c.region )
        univ = openmc.Universe([c,c_prime])
        return univ

    def burst_rod_csg_univ(self):
        br_mat = [m for m in self.materials if m.name=='burst_rod']
        helium = [m for m in self.materials if m.name=='helium']
        co = openmc.ZCylinder(r0=2.1844/2.0)
        ci = openmc.ZCylinder(r0=0.9525/2.0)
        tpo = openmc.ZPlane(z0=0)
        tpi = openmc.ZPlane(z0=-3.175)
        bpo = openmc.ZPlane(z0=-12.70)
        bpi = openmc.ZPlane(z0=-10.795)
        c = openmc.Cell( fill=cr_mat, region =( (+ci & -co & -top & +bpo) | (-ci & -tpi & +bpi) ) )
        c_prime = openmc.Cell( fill=helium, region=~c.region )
        univ = openmc.Universe([c,c_prime])
        return univ

    def bld_geometry(self):
        if self.materials is None:
            self.bld_materials()
        helium=[m for m in self.materials if m.name=='helium']
        #import dagmc-geometry defining the reactor core
        du=openmc.DAGMCUniverse(self.core_h5m, auto_geom_ids=True)


        #rod universe geometries
        control_rod_du = openmc.DAGMCUniverse(self.control_rod_h5m, auto_geom_ids=True)
        burst_rod_du = openmc.DAGMCUniverse(self.burst_rod_h5m, auto_geom_ids=True)

        dr=6.6675

        BRz, CR1z, CR2z = self.BRz, *self.CRz

        posBR=np.array( (-dr, 0.0, BRz-1e-2) )
        posCR1=np.array( (dr*math.cos(math.pi/3),dr*math.sin(math.pi/3), CR1z) )
        posCR2=np.array( (dr*math.cos(-math.pi/3),dr*math.sin(-math.pi/3), CR2z) )

        BR_cyl=openmc.ZCylinder(r=2.34/2.0,x0=posBR[0],y0=posBR[1])
        CR1_cyl=openmc.ZCylinder(r=2.34/2.0,x0=posCR1[0],y0=posCR1[1])
        CR2_cyl=openmc.ZCylinder(r=2.34/2.0,x0=posCR2[0],y0=posCR2[1])

        BR_upper=openmc.ZPlane(z0=posBR[2])
        CR1_upper=openmc.ZPlane(z0=posCR1[2])
        CR2_upper=openmc.ZPlane(z0=posCR2[2])

        #find lowest z-coordinate
        bottomz = np.min( [U.bounding_box.lower_left[2] for U in [du,burst_rod_du,control_rod_du]])
        rod_lower = openmc.ZPlane(z0=bottomz)
        rod_cyl_top = openmc.ZPlane(z0=4.64820)

        BR=openmc.Cell(region=(-BR_cyl & -rod_cyl_top & +rod_lower), fill=burst_rod_du)
        BR.translation=posBR
        #BR_void=openmc.Cell(region=(-BR_cyl & +BR_upper & -rod_cyl_top), fill=helium)

        CR1=openmc.Cell(region=(-CR1_cyl & -rod_cyl_top & +rod_lower), fill=control_rod_du)
        CR1.translation=posCR1
        #CR1_void=openmc.Cell(region=(-CR1_cyl & +CR1_upper & -rod_cyl_top), fill=helium)

        CR2=openmc.Cell(region=(-CR2_cyl & -rod_cyl_top & +rod_lower), fill=control_rod_du)
        CR2.translation=posCR2
        #CR2_void=openmc.Cell(region=(-CR2_cyl & +CR2_upper & -rod_cyl_top), fill=helium)

        #core cell definition
        #core_region=self.bc & ~BR.region & ~BR_void.region & ~CR1.region & CR1_void.region & ~CR2.region & ~CR2_void.region
        self.bc=du.bounding_region(padding_distance=5)

        model_region=self.bc & ~BR.region & ~CR1.region & ~CR2.region
        scene = openmc.Cell(region=model_region, fill=du)

        #define geometry
        root = openmc.Universe()
        #root.add_cells([core,BR,BR_void,CR1,CR1_void,CR2, CR2_void])
        root.add_cells([scene,BR,CR1,CR2])

        self.geometry = openmc.Geometry(root)

    def bld_settings(self):
        s=openmc.Settings(particles=int(self.particles), batches=self.batches)
        s.inactive=self.inactive
        s.rel_max_lost_particles=0.01
        s.source=self.bld_source()
        s.temperature={'method':'interpolation', 'range':[250,350]}
        self.settings=s

    def bld_source(self):
        src = openmc.Source()
        src.space = openmc.stats.Box(lower_left=[-10,-10,-8], upper_right=[10,10,8])
        src.angle = openmc.stats.Isotropic()
        src.energy = openmc.stats.Watt()
        return src

    def add_meshxy_tally(self, dims=[256,256,1],ll=None, ur=None, name='xy_nflux', scores=['flux','fission']):
        if ll is None:
            ll=self.geometry.bounding_box.lower_left
        if ur is None:
            ur=self.geometry.bounding_box.upper_right
        meshxy=openmc.RegularMesh()
        meshxy.n_dimensions=2
        meshxy.dimension=dims
        meshxy.lower_left=ll
        meshxy.upper_right=ur
        flt=openmc.MeshFilter(meshxy)
        npf=openmc.ParticleFilter('neutron')
        tally=openmc.Tally(name=name)
        tally.filters=[flt,npf]
        tally.scores=scores
        return tally

    def add_material_tally(self,ll=None, ur=None, name='mat',scores=['flux']):
        if ll is None:
            ll=self.geometry.bounding_box.lower_left
        if ur is None:
            ur=self.geometry.bounding_box.upper_right
        matfilter = openmc.MaterialFilter(self.materials)
        tally=openmc.Tally(name=name)
        tally.filters=[matfilter]
        tally.scores=scores
        return tally

    def add_energy_tally(self,erange=None, name='energy'):
        if erange is None:
            ebins=np.logspace(1,7,200)
        ef=openmc.EnergyFilter(ebins)
        tally=openmc.Tally(name=name)
        tally.filters=[ef]
        tally.scores=['flux']
        return tally

    def bld_tallies(self):
        #if not specified add default tallies
        if self.geometry is None:
            self.bld_geometry()
        print(self.geometry.bounding_box)
        if self.tallies is None:
            self.tallies = openmc.Tallies()
            self.tallies.append(self.add_meshxy_tally())
            self.tallies.append(self.add_meshxy_tally(dims=[256,1,256], name='xz_nflux'))
            self.tallies.append(self.add_meshxy_tally(dims=[1,256,256], name='yz_nflux'))
            self.tallies.append(self.add_material_tally(name='heating', scores=['heating']))
            self.tallies.append(self.add_energy_tally())

    def bld_plots(self):
        if self.geometry is None:
            self.bld_geometry()
        #if not specified add default plots
        if self.plots is None:
            colordict={m:c for m,c in zip(self.materials,it.cycle(cm['Paired'].colors))}
            if self.verbose:
                for k,v in colordict.items():
                    print(k.name,v)

            bb=self.geometry.bounding_box
            width=np.array([300,300,300])
            pixels=1024
            p1=openmc.Plot().from_geometry(self.geometry)
            p1.basis='xy'
            p1.width=width[:2]
            p1.origin=(0,0,0)
            p1.color_by='material'
            #p1.colors=colordict
            p1.pixels=(pixels,pixels)
            p1.filename='reactor_xy'

            p2=openmc.Plot().from_geometry(self.geometry)
            p2.basis='xz'
            p2.width=width[[0,2]]
            p2.origin=(0,0,0)
            p2.color_by='material'
            #p2.colors=colordict
            p2.pixels=(pixels,pixels)
            p2.filename='reactor_xz'

            p3=openmc.Plot().from_geometry(self.geometry)
            p3.basis='yz'
            p3.width=width[1:3]
            p3.origin=(0,0,0)
            p3.color_by='material'
            #p3.colors=colordict
            p3.pixels=(pixels,pixels)
            p3.filename='reactor_yz'
            self.plots=openmc.Plots([p1,p2,p3])

    def bld_model(self):
        if self.materials is None:
            self.bld_materials()
        if self.geometry is None:
            self.bld_geometry()
        if self.tallies is None:
            self.bld_tallies()
        if self.settings is None:
            self.bld_settings()
        self.model = openmc.Model(geometry=self.geometry, materials=self.materials, tallies=self.tallies, plots=self.plots, settings=self.settings)
        return self.model

    def keff_run(self):
        if self.model is None:
            self.bld_model()
        openmc.lib.init()

    def export_to_xml(self,cwd=None):
        if self.model is None:
            self.bld_model()
        if cwd is not None:
            pp=pl.Path(cwd)
            pp.mkdir(parents=True, exist_ok=True)
            self.model.cwd=cwd
        self.model.export_to_model_xml()

    def run_volume_calc(self,trigger=None,cwd=None):
        if self.model is None:
            self.bld_model()
        if cwd is not None:
            pp=pl.Path(cwd)
            pp.mkdir(parents=True, exist_ok=True)
            self.model.cwd=cwd
        bb=self.bc.bounding_box
        vol_calc = openmc.VolumeCalculation(self.materials, self.particles, lower_left=bb[0], upper_right=bb[1])
        if (not trigger is None):
            vol_calc.set_trigger(trigger,'std_dev')
        self.model.settings.volume_calculations = [vol_calc]
        self.model.settings.run_mode='volume'
        self.model.export_to_model_xml()
        openmc.run()
        exit(0)

if __name__=='__main__':
    import argparse

    ap=argparse.ArgumentParser(prog='build_model')
    ap.add_argument('--case','--c',default=1,type=int,help="Experimental case to use (1..5)")
    ap.add_argument('--core',default='GIV_core_simplified.h5m',help='h5m file to use for the core')
    ap.add_argument('--BR',default='GIV_BR_simplified.h5m',help='h5m file to use for the burst rod')
    ap.add_argument('--CR',default='GIV_CR_simplified.h5m',help='h5m file to use for the control rod')
    ap.add_argument('--particles','-p',default=1000,help='number of particles per batch')
    ap.add_argument('--batches','-b',default=10,type=int,help='number of batches')
    ap.add_argument('--inactive','-i',default=2,type=int,help='inactive batches')
    ap.add_argument('--vol',action='store_true', help='if given perform a volume calculation.')
    ap.add_argument('--volprec', nargs='?',type=float, const=0.01, help='the wanted precision for volume calculations')
    args = ap.parse_args()
    GIV=GIV_reactor(args.core,args.BR,args.CR, case=args.case,
                    particles=args.particles,batches=args.batches, inactive=args.inactive, verbose=True)
    if (args.vol):
        GIV.run_volume_calc(args.volprec)
    GIV.bld_plots()
    GIV.export_to_xml()
