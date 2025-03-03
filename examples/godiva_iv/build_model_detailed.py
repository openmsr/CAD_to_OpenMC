import openmc
import matplotlib.colors as mcolors
from matplotlib import colormaps as cm
import sys
import math
import pathlib as pl
import numpy as np
import itertools as it
in2cm = 2.54

fuels={
"fuel_ring_101": ( ('U233',4.6624e-6),('U234',4.7311e-4),('U235',4.3069e-2),('U236', 3.1300e-4),('U238',2.3373e-3),('Mo',1.7013e-3) ),
"fuel_ring_102": ( ('U233',4.6953e-6),('U234',4.7635e-4),('U235',4.3364e-2),('U236', 3.1521e-4),('U238',2.3631e-3),('Mo',1.7366e-3) ),
"fuel_ring_103": ( ('U233',4.6795e-6),('U234',4.7439e-4),('U235',4.3186e-2),('U236', 3.1416e-4),('U238',2.3876e-3),('Mo',1.7209e-3) ),
"fuel_ring_104": ( ('U233',4.6634e-6),('U234',4.7327e-4),('U235',4.3083e-2),('U236', 3.1308e-4),('U238',2.3332e-3),('Mo',1.6974e-3) ),
"fuel_ring_105": ( ('U233',4.6427e-6),('U234',4.7116e-4),('U235',4.2892e-2),('U236', 3.1168e-4),('U238',2.3228e-3),('Mo',1.7104e-3) ),
"fuel_ring_106": ( ('U233',4.6903e-6),('U234',4.7581e-4),('U235',4.3313e-2),('U236', 3.1488e-4),('U238',2.3652e-3),('Mo',1.7319e-3) ),
"up_in_subassy_plate" : ( ('U233',4.7059e-6),('U234',4.7743e-4),('U235',4.3462e-2),('U236', 3.1593e-4),('U238',2.3684e-3),('Mo',1.7319e-3) ),
"int_in_subassy_plate" : ( ('U233',4.5562e-6),('U234',4.6239e-4),('U235',4.2092e-2),('U236', 3.0587e-4),('U238',2.2795e-3),('Mo',1.7319e-3) ),
"safety_block" : ( ('U233',4.6411e-6),('U234',4.7095e-4),('U235',4.2873e-2),('U236', 3.1158e-4),('U238',2.3266e-3),('Mo',1.7319e-3) ),
"control_rod" : ( ('U233',4.6944e-6),('U234',4.7632e-4),('U235',4.3360e-2),('U236', 3.1516e-4),('U238',2.3580e-3),('Mo',1.7319e-3) ),
"burst_rod" : ( ('U233',4.7807e-6),('U234',4.8507e-4),('U235',4.4157e-2),('U236', 3.2095e-4),('U238',2.4013e-3),('Mo',1.7319e-3) ),
}

steel303_composition=( ('C',3.0083e-4),('Si',1.7154e-3),('P',1.5554e-4,),('S',4.5067e-4),('Cr',1.6678e-2),('Mn',8.7693e-4),('Fe',6.0580e-2),('Ni',7.3878e-3),('Mo',1.5065e-4) )
s303_bits={ nm:steel303_composition for nm in ["spindle","spindle_nut","safety_block_base" ,"clamp_support","long_align_pin", "belly_band", "belly_band_bolt", "locking_bolt", "steel_screen", "steel_invertedt", "steel"]}

vascomax300_composition=(
('C',8.0221e-5),('Al',1.7855e-4),('Si',8.5768e-5),('P',7.7770e-6),('S',7.5112e-6),('Ti',7.3453e-4),('Mn',4.3847e-5),
('Fe',5.7746e-2),('Co',7.1938e-3),('Ni',1.5186e-2),('Mo',2.4103e-3) )

vmax_bits={ nm:vascomax300_composition for nm in ["clamp"]}

SAE_4340_composition=(
('C',1.5940e-3),
('Si',3.7872e-4),
('P',2.7472e-5),
('S',2.9481e-5),
('Cr',7.2734e-4),
('Mn',6.2385e-4),
('Fe',8.1057e-2),
('Ni',1.4499e-3),
('Mo',1.2319e-4),
)
sae_bits={nm:SAE_4340_composition for nm in ["bearing_ring", "sub_assy_cover_plate", "support_pad","burst_rod_pin"]}

AISI_1019_Carbon_Steel_composition=(
('C',6.7080e-4),
('P',3.0603e-5),
('S',3.6946e-5),
('Mn',7.3328e-4),
('Fe',8.3960e-2),
)

aisi_bits={nm:AISI_1019_Carbon_Steel_composition for nm in ["steel_bottom_cont_shield"]}

al_6061_t6_composition=(
('Mg',6.6898e-4),
('Al',5.8593e-2),
('Si',3.4736e-4),
('Ti',2.5469e-5),
('Cr',6.0978e-5),
('Mn',2.2197e-5),
('Fe',1.0190e-4),
('Cu',7.0365e-5),
('Zn',3.1082e-5),
)

al_6061_bits = {nm:al_6061_t6_composition for nm in ["mounting_plate", "core_cover" ]}

lexguard_plastic_composition=(
('C',75.48),
('H',5.55),
('O',18.87)
)
lexguard = {nm:lexguard_plastic_composition for nm in [ "plastic_cont_shield" ]}


def create_materials(temp):
    ms=openmc.Materials()
    for grp in [vmax_bits,aisi_bits,al_6061_bits,fuels,sae_bits,s303_bits,lexguard]:
        for k,v in grp.items():
            print(k)
            m=openmc.Material(name=k, temperature=temp)
            rho=0
            if "plastic" in k:
                for nuc in v:
                    m.add_element(nuc[0],nuc[1],'wo')
                m.set_density('g/cc',1.2)
                ms.append(m)
                continue
            for nuc in v:
                rho+=nuc[1]
                if nuc[0][-1] in set('0123456789'):
                    m.add_nuclide(nuc[0],nuc[1],'ao')
                else:
                    m.add_element(nuc[0],nuc[1],'ao')
            m.set_density('atom/b-cm',rho)
            ms.append(m)
    return ms

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

class GIV_reactor_wrod(GIV_reactor):

    def bld_geometry(self):
        if self.materials is None:
            self.bld_materials()
        helium=[m for m in self.materials if m.name=='helium']
        #import dagmc-geometry defining the reactor core
        case_h5m = self.core_h5m.replace('core.h5m',f'case{self.case}.h5m')
        print(f'going for this input file: {case_h5m}')
        du=openmc.DAGMCUniverse(case_h5m, auto_geom_ids=True)

        self.bc=du.bounding_region()

        #core cell definition
        core_region=self.bc
        core = openmc.Cell(region=core_region, fill=du)

        #define geometry
        root = openmc.Universe()
        root.add_cells([core])

        self.geometry = openmc.Geometry(root)

if __name__=='__main__':
    import argparse

    ap=argparse.ArgumentParser(prog='build_model')
    ap.add_argument('--case','--c',default=1,type=int,help="Experimental case to use (1..5)")
    ap.add_argument('--core',default='GIV_core_detailed.h5m',help='h5m file to use for the core')
    ap.add_argument('--BR',default='GIV_BR_detailed.h5m',help='h5m file to use for the burst rod')
    ap.add_argument('--CR',default='GIV_CR_detailed.h5m',help='h5m file to use for the control rod')
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
