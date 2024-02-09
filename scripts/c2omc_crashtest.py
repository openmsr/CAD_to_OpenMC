#!/usr/bin/env python
import pathlib as pl
import subprocess as sp
import sys
import h5py
import numpy as np

import openmc

class h5mF(h5py.File):
    def __init(self,*args,**kwargs):
        super().__init__(args,**kwargs)
        taggroup=h5f['tstt']['tags']
        id_dataset=taggroup['NAME']['id_list']
        mat_dataset=taggroup['NAME']['values']
        cat_dataset=taggroup['CATEGORY']['values']
        cat_ids=taggroup['CATEGORY']['id_list']

class MeshPlot():
    def __init__(self,stepf='cad.step', init_openmc=True, pxl=2048, cyl=False, vol=False, origin=[0.0,0.0,0.0], width=[0.0,0.0,0.0], particles=100000, **kwargs):
        self.stepf=stepf
        self.h5mf=pl.Path(stepf).with_suffix('.h5m')
        self.pixels=pxl
        self.cyl=cyl
        self.vol=vol
        self.origin=np.array(origin)
        self.width=np.array(width)
        self.particles=int(float(particles))
        if (init_openmc):
            self.init_openmc()

    def mesh(self,tol=1e-3,angtol=0.01):
        self.amb.delete_intermediate_files=True
        self.amb.tolerance=tol
        self.amb.angular_tolerance=angtol
        self.amb.run(backend='stl2', h5m_filename=self.h5mf, merge=True)

    def init_openmc(self):
        self.bld_mats()
        self.bld_geom()
        self.bld_settings()
        if (self.vol):
            self.bld_plot3d()
        else:
            self.bld_plots()

    def plot(self):
        self.bld_geom()
        self.bld_plots()
        sp.run(['openmc','-p'])

    def plot3d(self):
        self.bld_geom()
        self.bld_plot3d()
        sp.run(['openmc','-p'])

    def bld_geom(self):
        du=openmc.DAGMCUniverse(self.h5mf, auto_geom_ids=True)
        bb=du.bounding_box
        print(f"bounding box of geometry: lower left={bb[0]}, upper_right={bb[1]}")
        if(self.cyl):
            rad=np.linalg.norm(0.5*(bb[1]-bb[0])[:2])
            ztop=openmc.ZPlane(z0=bb[1][2],boundary_type='vacuum')
            zbot=openmc.ZPlane(z0=bb[0][2],boundary_type='vacuum')
            bcyl=openmc.ZCylinder(r=rad,x0=0.5*(bb[1]+bb[0])[0],y0=0.5*(bb[1]+bb[0])[1],boundary_type='vacuum')
            uc=openmc.Cell(region=(-bcyl & -ztop & +zbot), fill=du)
        else:
            min_r=np.linalg.norm(0.5*(bb[1]-bb[0]))
            bndr=openmc.Sphere(r=min_r,boundary_type='vacuum')
            uc=openmc.Cell(region=-bndr, fill=du)
        root=openmc.Universe()
        root.add_cell(uc)
        self.geometry=openmc.Geometry(root)
        self.bb=bb
        if np.sum(self.origin)==0.0:
            self.origin=np.mean(bb,axis=0)
        if np.sum(self.width)==0.0:
            self.width=((bb[1]-bb[0]))
        self.geometry.export_to_xml()

    def retag_h5m_mats(self,retag: dict):
        h5mf(self.h5mf)
        mtags_in=[ bytes(v).decode().strip('\x00') for v in h5mf.mat_dataset ]
        mtags_out=mtags_in
        for k,v in retag:
            if k in mtags_in:
                h5mf.mat_dataset[k]=v

    def bld_mats(self):
        h5f=h5py.File(self.h5mf)
        taggroup=h5f['tstt']['tags']
        id_dataset=taggroup['NAME']['id_list']
        mat_dataset=taggroup['NAME']['values']
        #cat_dataset=taggroup['CATEGORY']['values']
        #cat_ids=taggroup['CATEGORY']['id_list']

        mats=[]
        ids=[]

        #for i,v,c in zip(id_dataset,mat_dataset,cat_dataset):
        for i,v in zip(id_dataset,mat_dataset):
          s=bytes(v).decode()
          ii=int(i)
          ids.append(ii)
          mats.append(s.strip('\x00'))
          ss=s.strip('\x00')
          #cc=bytes(c).decode().strip('\x00')

        matlist=openmc.Materials()
        for matraw in mats:
            mat=matraw.replace("mat:","")
            if mat.endswith('_comp'):
                mat=mat.replace('_comp','')
            if mat in [n.name for n in matlist]:
                continue
            m=openmc.Material(name=mat)
            m.add_nuclide('He3',1.0,'ao')
            matlist.append(m)
        matlist.export_to_xml()
        self.matlist=matlist

    def bld_plots(self):
        colors=['blue','green','magenta','red','cyan','steelblue']
        colordict={m:c for m,c in zip(self.matlist,colors)}
        bb=self.bb
        p1=openmc.Plot().from_geometry(self.geometry)
        p1.basis='xy'
        p1.width=self.width[:2]
        p1.origin=self.origin
        p1.color_by='material'
        p1.colors=colordict
        p1.pixels=(self.pixels,self.pixels)
        p1.filename=pl.Path(self.h5mf).stem+'_xy'

        p2=openmc.Plot().from_geometry(self.geometry)
        p2.basis='xz'
        p2.width=self.width[[0,2]]
        p2.origin=self.origin
        p2.color_by='material'
        p2.colors=colordict
        p2.pixels=(self.pixels,self.pixels)
        p2.filename=pl.Path(self.h5mf).stem+'_xz'

        p3=openmc.Plot().from_geometry(self.geometry)
        p3.basis='yz'
        p3.width=self.width[1:3]
        p3.origin=self.origin
        p3.color_by='material'
        p3.colors=colordict
        p3.pixels=(self.pixels,self.pixels)
        p3.filename=pl.Path(self.h5mf).stem+'_yz'
        plts=openmc.Plots([p1,p2,p3])
        plts.export_to_xml()

    def bld_plot3d(self):
        colors=['blue','green','magenta','red','cyan','steelblue']
        colordict={m:c for m,c in zip(self.matlist,colors)}
        pp=openmc.Plot().from_geometry(self.geometry)
        pp.type='voxel'
        pp.width=self.width
        pp.origin=self.origin
        pp.color_by='material'
        pp.colors=colordict
        pp.pixels=(self.pixels,self.pixels,self.pixels)
        pp.filename=pl.Path(self.h5mf).stem+'_xyz'
        plts=openmc.Plots([pp])
        plts.export_to_xml()


    def bld_settings(self):
        settings=openmc.Settings()
        settings.particles=self.particles
        settings.batches=10
        settings.run_mode='fixed source'
        source=self.bld_source()
        settings.source=self.bld_source()
        self.settings=settings
        settings.export_to_xml()

    def bld_source(self):
        source=openmc.Source()
        source.space=openmc.stats.Point(self.origin)
        source.angle=openmc.stats.Isotropic()
        source.energy=openmc.stats.Discrete([2e6],[1.0])
        return source

if __name__=='__main__':
    import argparse
    prs=argparse.ArgumentParser()
    prs.add_argument('stepfile')
    prs.add_argument('--pxl', type=int, default=256)
    prs.add_argument('--cyl', action='store_true')
    prs.add_argument('--vol', action='store_true')
    prs.add_argument('--origin',nargs=3,type=float, default=[0.0,0.0,0.0])
    prs.add_argument('--width',nargs=3,type=float, default=[0.0,0.0,0.0])
    args=prs.parse_args()
    try:
        mp=MeshPlot(args.stepfile, True, **vars(args))
    except KeyError:
        print('please supply as filename as first argument')
        exit(1)
    #mp.plot()
    #mp.plot3d()
