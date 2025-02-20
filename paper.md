---
title: "CAD_to_OpenMC: from CAD design to particle transport"
tags:
    - Python
    - Monte Carlo
    - Model preparation
    - CAD
    - Mesh

authors:
    - name: Erik Bergbäck Knudsen
      orcid: 0000-0003-1827-9188
      corresponding: true
      affiliation: 1
      equal contribution: true
    - name: Lorenzo Chierici
      affiliation: 2

affiliations:
    - name: United Neux, Denmark
      index: 1
    - name: Copenhagen Atomics A/S, Denmark
      index: 2

date: XX September 2024
bibliography: bibliography.bib
---

# Summary
<!-- Here's a citation @knudsen2013mcxtrace , or in a different format [@knudsen2013mcxtrace].-->
We present CAD_to_OpenMC - a practical tool for allowing geometries designed in
CAD-systems to be used in Monte Carlo particle transport simulation systems.
Written in python, and easily installable through pip and/or conda channels,
CAD_to_openMC does so by importing a geometry in the form of a step-file, a
common export format among many CAD-systems, identifying objects within that
file, and running a triangularization algorithm on all surfaces contained within
the geometry.
The geometry is written to a HDF5-file in a format defined by the
DAGMC-library, which is usable by the main particle-transport packages.


# Statement of need
Most Monte Carlo neutronics-, and general particle transport (e.g. OpenMC, GEANT4, MCNP, Fluka [@openmc_2013;
@geant4_2003; @mcnp_2022; @fluka_2005; @fluka_2014]) software packages use the
concept of Constructive Solid Geometry (CSG) to the describe the geometry or
scene in which a computations are to performed. Computations generally amount
to tracing a large number of rays through the geometry. Such models range in
complexity from simple scenes of a only a few geometrical objects to several
million, but are generally considered tricky to work with - in that all object
boundaries must be constructed from analytical descriptions of geometrical
primaries (cylinders , planes, 2nd order surfaces, etc...).
Most engineering design, however, is performed using CAD-tools, which are
optimized for ease of design.

These facts form a compelling argument for bridging the gap between these
worlds. So far a number of solutions have been proposed which can mostly be put
two categories.

1. tools that convert CAD-models into CSG [@du_inversecsg_2018; @catalan_geouned_2024], and
2. tools that allow ray-tracing on CAD-geometries, or close derivatives thereof.

In the second category, one of the most wide-spread open source solutions is DAGMC [@wilson_acceleration_2010],
a library with couplings to all of OpenMC, MCNP, fluka, and Geant4.
DAGMC is proven to be fast and reliable, yet one hurdle remains: DAGMC
requires a discretized (meshed) approximation of the CAD-geometry. In general
not a trivial requirement.

CAD-software packages are for the most
part proprietary packages with little to no focus on being interchangeable. One
upshot of this is a notable lack of overlapping standard formatting between
packages, the exception being the STEP-format which is backed by an
international standard [@stepStd_2024]. Thus, a tool should support this format, if
it is not to be tied to any particular CAD-engine.

The CAD_to_OpenMC tool eases the process of generating a meshed description
of a CAD-generated geometry (in the form of a step-file) ready for inclusion in
transport codes through DAGMC.

# Method
CAD_to_OpenMC uses OCCT [@occt3d] to interact with CAD-geometries
and its' handling of step-files. Once the geometry has been imported, one of
several meshing back-ends may be called to create a discretized version of the
geometry. In the end the so generated geometry is exported into the DAGMC-expected format.

## Triangularization / Surface Meshing
<!-- Describe the imprinting and why it is very important especially for this case -->
When a curved geometry is to be discretized the result is an approximation.
This leads to potential problems with overlaps between regions if

a. objects are close to each other, such as for a cylindrical can with thin walls
b. objects share a surface such as the case for a liquid inside a can.
c. when objects have surfaces that touch each other.

The first case puts a constraint on the absolute tolerance of the
discretization, i.e. triangles have to be small enough not to cause crossing
surfaces.

In terms of case b, the discretization process must take care not to
re-evaluate any surface that is shared. Instead it must reuse the triangles of
the former instance. This means that objects cannot be independently processed.
The whole process is generally handled by the underlying geometry engine [@occt3d] through a
hierarchical model, in which objects consist of a set of surfaces,
formed by a set of curves, which themselves consist of line segments connecting points.
Object simply store links to surfaces
and their triangularizations (similar for curves --- line segments).
CAD_to_OpenMC handles this by generating a hash-code for each surface upon
processing. Each time a surface is encountered, it is examined, and if the
surface has already been encountered, it will be re-used, not recomputed.

In the latter case (c), imprinting has to be performed. This is the process
where the boundary curves of the smaller surface is projected onto the larger,
splitting the larger into two or more sub surfaces.

<!-- meshing backends -->
### Meshing backends
CAD_to_OpenMC is constructed such that it is flexible in terms of the
algorithms used to generate triangularized surfaces that approximate the
geometries, known as back-ends.
The present release supports a set of back-ends: {'gmsh','stl','stl2','db'}. In
most cases we have encountered 'stl2' performs well. It uses a very basic
algorithm, with no restriction on triangle aspect ratio, but also produces
discretized geometry-files of moderate size. In practice the aspect ratio does
not appear to be a problem for particle transport. In case the same
discretization is to be used for other types of application care should be
taken.

In some cases, the simple algorithm fails to create a watertight model. In such
cases we recommend using the 'db' back-end if available. The 'gmsh'-backend consistently
produces practical model, but often has the problem that it is memory hungry. E.g.
neither the ARE, nor the MSRE (see below) models could be run this way on our
available hardware (64GB workstation).


### Material tags
<!-- The way material tags are extracted here -->
If wanted CAD_to_OpenMC is able to use the CAD-generated part names, as material tags.
The default behavior is to use the first part of the part name as material
tag. This may be changed by supplying a hash-table, i.e. python dictionary, as
the tags-argument. Here the keys to the table are interpreted as regular
expressions and the values are taken to be material tags to use. Parts, whose
names do not match any of the regexp-keys, are tagged as vacuum by default. If
however a flags is set, then the material tag is extracted from the part name,
as if no list of tags is supplied. This allows to re-tag only a subset of the
parts.
The below example shows how to tag all parts with the name "wall" in them with "concrete"
and all parts ending with bellows with "steel".
```python
tag_dict = {".*bellows" : "steel", ".*wall.*" : "concrete"}
```

<!--
here's a test equation:
$$x= \int_a^{b+23} \frac{1}{1+gx} \mathrm{d}x \label{eq:integral}$$
which is nice

I am now going to refer to the equation: \autoref{eq:integral}.
-->



## Implicit complement
CAD_to_OpenMC can add an implicit complement to the output file. That is, the material tag
that will be applied to any part of the geometry which is *not* claimed by any
CAD part. This is done by simply assigning the name of the material as a string
to the implicit complement attribute of the base Assembly object.
If the attribute is set C2O assigns the extra material tag (with a suffix of \'\_comp\')
to the last part handled. That material tag gets picked up by the DAGMC-system and is used for any unclaimed volume.

# Results
We have chosen 3 reactor models as test systems. A table-top reactor and two full-scale molten salt reactors. The former
(GODIVA IV) model is included in the ICSBEP-benchmark project [@icsbep_2020] as case HEU-MET-FAST-086, the latter two were part of the motel salt reactor program at ORNL.

## GODIVA IV
This model was chosen since it exhibits moderate complexity and has a generous set of experimental data to benchmark against.
It is detailed enough to be cumbersome to model quickly using CSG but for with
CAD. Further, CSG models (for MCNP) may be found in the published benchmark,
along with detailed experimental validation data allowing for practical
comparison.
The reactor consists of a cylindrical core geometry (\autoref{fig:GIV_CAD}) held by 3 sets of clamps set
at 120 deg. offset. Additionally the core has three vertical holes (similarly
120 deg. apart), into which control and burst rods may be inserted from below
by vertical actuators. The rods themselves are similar in composition to the
fuel elements.

The benchmark includes 5 cases, which differ in terms of control- and burst-rod positions.
Also there are 3 geometries described:
1. detailed model which includes as built things (curved clamp etc.),
2. simplified model, where all surfaces are along principal axis, all corners 90 deg. etc., and
3. cylindrical symmetric model, to allow 2D-computations.
The benchmark reports only experimental results for the 2 first ones, but contains MCNP-geometries for
all 3 in [@icsbep_2020]. The 3rd model is created in order to use some legacy analysis tool which are purely 1D/2D.

Corresponding to \autoref{fig:GIV_CAD}, \autoref{fig:GIV_meshed} show the discretized version of the two reactor models used for our further analysis.

![CAD-drawings of the Godiva IV reactor, detailed (left) and simplified (right) versions. Note the rectangular clamps and supports in the simplified version. Visible are also the set of control and burst rods.\label{fig:GIV_CAD}](figs/GIV_both.pdf){#GIV_CAD width="50%"}

![Discretized Godiva IV-models, detailed (left) and simplified (right) versions.\label{fig:GIV_meshed}](figs/both_meshed.pdf){#GIV_meshed width=50%}

![Part-by-part comparison between volume calculations using stochastic volume estimation in OpenMC compared with direct volume calculations reported by CAD-software for the detailed model (green) and simplified benchmark model (magneta). The indictated intervals are the compuational error margins, almost exclusively stemming from the estimated error in the stochastic volume computation.\label{fig:voldiff}](figs/both_rel_voldiff.png)

Figure \ref{fig:voldiff} shows
differences in volumes between the discretized models and the exact CAD-model
for the various objects making up the parts. Volumes have been calculated using
a built-in feature of our CAD-package, whereas volumes from the discretized models have
been computed using the stochastic volume computation feature of OpenMC
[@openmc_2013]. In the latter case volumes are being computed by sampling a
number of points within a set boundary, while recording volume each point falls
within. The precision of the algorithm is directly governed by the number of
points sampled, meaning a direct dependence on run time. In practice, the
calculation can be run to an arbitrary set target tolerance.

Generally, differences in volume has a much bigger influence on the neutronics
of a reactor than do small boundary changes (with constant volume). Hence, this
is a useful measure for performance. Notably the errors found (fig. \ref{fig:voldiff}) are dominated by
the error in the stochastic volume estimator, not the volume error itself. This
is evidenced by the very small error in burst- and control-rod volume for the
detailed model, which were run with smaller tolerances.

The 5 cases considered, each have different settings for the control- and
burst-rods (see tables \ref{tab_bm_giv_rod_pos} and \ref{tab_det_giv_rod_pos}).

|case | CR 1 top | CR 2 top |BR top | $k_{eff}$ CAD0| $k_{eff}$ CAD | $k_{eff}$ CSG | $k_{eff}$ Lit.|
|-----|------|-----|-----|----|----|----|
| 1   | -4.001   | -0.449   | 0.0   | 0.98326       | 0.98026       | 0.98187       | 0.9865        |
| 2   | -1.998   | -1.666   | 0.0   | 0.98427       | 0.98101       | 0.98185       | 0.9867        |
| 3   | -0.493   | -3.794   | 0.0   | 0.98460       | 0.98124       | 0.98297       | 0.9878        |
| 4   | -0.469   | -0.447   | -2.970| 0.98446       | 0.98745       | 0.98359       | 0.9883        |
| 5   | -0.319   | -0.656   | 0.0   | 0.98980       | 0.98706       | 0.98844       | 0.9933        |
Table: Control rod (CR) and burst rod (BR) positions for the 5 cases of the Godiva-IV benchmark/simplified model from HEU-MET-FAST-086
[@icsbep_2020; @hagopian2018updating]\label{tab_bm_giv_rod_pos}. Measures in inches
withdrawn from fully inserted position. The two rightmost columns contain
criticality numbers for the device. MC refers to simulated Monte Carlo estimates, whereas Lit. refers
to numbers drawn from the benchmark report. CAD0 is a single-mesh model whereas CAD refers to a model where core, control rods, and burst rods have been discretized separately, then put together as separate universes in OpenMC.

|case | CR 1 top | CR 2 top |BR top | $k_{eff}$ CAD | $k_{eff}$ CSG | $k_{eff}$ Lit.|
|-----|------|-----|-----|----|----|----|
| 1   | -4.001   | -0.449   | 0.0   | 0.97905      | 0.98303       | 0.9880        |
| 2   | -1.998   | -1.666   | 0.0   | 0.98390      | 0.98275       | 0.9880        |
| 3   | -0.493   | -3.794   | 0.0   | 0.97885      | 0.98330       | 0.9887        |
| 4   | -0.469   | -0.447   | -2.970| 0.98352      | 0.98426       | 0.9897        |
| 5   | -0.319   | -0.656   | 0.0   | 0.98390      | 0.98969       | 0.9945        |
Table: Control rod (CR) and burst rod (BR) positions for the 5 cases of the detailed Godiva-IV model from HEU-MET-FAST-086
[@icsbep_2020; @hagopian2018updating]\label{tab_det_giv_rod_pos}. Measures in inches
withdrawn from fully inserted position. The two rightmost columns contain
criticality numbers for the device. MC refers to simulated Monte Carlo estimates, whereas Lit. refers
to numbers drawn from the benchmark report. In the detailed case, for simplicity, only a separately discretized model was run.

Each experimental case was modelled in CAD in two ways. First, the entire
reactor, including burst- and control-rods in respective positions was drawn in
CAD, exported and converted to a transport-compatible geometry. Second, the
core geometry, burst, and control-rod was exported and converted individually.
In the latter case the full reactor model is assembled when the transport
geometry is described within the confines of an OpenMC-geometry. This yields
the important advantage of a model in which the rods can be moved dynamically,
i.e. not restricted to predefined cases locked in by the discretized geometry.
In the case of the detailed model minor adjustments had to be made to the stated
measurements, in order for the model to fit.: The locking bolts had to be shortened slightly and the height of the shelf
of the inner intermediate sub-assembly plate reduced. In both cases the edits were ~= 1 mm and will not affect the results.
We assume the errors are mere misprints in the drawings.

It is clear from tables \ref{tab_bm_giv_rod_pos} and \ref{tab_det_giv_rod_pos} that the agreement between model and
experiment is not perfect, yet the difference between CAD-based and CSG-based
models is small (generally on the order of 5e-3). We take this as proof that the geometry generation
works as intended.


## Molten Salt Reactors
![Rendering of a triangularized models of the ARE- and MSRE-reactors as generated using CAD_to_OpenMC. The MSRE model (left) includes both reactor core, liquid fuel contained therein, graphite moderator stringers, as well as thermal shielding and reactor enclosure. The reactor pit _is_ included in the model, but we have excluded it from the image to make the core more visible. The ARE model (right) includes the set of three safety/shim rods and the regulating rod in the center.\label{fig:msreAre}](figs/msre_are.png)


\autoref{fig:msreAre} (left) shows pictures of the meshed MSRE and ARE geometries, including reactor enclosure control rods etc.
The Aircraft Reactor Experiment (ARE) and Molten Salt Reactor Experiment (MSRE) were two molten salt reactor experiments carried out at Oak Ridge National laboratory; ARE in November '59 and MSRE was run between '65 and 70.
In the confines of this article, this pair serves as examples of complex reactor geometries that can be handled by CAD_to_OpenMC. Very detailed CAD-models are available for the set which we have used as inputs.
\autoref{tab_bmark} tabulates the $k_{eff}$-values computed using a materials composition set to mimic the reported values as closely as possible.
It is clear that, similar to the case for the GODIVA IV-device the modelled
values are not in complete agreement with the reported ones, yet this may be
explained by possible discrepancies between the drawings in accessible reports and the actual experiment. Any engineering realities not written down in these reports, are likely lost.
This is particularly true for ARE where the discrepancy is largest, and details scarce.
It should be noticed though that the calculated keff for the MSRE case has better agreement  than any other MSRE criticality benchmarks known to the authors.

| model\label{tab_bmark} | $k_{eff}$ | $k_{eff,lit}$ |
|-------|-----------|-------|
| ARE  :| 1.07      | 1.0   |
|MSRE   | 1.00549  | 1.0|

Table: Criticality numbers, $k_{eff}$, for a set of three molten salt reactors
as computed by the combination of CAD-model, CAD_to_OpenMC, and OpenMC, using
ENDF v8.0 cross-sections. Also tabled are corresponding values found in
literature,$k_{eff,lit}$ [@cottrell_are_operation_1955; @robertson_msre_1965].

<!-- meshing timings?-->


# Discussion and Conclusion
We submit that have shown that the tool presented in this paper is a convenient tool for making CAD geometries available for particle transport. By utilizing the DAGMC-layer the resulting geometries are not restricted to OpenMC, but in fact may be used also in MCNP, fluka, etc. Experience has shown that a particularly useful feature is to extract tags from CAD-defined parts and interpret them as material tags for transport. This enables a consistent material naming scheme throughout the entire modeling procedure.



# References
