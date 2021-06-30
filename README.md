# 2D Simulator

## Quick start

To run the simulator, use the ``makeSim`` command, for example:

    makeSim --visit 35 --detector b1 --pfsDesignId 3 --exptime 30 --type flat --imagetyp flat_odd

The above example generates a 30-second simulated flat exposure corresponding to visit 35,
for detector ``b1``, and for odd-numbered fibers, with the active fibers provided
in pfsDesign ID 3.

    makeSim --visit 58 --detector m1 --pfsConfig --pfsDesignId 1 --exptime 300 --type object --objSpectraDir /projects/HSC/PFS/simulator/

The next example above a 300-second simulated object exposure corresponding to visit 58,
for detector ``m1``, with the objects associated with each active fiber listed
in pfsDesign ID 1. Note in the case of object exposures, object spectra are required.
Their location are provided by the ``objSpectraDir`` argument. In this case, `/projects/HSC/PFS/simulator`.

See the simulator [User Guide](https://sumire.pbworks.com/w/file/fetch/131637714/PFS-DRP-PRU060011-01_PFS2DDRPSimulatorUserDocumentation.pdf) for more information regarding usage.

## Object Spectra catalogs

For running science object simulations (ie where the ``--type`` is ``object``), object spectra are required,
in the form of `pfsSimObject` FITS files. Each object spectrum (``pfsSimObject``) file is identified by the
catalog ID or ``catId``, and the object identifier for the object within that catalog ``objId``.
The object spectra files needed by the simulator must be organized under a single root directory, which
is the value of the ``--objSpectraDir`` argument.

This root directory must contain a file ``catalog_config.yaml``, which specifies the available catalogs,
their unique ``catId`` s, and their location relative to the root directory.
Specifically, the YAML file must contain a sequence of one or more mapping nodes,
where each mapping node consists of three key-value pairs: the `catId`, the catalog `name`, and
the `rel_location` (relative location).

Below is an example of a
``catalog_config.yaml`` file that is used by the Princeton group:

```
- catId: 1
  description: GE z<2 galaxy spectra built on the COSMOS catalog. Galaxies only
  rel_location: pfsSimObjects/lowz_COSMOS_2020_12_14
- catId: 5
  description: 10,000 GA stars from Laszlo from 2021-04
  rel_location: pfsSimObjects/GA_Laszlo_2021_04_01
- catId: 6
  description: GE combined z<2, tomography, LAEs, and galaxies
  rel_location: pfsSimObjects/GE_combined_2021_08_05
```

And an example of the directory structure below the ``--objSpectraDir`` root directory is as follows:
```
├── catalog_config.yaml
├── pfsDesign
│   ├── pfsDesign-0x0000000000000009.fits
│   ├── pfsDesign-0x000000000000000a.fits
│   ├── pfsDesign-0x000000000000000b.fits
│   ├── pfsDesign-0x000000000000000c.fits
│   ├── pfsDesign-0x000000000000000d.fits
│   └── pfsDesign-0x000000000000000e.fits
└── pfsSimObjects
    ├── lowz_COSMOS_2020_12_14
    |   ├── pfsSimObject-00001-00000-0,0-000000000002098f.fits
    |   ├── pfsSimObject-00001-00000-0,0-000000000004011f.fits
    |   └── [...]
    ├── GA_Laszlo_2021_04_01
    |   ├── pfsSimObject-00005-00000-0,0-000000000000270e.fits
    |   ├── pfsSimObject-00005-00000-0,0-00000000000025be.fits
    |   └── [...]
    └── GE_combined_2021_08_05
        ├── pfsSimObject-00006-00000-0,0-0000000000012345.fits
        └── [...]
```
