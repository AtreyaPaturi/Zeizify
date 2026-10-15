# Methods

This document describes what Zeizify computes, in enough detail to report it in a methods section and to
reproduce it. Version-specific parameters for any given run are recorded in that run's `run_manifest.json`.

## Input

A Zeiss CZI file from a multi-well-plate acquisition in ZEN, in which every scene (`S` index) is one imaged field
of one well, plus a plate map CSV assigning each well a cell type, a treatment, and the antibody detected in
each fluorophore channel (`Antibody_<fluorophore>` columns). A well may carry several antibodies in different
channels (multiplexed or co-localization staining). Its *staining panel* is the ordered set of its
per-channel antibodies (e.g. `CD7+HLA`).

## 1. Scene-to-well assignment

**Metadata (primary).** ZEN writes a `<Scene Index=… Name=…>` element per scene under
`Metadata/Information/Image/Dimensions/S/Scenes`, which records the scene's well as `<ArrayName>` (plate-carrier
acquisitions) or as the name of its `<Shape Name="B7">`, together with `<Shape><RowIndex/><ColumnIndex/></Shape>`.
When every scene carries a well ID (`ArrayName`, else a Shape name of the form letter + number), it is used
as the assignment. Scenes are keyed by their `Index` attribute, which is the value referenced by the image subblocks,
not by their order in the XML. As an independent consistency check, `ArrayName` is compared against
`RowIndex`/`ColumnIndex` under 0- and 1-based conventions, and any disagreement is reported.

**Geometry (fallback).** If any scene lacks a well ID, assignments are inferred from each scene's stage
position, taken as the (Y, X) start coordinate of its first subblock in the subblock directory (no pixel data
is decoded). Only part of the plate may have been imaged, so the number of wells is not assumed.
1. Y and X coordinates are each sorted and split into bands (rows and columns) at gaps wider than 0.35 × the
   nominal well pitch (SBS pitch for the plate format, e.g. 9 mm for 96 wells, divided by the pixel size).
   Without a pixel size, the split is at the gaps that clearly exceed the rest (> half the largest gap,
   provided it is > 4× the median gap).
2. The working pitch is the median spacing between neighbouring band centres (mid-ranges) that are one well
   apart, which absorbs small differences between stage scaling and the nominal pitch.
3. Every band must be narrower than 0.85 pitch, and every band centre must lie within 0.35 pitch of the
   integer grid. Otherwise the mapping is refused, so neighbouring rows or columns (e.g. under plate tilt) are
   never merged or misnumbered.
4. Bands are numbered by their distance from the first band in pitches, so skipped rows and columns are
   preserved, and labelled relative to the operator-supplied first imaged well (the imaged well closest to A1;
   default A1). Because geometry alone cannot tell which physical corner is A1, the row and column order can be
   reversed (`flip_rows`, `flip_cols`).

The operator confirms the first well and orientation on a QC plot. If any check fails, no mapping is used, and
a diagnostic plot colored by row band is written.

In both modes, a QC table (`scene_to_well.csv`) and plot (`scene_to_well_qc.png`) are written. Extraction
does not start until the mapping has been confirmed: interactively in the GUI, by a prompt in the CLI, or
explicitly with `--yes`.

## 2. Projection and export

For each labeled well and each of its scenes (fields), only that scene's level-0 subblocks are read
(`czifile`, pyramid levels excluded). Subblocks are grouped by channel (`C`) and reduced by a pixel-wise
maximum over Z (maximum-intensity projection). Scenes containing more than one tile position (`M`/`X`/`Y`) or
more than one time point (`T`) are rejected, so that tiles or time points are never merged. Each channel's
projection is written once as a single-plane TIFF (`tifffile`), in the original bit depth and data type, with
the X pixel size from `Metadata/Scaling/Items/Distance[@Id='X']` (meters per pixel, converted to µm) embedded as
ImageJ resolution metadata.

Channel names are taken from `Metadata/Information/Image/Dimensions/Channels/Channel/@Name` in `C`-index
order. Field names are the scene `Name` attributes (e.g. `P1`).

A field/channel is flagged as saturated when its maximum equals or exceeds 2<sup>*b*</sup> − 1, where *b* is
`Information/Image/ComponentBitCount` (or the unique per-channel value).

## 3. Reference contrast

Each extracted channel is tagged with the antibody the plate map assigns to that channel in that well. CZI
channel names are matched to plate-map fluorophore names, case- and punctuation-insensitively, falling back
to a unique match on the three-digit wavelength. Fields are grouped by staining panel, so only identically
stained wells share a reference. For each panel, one field is designated as the reference. By default this is
the first field whose treatment contains a control keyword (mock, control, ctrl, DMSO, untreated, vehicle,
PBS), otherwise the first field; the operator can choose a different one. The reference field's per-channel
minimum and maximum define the display range for every field of that panel. Optionally, a custom (min, max)
can be entered per panel and channel, which overrides the reference for that channel only. Ranges are
recorded in `intensity_report.csv` and are not applied to the TIFFs.

## 4. Layout images

For visualization only, each field's channels are linearly rescaled to 8 bit using the panel's display range
(values outside it are clipped), pseudocolored and combined additively, so overlapping signal takes the
mixed color (e.g. green + magenta → white). The operator selects the panels and each channel's color (red,
green, blue, magenta, cyan, yellow or gray). Defaults are inferred from fluorophore names. The default panels
are each signal channel alone, the nuclear stain alone, each signal channel merged with the nuclear stain,
and, when two or more signal channels exist, all signal channels merged with and without the nuclear stain.
Panels are placed side by side with a 24-pixel white gap and captioned with antibody names. The chosen panels
and colors are recorded in `run_manifest.json`. Zeizify does not compute co-localization coefficients.

## Software

Python ≥ 3.12, czifile, tifffile, imagecodecs, NumPy, SciPy (hierarchical clustering), Matplotlib (QC
plots), Pillow (layouts) and customtkinter (GUI). Exact versions for each run are listed in
`run_manifest.json → software`.
