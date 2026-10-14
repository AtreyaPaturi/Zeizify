# Methods

This document describes what Zeizify computes, in enough detail to report it in a methods section and to
reproduce it. Version-specific parameters for any given run are recorded in that run's `run_manifest.json`.

## Input

A Zeiss CZI file from a multi-well-plate acquisition in ZEN, in which every scene (`S` index) is one imaged field
of one well, plus a plate map CSV assigning cell type, treatment and antibody to each well.

## 1. Scene-to-well assignment

**Metadata (primary).** For plate acquisitions, ZEN writes a `<Scene Index=… Name=…>` element per scene under
`Metadata/Information/Image/Dimensions/S/Scenes`, containing `<ArrayName>` (the well ID, e.g. `B7`) and
`<Shape><RowIndex/><ColumnIndex/></Shape>`. When every scene carries an `ArrayName`, it is used as the
assignment. Scenes are keyed by their `Index` attribute, which is the value referenced by the image subblocks,
not by their order in the XML. As an independent consistency check, `ArrayName` is compared against
`RowIndex`/`ColumnIndex` under 0- and 1-based conventions, and any disagreement is reported.

**Geometry (fallback).** If any scene lacks an `ArrayName`, assignments are inferred from each scene's
stage position, taken as the (Y, X) start coordinate of its first subblock in the subblock directory (no
pixel data is decoded).
1. The scene positions are agglomeratively clustered (average linkage, Euclidean) into *n*<sub>wells</sub>
   clusters. Fields of one well lie much closer together than adjacent wells, and wells may have unequal field
   counts.
2. Well centroids are clustered on Y (1-D average linkage) into *n*<sub>rows</sub> bands. Every band must
   contain exactly *n*<sub>cols</sub> wells.
3. Within each row band, centroids are clustered on X into *n*<sub>cols</sub> bands. Every band must contain
   exactly one well.
4. Bands are ordered by coordinate and mapped to row letters and column numbers. Because geometry alone
   cannot tell which physical corner is A1, the row and column order can be reversed (`flip_rows`,
   `flip_cols`), and the operator confirms the orientation on a QC plot.

If any of steps 1–3 does not produce a complete grid, the mapping is refused (no partial mapping is ever used)
and a diagnostic plot colored by row band is written.

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

For each antibody, one field is designated as the reference. By default this is the first field whose
treatment contains a control keyword (mock, control, ctrl, DMSO, untreated, vehicle, PBS), otherwise the first
field; the operator can choose a different one. The reference field's per-channel minimum and maximum define the
display range for every field of that antibody. Optionally, a custom (min, max) can be entered per
antibody and channel, which overrides the reference for that channel only. Ranges are recorded in
`intensity_report.csv` and are not applied to the TIFFs.

## 4. Layout images

For visualization only, each field's channels are linearly rescaled to 8 bit using the antibody's display
range (values outside it are clipped), pseudocolored by fluorophore and combined additively. Panels are each
non-nuclear channel alone, the nuclear stain alone, and each non-nuclear channel merged with the nuclear
stain, separated by a 24-pixel white gap.

## Software

Python ≥ 3.12, czifile, tifffile, imagecodecs, NumPy, SciPy (hierarchical clustering), Matplotlib (QC
plots), Pillow (layouts) and customtkinter (GUI). Exact versions for each run are listed in
`run_manifest.json → software`.
