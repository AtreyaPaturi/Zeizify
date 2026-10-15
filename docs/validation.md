# Validation

## Cross-validation against existing tools

No existing open-source tool combines plate-map annotation, scene-to-well mapping, per-channel projection
and control-referenced contrast for Zeiss plate CZIs, so the individual steps were checked against the
tools that implement each part.

| Step | Reference implementation | Zeizify | Agreement |
|---|---|---|---|
| Well ID per scene | [czitools](https://github.com/sebi06/czitools) by ZEISS, `CziSampleInfo`: reads `ArrayName` from `Information.Image.Dimensions.S.Scenes.Scene`, falling back to `Shape.Name` | Same element path, same fields and order | ✔ (the `Shape.Name` fallback was added in 1.2.1 after a real acquisition without `ArrayName`) |
| Well row/column | czitools: `Shape.RowIndex` / `Shape.ColumnIndex` (`well_rowID`, `well_colID`) | Used as an independent cross-check of `ArrayName` | ✔ (adopted from czitools) |
| Scene identity | czitools maps scenes by list position | Maps by the scene's `Index` attribute, which is what subblocks reference | Equivalent for ZEN-written files; more robust otherwise |
| Missing well IDs | czitools: scenes with no `ArrayName` are left unmapped | No partial mapping; falls back to geometry for the whole file, with human confirmation | Stricter |
| Channel names | czitools `CziChannelInfo`: `Information.Image.Dimensions.Channels.Channel/@Name`, list order = C index | Same path and ordering | ✔ (adopted; the original script took the first `<Channel Name>` anywhere in the XML, which can follow acquisition-track order instead) |
| Pixel size | czitools `CziScaling`: `Scaling.Items.Distance.Value` × 10<sup>6</sup> → µm | `Distance[@Id='X']/Value` × 10<sup>6</sup> | ✔ (selected by `Id` rather than list position) |
| Subblock selection | [czifile](https://github.com/cgohlke/czifile) `filtered_subblock_directory`: pyramid level 0 only | Uses it directly | ✔ |
| Projection | Standard maximum-intensity projection over Z | `np.max` over Z per channel | ✔ |

## Changes made during validation

The original scripts were reviewed against the references above. The following were changed:

1. **Channel names** are read from the `Information` block rather than from the first `<Channel Name=…>` in the
   document, which could follow acquisition-track order (covered by `tests/test_metadata.py`).
2. **Scenes are parsed as XML** rather than by regular expression, so attribute order (`Name` before `Index`)
   no longer matters.
3. **ArrayName ↔ RowIndex/ColumnIndex cross-check** was added, following czitools.
4. **Mosaic and time-lapse guard.** Previously, every subblock of a scene was stacked and maximum-projected
   together, so a scene with several tiles or time points would have been silently merged. Such scenes are
   now rejected.
5. **Projection output is squeezed to 2-D** before writing, so TIFFs are single-plane regardless of the
   subblock's singleton dimensions.
6. **Plate-map labels are re-slugged on load**, so a hand-edited CSV cannot put path separators into output
   file names.
7. **CZI file handles are closed** after use (context manager).
8. **Configurable plate format** (6–384 wells) and **configurable channels**. The hard-coded
   AF488/AF594/DAPI assumptions are now defaults inferred from fluorophore names.

## Automated tests

`pytest` runs 50 tests on synthetic data:

- metadata parsing, including a decoy channel order and reversed attribute order;
- geometric mapping on a simulated 96-well plate with stage tilt (0.4 %) and positional jitter (σ = 40 px),
  in all four orientations, with unequal field counts, and on a 384-well plate;
- partially imaged plates: a block starting at an arbitrary first well, a skipped column, and a simulation of
  the columns-2–5 acquisition above (wide field spread, stage pitch 4 % below nominal);
- well IDs from `Shape/@Name`, and rejection of non-well shape names (e.g. tile region `TR1`);
- refusal on ambiguous geometry (extreme tilt), with a diagnostic plot written;
- metadata precedence, fallback when a single `ArrayName` is missing, and an out-of-plate well ID error;
- Z projection values, mosaic refusal, ImageJ pixel-size metadata, saturation flags and resume without rewriting;
- control selection, custom ranges, intensity report contents and layout PNG geometry;
- multi-stain (co-localization) plate maps: per-channel antibody columns, panel grouping, stain-to-channel
  matching, co-localization merge panels, user colors and migration of older saved maps;
- an end-to-end `zeizify extract` run in both mapping modes (metadata and geometry) that verifies every
  scene's assigned well against ground truth.

## Validation on real data

The metadata and geometric methods were run independently on two acquisitions:

| Acquisition | Well IDs stored as | Imaged | Metadata vs geometry |
|---|---|---|---|
| 96-well, 288 scenes, 3 fields/well | `ArrayName` | full plate (A1–H12) | 288 / 288 scenes agree |
| 96-well co-localization, 160 scenes, 5 fields/well, 0.086 µm/px, 133 GB | `Shape/@Name` (no `ArrayName`) | partial: columns 2–5 (A2–H5) | 160 / 160 scenes agree (geometry with first well A2); well IDs also agree with 1-based Shape Row/ColumnIndex for all scenes |

The second file exposed two defects in version 1.2.0, both fixed in 1.2.1. First, its well IDs (in
`Shape/@Name`) were ignored. Second, the geometric fallback assumed a fully imaged plate, forced 32 wells
into 96 clusters, and failed. The measured well pitch on this stage was 103,829 px, 0.5 % from the nominal
9 mm / 0.08625 µm = 104,348 px. To repeat this on
your own file:

```bash
zeizify map plate.czi -o check/            # uses ArrayName if present
```

To force the geometric path for comparison, use the Python API:

```python
from czifile import CziFile
from zeizify.mapping import scene_positions, geometric_mapping, mapping_from_metadata
from zeizify.czi_metadata import read_scenes

with CziFile("plate.czi") as czi:
    pos = scene_positions(czi.filtered_subblock_directory)
    meta = mapping_from_metadata(read_scenes(czi.metadata()))
geo, _ = geometric_mapping(pos)          # add flip_rows/flip_cols if needed
print(sum(geo[s] == meta[s] for s in pos), "/", len(pos), "scenes agree")
```
