# Changelog

All notable changes are listed here. Versions follow [Semantic Versioning](https://semver.org/):
MAJOR for incompatible output or CLI changes, MINOR for new features, PATCH for fixes.

## [1.4.0] - 2026-10-15

### Added
- All cell types (e.g. WT and KO) share each SuperPlot, side by side in their own color, one overview per
  staining panel (`superplots/<panel>_overview.png`, per-metric PNGs in `superplots/<panel>/`).
- Cell types are compared directly within each treatment: `stats_pairwise.csv` has new `compare`
  (`treatment` / `cell_type`) and `treatment` columns.
- `superplots/<panel>_spots_per_nucleus.png`: spots per nucleus for each channel, as a SuperPlot and as the
  fraction of nuclei per count class (0, 1-5, 6-10, ...; classes scale to the data). The fractions are in
  `spot_count_distribution.csv`.
- Treatment and cell-type order: control / wild type first, then plate-map order (e.g. MOCK, ML792,
  MG132, ML792+MG132). Override with *Treatment order* / *Cell-type order* on the Analyze page, or
  `--treatment-order` / `--cell-type-order`.

### Changed
- Treatments were previously ordered control first, then alphabetically; per-cell-type SuperPlot files
  (`<panel>_<cell type>_overview.png`) are replaced by one file per panel.

## [1.3.0] - 2026-10-15

### Added
- **Zeizify Analyze**, a per-nucleus quantification of nuclear bodies and colocalization that runs on the
  extracted TIFFs (new *Analyze* page, `zeizify analyze`, `zeizify.analysis` package; methods in
  `docs/analysis.md`).
  - Nuclei: Cellpose (GPU) or classical Otsu + watershed segmentation, size window, border exclusion, and
    QC flags for condensed (mitotic/apoptotic), dim and saturated nuclei. Masks are cached.
  - Spots: per-nucleus white top-hat, a per-object threshold with a robust noise floor, declumping, and a µm
    size window.
  - Per-nucleus counts, sizes, shapes, background-paired intensities, and threshold-free dispersal (CV,
    enrichment, fraction in spots).
  - Colocalization: Pearson, Spearman, Li ICQ, Manders M1/M2 and Costes tM1/tM2 (Coloc 2 definitions), an
    optional Costes significance test, and object overlap / enrichment / nearest-neighbour distances.
    Single-stain controls are labelled.
  - Statistics on replicate medians: Welch ANOVA, Games-Howell, BH q-values, Hedges' g, Cliff's delta, and
    SuperPlots.
  - Parallel measurement over CPU cores; a one-field Preview for tuning; a progress panel with Stop.
- After *Apply & export*, an *Analyze* shortcut on the Run page.
- `scikit-image` and `pandas` dependencies; optional extra `zeizify[cellpose]`.

### Fixed
- `python -m zeizify` is guarded for multiprocessing on Windows.

## [1.2.2] - 2026-10-15

### Added
- Export progress on the Contrast page. *Apply & export* opens a panel with the current step (intensity
  report → layout PNGs → run manifest), a progress bar with percentage, the field being written, counts
  (written / kept / skipped), elapsed time and an estimate of the time left, a live log, and a *Stop*
  button. Failures are shown in the panel with the error.
- `generate_png_layouts` / `pipeline.finalize` accept `progress`, `log` and `cancel` callbacks.

### Fixed
- The export creates the output folder if it is missing.

## [1.2.1] - 2026-10-15

### Fixed
- **Well IDs stored as Shape names are now read.** ZEN acquisitions that do not write `<ArrayName>` record each
  scene's well as `<Shape Name="A2">`. Zeizify now reads it, as czitools does, instead of falling back to
  geometry.
- **Partially imaged plates no longer fail in the geometry fallback.** The fallback assumed every well of the
  plate was imaged and forced the scenes into 96 clusters, so imaging e.g. only columns 2–5 failed with
  "row-band sizes … expected 12". It now detects the imaged rows and columns from the gaps between wells,
  accounts for skipped rows/columns using the well pitch, and refuses (rather than guesses) when a band is
  wider than a well or off the grid.

### Added
- "First imaged well" setting for the geometry fallback (GUI field on the mapping card, CLI `--first-well`,
  defaulting to the labeled well closest to A1 in the plate map), recorded in `run_manifest.json`.

## [1.2.0] - 2026-10-14

### Added
- **Co-localization / multi-stain plate maps.** Declare any number of fluorophore channels in the plate map
  editor (quick-add for common dyes). Each gets its own antibody category, so a well can carry a different
  antibody in each channel. Wells show one wedge per category.
- Plate map CSV gains one `Antibody_<fluorophore>` column per channel. `Antibody` becomes the staining panel
  (e.g. `CD7+HLA`), filled in automatically if missing, as is `Label`.
- Every extracted channel is tagged with its antibody (`Stain` column in `intensity_report.csv`), and wells are
  grouped by panel for output folders and reference contrast.
- **Figure layout options** on the Contrast page: choose any single-channel, nuclear-merge or co-localization
  merge panels and each channel's pseudocolor (green/magenta preset). Panels are captioned with antibody names.
- CLI: `--panel CH+CH` and `--color CH=COLOR`.
- `examples/plate_map_colocalization.csv`.

### Changed
- Default layouts add `AF488+AF594` and `AF488+AF594+DAPI` merges when two or more signal channels exist.
- Saved single-antibody plate maps are migrated, with the antibody assigned to the AF594 channel (rename it in
  the editor if needed).

## [1.1.0] - 2026-10-14

### Added
- Version shown in the window title, as a badge next to the logo, and in the sidebar footer.
- About dialog (click the version badge or *About*) listing the Zeizify, Python and library versions, with
  *Copy version info* for bug reports and methods sections.
- Help panel for every page (Plate map, Run, Contrast): click *Help* in the page header or press F1.

## [1.0.0] - 2026-10-14

### Added
- First public release: `src/` package layout, desktop GUI (plate map editor, run timeline, contrast),
  CLI (`zeizify info | map | extract`), run manifest, test suite, methods and validation docs.
