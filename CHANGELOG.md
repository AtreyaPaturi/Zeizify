# Changelog

All notable changes are listed here. Versions follow [Semantic Versioning](https://semver.org/):
MAJOR for incompatible output or CLI changes, MINOR for new features, PATCH for fixes.

## [1.10.0] - 2026-10-22

### Added
- Figures page: *Choose on image…* for every condition. It opens the whole field image, with nuclei
  outlined (analysed / QC-flagged / chosen) and the ranked suggestions labelled. You can:
  - step through the condition's fields and switch between merge and single-channel views;
  - click a nucleus to see its values against the condition median, the crop box and any warnings
    (condensed DAPI, crop off the image);
  - use it in the figure, where it is recorded as a manual pick.
- `cellfigure.field_overview`, `nucleus_at`, `cell_summary`, `condition_fields`, `preview_ranges`.

## [1.9.0] - 2026-10-22

### Added
- **Single-cell figures**: a new *Figures* page, `zeizify figure`, and the `zeizify.cellfigure` module.
  - One representative cell per cell type × treatment, picked automatically or chosen by hand.
    Automatic picks are ranked against the condition's medians (spot counts, double-positive spots,
    Pearson's r, spot size), among interphase, isolated nuclei whose signal is nuclear. Each condition has a
    ranked list with thumbnails for manual choice.
  - Custom crop size in µm.
  - Layout: the Contrast page's panel and color choices, panel order, treatments as rows or columns, cell
    types side by side or stacked, renamed labels, and journal widths (89 / 120 / 183 mm), a 16:9 slide or
    a custom width. Font, gaps, scale bar and nucleus outline are adjustable.
  - Display: linear, background-subtracted, one range per channel shared by every condition (or per
    panel), with an optional fixed maximum.
  - Export: RGB PNG / TIFF with dpi, and PDF / SVG with editable text, plus `cells.csv`, a draft
    `legend.txt` and `figure_settings.json`.
- After an analysis, a *Single-cell figure* button on the Analyze page.

## [1.8.0] - 2026-10-19

### Added
- `superplots/<panel>_spot_categories.png`, per channel pair:
  - a 100 % stacked bar per condition: the share of spots that are double-positive or single-positive;
  - bar charts of mean double-positive, A-only and B-only spots per nucleus, each on its own axis, with
    error bars: SD of replicate means with >= 2 replicates per condition, otherwise SEM over nuclei
    (labelled).
- `spot_categories.csv` with the plotted values.

## [1.7.1] - 2026-10-19

### Changed
- Single-positive spots are also in the overview figure and as per-metric PNGs (`<pair>_n_a_only.png`,
  `<pair>_n_b_only.png`).
- Channel-pair graph titles name the stains instead of A / B (e.g. "PML-only spots per nucleus").
- `superplots/` is cleared at the start of each run, so figures from earlier versions or settings no
  longer linger.

## [1.7.0] - 2026-10-19

### Added
- Single-positive spots per nucleus for every channel pair: A spots touching no B spot (`n_a_only`) and B
  spots touching no A spot (`n_b_only`) in `colocalization.csv`, summarised and tested like the other
  metrics.
- `superplots/<panel>_single_double_positive.png`: A-only, B-only and double-positive spots per nucleus, each
  as a SuperPlot and a count distribution.

### Changed
- `spot_count_distribution.csv` has a `metric` column (`spots_n`, `n_double`, `n_a_only`, `n_b_only`).

## [1.6.2] - 2026-10-16

### Added
- `superplots/<panel>_spot_area_per_nucleus.png`: average spot area per nucleus (size at half maximum), one
  SuperPlot per channel, all cell types side by side. Nuclei without spots are not included.

## [1.6.1] - 2026-10-16

### Added
- Brightness-independent spot size: area and diameter at half of each spot's own peak (`area_halfmax_um2`,
  `diam_halfmax_um` in `spots.csv`; `<ch>_spot_mean_area_halfmax_um2` per nucleus), summarised, tested and
  plotted in the overview.

### Changed
- The thresholded spot area is labelled "thresholded"; it grows with spot brightness (RPE1 PML: Spearman
  0.63 with brightness). Untreated KO1 PML bodies, about 40 % brighter than WT, measured 12 % larger by
  thresholded area but identical at half maximum. Compare sizes with the half-maximum measure.

## [1.6.0] - 2026-10-15

### Added
- Fixed (experiment-wide) spot threshold per channel: one brightness floor for every well, automatic
  (Otsu on log spot brightness pooled over all wells) or a typed value. *Threshold* per channel on the
  Analyze page; `--fixed-threshold CHANNEL[=VALUE]`; `spot_fixed_threshold` parameter. The floors used and
  their separability are recorded in `analysis_manifest.json`.
- `qc/spot_brightness.png`: spot brightness per treatment and cell type for every channel, with the floor
  applied or suggested, and a warning in the log when a fixed floor is used on a channel without two
  separated populations.

### Fixed
- The per-nucleus threshold alone called grainy diffuse signal "spots" when a treatment dissolves the bodies
  (RPE1 SUMO2/3 after ML792: about 9-13 false spots per nucleus, about 5x dimmer than real SUMO2/3 bodies).
  Use the fixed threshold for such stains. Defaults are unchanged.

## [1.5.0] - 2026-10-15

### Added
- Double-positive spots per nucleus for every channel pair (`n_double` in `colocalization.csv`): connected
  spot regions containing spots of both channels, so a body overlapping two foci of the other channel
  counts once. Also `n_a_overlap` / `n_b_overlap`, the number of each channel's spots overlapping the other.
- `n_double` is summarised and tested like the other metrics, plotted in the overview, and added as a row
  (SuperPlot + count distribution) to `<panel>_spots_per_nucleus.png`.

## [1.4.1] - 2026-10-15

### Changed
- SuperPlot style: nuclei as a beeswarm colored by replicate (Okabe-Ito colors), replicate medians as large
  dots of the same color with area proportional to the number of nuclei, a boxed frame, a two-level x axis
  (cell type under each column, treatment underlined beneath), and a legend with replicate colors and an
  *n* size key. Without a Replicate column, replicates (wells) are numbered 1, 2, 3 within each condition.

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
