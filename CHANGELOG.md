# Changelog

All notable changes are listed here. Versions follow [Semantic Versioning](https://semver.org/):
MAJOR for incompatible output or CLI changes, MINOR for new features, PATCH for fixes.

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
