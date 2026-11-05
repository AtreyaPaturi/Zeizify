# Zeizify Analyze - methods

Zeizify Analyze quantifies nuclear bodies (nuclear speckles, PML bodies, foci) and the colocalization of two
or more proteins, nucleus by nucleus, from the per-field TIFFs Zeizify extracts. Each step adopts the approach
of a published tool or study (listed below). Every parameter of a run is recorded in
`analysis/analysis_manifest.json`.

## Input

A Zeizify output folder: per-field, per-channel maximum-intensity projections, `intensity_report.csv`
(well, field, treatment, cell type, staining panel and the antibody in each channel), and
`run_manifest.json` (pixel size and sensor bit depth). Analysis settings are given in µm and converted with
the pixel size.

## 1. Nuclear segmentation

**Cellpose** (Stringer et al., *Nat Methods* 18:100–106, 2021) is used when installed: the generalist `cpsam`
model in Cellpose ≥ 4, or the `nuclei` model in v2/v3. The DAPI image is downscaled so a typical nucleus is
~30 px across, segmented, and the labels are upscaled (nearest neighbour). Without Cellpose, a **classical**
segmentation follows CellProfiler's IdentifyPrimaryObjects:
1. Gaussian smoothing (σ = 1 px at ~40 px per nuclear diameter);
2. Otsu threshold and hole filling;
3. removal of objects below 5 % of the typical nuclear area;
4. declumping by watershed on the smoothed distance transform, with seeds ≥ 0.3 diameters apart;
5. boundary refinement at full resolution.

With both methods, nuclei are then filtered by a size window (0.25–4 × the typical nuclear area), and nuclei
touching the image border are removed, as in the Aktas lab speckle workflow (*Cell* 2026). Masks are cached
and reused while nuclear settings are unchanged.

**Nucleus QC.** Every segmented nucleus within the size window is analysed, whatever its DAPI intensity:
DAPI is used only to find the nuclei. `dapi_ratio`, a nucleus's mean DAPI over the median of its field, is
reported in `nuclei.csv` but not used. A nucleus is flagged (`qc_pass = False`) only when > 1 % of its pixels
are saturated in a measured signal channel, because its intensities would then be clipped. Flagged nuclei
stay in `nuclei.csv` and are excluded from summaries and statistics.

## 1b. Cells and cytoplasm (optional)

Wells stained for a cell-body or cell-edge marker get whole-cell boundaries. The marker is found in the
plate map: F-actin / phalloidin, CellMask, tubulin, WGA or a membrane / cytoplasm stain. A channel can also
be chosen explicitly. Other wells are analysed on nuclei only.

- **Cellpose-SAM** (`cpsam`, Cellpose ≥ 4; Pachitariu, Rariden & Stringer, *bioRxiv* 2025) is used when
  Cellpose is installed: the generalist model that separates crowded, touching cells. It is given the
  marker and the nuclear stain together, downscaled so a typical cell (default 45 µm) is about 90 px
  across. With Cellpose 2/3, the `cyto3` model is used.
- **Classical** (without Cellpose): a seeded watershed from the nuclei on the marker's gradient, within the
  marker's Otsu foreground and one typical cell diameter of a nucleus, like CellProfiler's
  IdentifySecondaryObjects.
- **Assignment.** Each cell is tied to its nucleus (cell id = nucleus id):
  - a nucleus belongs to the cell holding ≥ 50 % of it;
  - a cell holding several nuclei (touching cells found as one) is split between them by nearest
    nucleus;
  - a cell with no analysed nucleus is dropped.
- **Cytoplasm** = cell minus nuclei. Cells cut by the image edge get no cytoplasm measurements.
- **Background.** Where cells were segmented, background is taken > 1 µm outside the cells, so cytoplasmic
  signal is not mistaken for background. In confluent fields (< 1 % cell-free) it falls back to outside
  the nuclei.

Per cell (`nuclei.csv`): `cell_area_um2`, `cyto_area_um2`, `cell_marker`, and per channel
`<ch>_cyto_mean(_bgsub)`, `<ch>_cyto_integrated`, `<ch>_cell_mean_bgsub` and `<ch>_nc_ratio` (nuclear /
cytoplasmic mean, both background-subtracted). The cell area, cytoplasmic intensity and N/C ratio are
summarised, tested and plotted like the nuclear metrics. The QC overlays show cell boundaries in cyan.

Validation:
- synthetic fields: cell masks match the true cell bodies (IoU > 0.85), touching cells found as one are
  split between their nuclei, and N/C ratios are recovered;
- RPE1 cells stained for F-actin and cyclin A (24 fields, 1,361 nuclei): Cellpose-SAM gave 80–92 % of
  nuclei a complete cell (most of the rest were cut by the image edge), in about 2 s per field on a GPU.

## 2. Spot (speckle / body) detection

The approach is "nuclear mask, then spot objects", as in Alexander et al. (*Nat Cell Biol* 2024), the Aktas
lab (*Cell* 2026) and the per-nucleus CellProfiler workflow of the condensate-morphology study (*Cell* 2026).
For each nucleus and channel:
1. **Background removal.** A white top-hat (disk radius 0.75 × the maximum spot diameter) is applied to the
   lightly smoothed signal (σ = ¼ of the minimum spot diameter), removing the smooth nucleoplasmic component.
2. **Per-nucleus threshold** = max(Otsu of the top-hat values in that nucleus, median + *k* · 1.4826 · MAD).
   Thresholding per object adapts to cell-to-cell staining differences. The robust noise floor (*k* = 6 by
   default) stops dispersed nuclei from being split into noise spots. This matters because speckles are rarely
   truly absent (Ilik et al., *eLife* 2020), so a count-based readout degrades at the dispersed end. The
   default was chosen on RPE1 SUMO2/3–PML images: at *k* = 3 grainy nucleoplasm was called as spots (e.g.
   36 SUMO "spots" in a nucleus with ~12 bodies); at *k* = 6 only discrete bodies remained, and PML counts
   were identical for *k* = 3–8.
3. **Declumping** of touching spots by watershed on the top-hat image, seeded at local maxima (optional).
4. **Fixed (experiment-wide) floor**, optional per channel. A per-nucleus threshold always finds the brightest
   part of a nucleus, so when a treatment dissolves the bodies it calls grainy diffuse signal spots. On RPE1
   cells treated with the SUMO E1 inhibitor ML792, it found about 9-13 SUMO2/3 "spots" per nucleus. Their
   mean top-hat intensity was about 600, against about 3,000-3,800 for SUMO2/3 bodies in untreated cells,
   and their area of about 0.10 µm² was at the size limit. For such stains, spots whose mean top-hat
   intensity is below one floor are dropped in every well, as with a global threshold in CellProfiler. The
   floor is either typed in or automatic:
   - the per-nucleus candidates of all fields are pooled first;
   - Otsu's threshold on their log brightness gives the floor;
   - Otsu's separability (between-class / total variance) is reported with it. A single bell-shaped
     population scores about 0.64, so below 0.7 a warning is logged.

   `qc/spot_brightness.png` shows every channel's spot brightness by treatment, with the floor applied or
   suggested. On the RPE1 data:
   - SUMO2/3: two populations (separability 0.80, floor 1,549). With the floor, SUMO2/3 spots per nucleus
     were 0 after ML792, versus 7-9 in untreated cells.
   - PML: one population (separability 0.62). PML bodies persist under ML792, so PML keeps the per-nucleus
     threshold; an automatic floor there would have removed most untreated PML bodies.
5. **Size window** on area-equivalent diameter (default 0.25–3 µm), applied before step 4.

Each spot is assigned to its parent nucleus, the equivalent of CellProfiler's RelateObjects.

**Compartments.** Spots can be detected in the nucleus, in the cytoplasm, or both (`spot_compartments`;
*Quantify spots in* on the Analyze page). Cytoplasmic spots need cell boundaries (section 1b). They are
detected in each cell's cytoplasm with the same steps as in a nucleus: top-hat, per-cell threshold with the
noise floor, declumping, size window and any fixed brightness floor. Per cell, the outputs are:
- `<ch>_cyto_spots_n`, `<ch>_cyto_spots_area_um2`, `<ch>_cyto_spots_area_frac`,
  `<ch>_cyto_spot_mean_area_halfmax_um2` and `<ch>_cyto_spot_mean_intensity`;
- with both compartments, `<ch>_cell_spots_n`: nucleus + cytoplasm.

`spots.csv` has a `compartment` column. With the nucleus not selected, the nuclear spot columns, and the
object colocalization that needs nuclear spots, are left empty. `<panel>_spots_by_compartment.png` shows
spots per cell in the nucleus, the cytoplasm and the whole cell as SuperPlots, one row per channel.

A diffuse cytoplasmic protein can produce many faint, near-minimum-size "spots": the per-cell threshold
finds the brightest grains, just as it does in a nucleus. Check the QC overlays and
`qc/spot_brightness.png`, and use the fixed brightness floor (section 2, step 4) when the counts are not
real bodies.

## 3. Measurements

**Per spot** (`spots.csv`): area, equivalent diameter, eccentricity, solidity, mean / max / integrated
intensity, integrated top-hat intensity, centroid, and the **size at half maximum** (`area_halfmax_um2`,
`diam_halfmax_um`).

**Spot size: use the half-maximum size to compare conditions.** The thresholded area depends on
brightness:
- a spot's outline sits where its signal crosses the nucleus threshold;
- that threshold is a smaller fraction of a brighter spot's peak, so the brighter spot gets the larger
  outline.

For spots at or near the resolution limit, as most PML bodies are, the thresholded area therefore partly
measures brightness. The half-maximum size does not depend on brightness. It is the area of the connected
region around the spot's peak where the top-hat signal is ≥ half that peak, without crossing into
neighbouring spots.

On RPE1 PML:
- thresholded area correlated with spot brightness (Spearman 0.63) and rose from 0.10 to 1.0 µm² from
  the dimmest to the brightest spots;
- untreated KO1 PML bodies, about 40 % brighter than WT, had a 12 % larger thresholded area, but an
  identical half-maximum size (0.231 µm² both).

Below about twice the resolution limit, the half-maximum size approaches the size of the microscope's
point-spread function rather than the body's physical size.

**Per nucleus** (`nuclei.csv`): area, eccentricity, solidity, mean and integrated DAPI. Then, per channel:
- **Intensity:** mean, median and integrated intensity. Background (`bg`) is the 10th percentile of pixels
  > 1 µm from any nucleus, paired as in SpeckleSPEC (github.com/Chalietia/CellProfiler); `mean_bgsub` is the
  background-subtracted mean.
- **Spots:** count, area, area fraction, and mean spot size, eccentricity, solidity and intensity.
- **Threshold-free readouts:**
  - `cv`: the coefficient of variation of the background-subtracted nuclear signal (high when the signal is
    concentrated in bodies, low when diffuse);
  - `enrichment`: spot vs nucleoplasm mean intensity;
  - `frac_in_spots`: the fraction of nuclear signal inside spots.

  These are reported alongside counts. CV is not, to our knowledge, a published, validated speckle-dispersal
  metric; treat it as exploratory.

## 4. Colocalization (per nucleus, every pair of analysed channels)

**Pixel-based**, on background-subtracted intensities (clipped at 0) inside the nucleus. Definitions match
Fiji's Coloc 2 (github.com/fiji/Colocalisation_Analysis):
- **Pearson's r**, **Spearman's ρ**, and **Li's ICQ** (fraction of pixels with (a − ā)(b − b̄) > 0, minus 0.5).
- **Costes automatic thresholds** (Costes et al., *Biophys J* 2004):
  - an orthogonal (total-least-squares) regression b = m·a + c;
  - the threshold on the stepped channel (a if |m| ≤ 1, otherwise b; the other channel follows the line) is
    found by bisection where Pearson's r of the below-threshold pixels (a < t_a **or** b < t_b) reaches 0.
- **Manders M1/M2**: Σa where b > 0 / Σa, and vice versa. Inside a nucleus these are ≈ 1 and rarely
  informative.
- **Thresholded tM1/tM2**: Σa where b ≥ t_b / Σa, and vice versa, as in Coloc 2.
- **Optional Costes significance**: the fraction of block-scrambled images (blocks ≈ PSF, default 0.25 µm)
  whose r ≥ the observed r.

**Object-based**, after SpotitPy (*BMC Bioinformatics* 2022) and DiAna (Gilles et al., *Methods* 2017):
- the number of A spots overlapping any B spot (`n_a_overlap`), and vice versa, with the corresponding
  fractions;
- the number of **single-positive spots** per nucleus: A spots touching no B spot (`n_a_only`), and B spots
  touching no A spot (`n_b_only`);
- the number of **double-positive spots** per nucleus (`n_double`): connected regions of A-or-B spot pixels
  that contain both an A and a B spot (>= 1 shared pixel), so a body overlapping two foci of the other
  channel counts once;
- object Manders (overlapping area / spot area);
- overlap relative to chance: obj_M1 divided by the fraction of the nucleus covered by B spots (> 1 means more
  overlap than random placement);
- the median nearest-neighbour distance between spot centroids.

**Bleed-through.** Coloc 2's documentation warns that emission bleed-through can look like perfect
colocalization. Wells with an antibody in only one channel are labelled `single-stain control` in
`colocalization.csv`; their apparent colocalization estimates crosstalk.

## 5. Statistics

Following **SuperPlots** (Lord et al., *J Cell Biol* 219:e202001064, 2020), nuclei are not treated as
independent samples. Each metric is summarised per replicate (the median of its QC-passed nuclei) in
`replicate_summary.csv`, and tests compare replicate medians:
- A replicate is a well by default, or the value of a `Replicate` column in the plate map (biological
  replicates). Fields can be selected as technical replicates, but they are not independent.
- **Welch's ANOVA** across treatments (SciPy `f_oneway(equal_var=False)`), run within each staining panel and
  cell type.
- **Games-Howell** pairwise comparisons, as used with CellProfiler output by Kono et al. 2023
  (github.com/Kimura-Lab/Kono-et-al.-2023).
- **Benjamini-Hochberg** q-values across all omnibus tests of a run.
- **Effect sizes**: Hedges' g on replicate medians, and Cliff's delta on nuclei (descriptive).
- **Cell types** (e.g. WT vs KO) are compared within each staining panel and treatment, with the same
  pairwise statistics (`compare = cell_type` in `stats_pairwise.csv`; wild type is the reference).

- **Reference condition.** One condition, chosen as cell type + treatment (default: the wild type and the
  control, e.g. WT · MOCK), is the common reference. Every other condition, of every cell type, is compared
  with it (`compare = vs_reference` in `stats_pairwise.csv`), in addition to the within-cell-type and
  within-treatment comparisons.
- **High cells.** Per channel, a cell is "high" when its nuclear, cytoplasmic or whole-cell intensity
  (background-subtracted) is above a percentile (default 95th) of the reference condition's cells of the
  same staining panel. The cut-offs are in `high_thresholds.csv`; per cell, `<ch>_high_nucleus`,
  `_high_cytoplasm` and `_high_cell` are 100 or 0. Per replicate, the outputs are:
  - the % of high cells (`high_*` metrics in the summaries, tested like the others);
  - `<panel>_high_cells.png`, a 100 % stacked bar per condition.

  This captures phenotypes that live in a subpopulation (e.g. more very bright cells), which a median hides.

With fewer than two replicates per treatment no test is performed and only effect sizes are reported.
**Figures** (`superplots/`), per staining panel, with all cell types in the same graph (treatments on the x
axis, cell types side by side, one color per cell type):
- `<panel>_spots_per_nucleus.png`: spots per nucleus for each channel, and double-positive spots per
  nucleus for each channel pair, each as a SuperPlot and as the fraction of
  nuclei in each count class (0, then classes scaled to the 95th percentile of counts, e.g. 1-5, 6-10, 11-20,
  21-30, ≥ 31), pooled over replicates (`spot_count_distribution.csv`);
- `<panel>_single_double_positive.png` (also in the overview and as `<pair>_n_a_only.png` /
  `<pair>_n_b_only.png` / `<pair>_n_double.png` in `<panel>/`): for each channel pair, A-only, B-only and double-positive spots per
  nucleus, each as a SuperPlot and as a count distribution. A-only + A spots overlapping B = all A spots;
  double-positive spots merge touching A and B spots into one object, so they can be fewer than the
  overlapping A spots;
- `<panel>_spot_categories.png` (`spot_categories.csv`), per channel pair:
  - top: a 100 % stacked bar per condition, the share of all spots that are double-positive, A-only or
    B-only (computed per replicate, averaged over replicates);
  - bottom: mean double-positive, A-only and B-only spots per nucleus, each on its own axis, with error
    bars. With >= 2 replicates per condition the error is the SD of the replicate means. With one
    replicate per condition it is the SEM over nuclei, which describes one well and is not a
    replicate-level error; the figure title says which is shown;
- `<panel>_spot_area_per_nucleus.png`: for each channel, the average spot area per nucleus (size at half
  maximum; nuclei without spots are not included);
- `<panel>_overview.png`, and one PNG per metric in `<panel>/`: nuclei as a beeswarm colored by replicate,
  each replicate's median as a large dot of the same color with area proportional to its number of nuclei,
  and mean ± SD of the replicate medians in black. The x axis has two levels: cell type under each column,
  treatment underlined beneath. Replicate labels from a Replicate column keep their color across
  conditions; with wells as replicates they are numbered 1, 2, 3 within each condition.

**Order.** Treatments: controls (MOCK, DMSO, ...) first, then in plate-map order (the first well, row by row,
in which each appears), so MOCK, ML792, MG132, ML792+MG132 laid out down a column plot in that order. Cell
types: wild type (WT, parental, ...) first, then plate-map order. Both can be set explicitly
(`treatment_order`, `cell_type_order`).

## Outputs (`<output folder>/analysis/`)

| File | Content |
|---|---|
| `nuclei.csv` | one row per nucleus: identity, QC, shape, DAPI, and per-channel intensity / spot / dispersal measures |
| `spots.csv` | one row per spot |
| `colocalization.csv` | one row per nucleus × channel pair: pixel and object colocalization |
| `fields.csv` | per field: nuclei counts, background, saturation |
| `spot_categories.csv` | per channel pair and condition: mean double-/single-positive spots per nucleus with error, share of all spots |
| `spot_count_distribution.csv` | fraction of nuclei per spot-count class, per channel, cell type and treatment |
| `replicate_summary.csv`, `condition_summary.csv` | medians per replicate; mean ± SD of replicate medians per condition |
| `stats_omnibus.csv`, `stats_pairwise.csv` | Welch ANOVA (+ BH q); Games-Howell, Hedges' g, Cliff's delta (treatments within cell type, cell types within treatment) |
| `superplots/` | regenerated on every run (figures of earlier runs are removed): spots-per-nucleus figure, SuperPlot overview and per-metric SuperPlots, per staining panel |
| `qc/` | per-field overlay: nuclei (white = analysed, red = QC-flagged) and spots per channel |
| `masks/` | cached nuclear label images |
| `analysis_manifest.json` | parameters, software versions, timings |

## Validation

- **Synthetic ground truth** (`tests/test_analysis.py`):
  - 9 nuclei with 3–11 spots each: nuclei and spot counts are recovered exactly;
  - nuclei without spots yield no spots;
  - the fraction of shared spots is recovered exactly;
  - colocalized vs independent channels are separated (r > 0.9 vs < 0.5);
  - border nuclei are excluded and the size window is applied.
- **Reference implementations:**
  - Pearson's r and Manders (M1, tM1, tM2 at the Costes thresholds) equal scikit-image's `pearson_corr_coeff`
    and `manders_coloc_coeff`;
  - the below-threshold correlation at the Costes threshold is ≈ 0;
  - Games-Howell (t, df, p) and Welch's ANOVA (F, p) match pingouin 0.7.0;
  - BH q-values match `scipy.stats.false_discovery_control`.
- **Real data:** 40 fields of RPE1 cells (SUMO2/3–AF594, PML–AF488, DAPI; 0.086 µm/px):
  - 2,019 nuclei segmented by Cellpose on a GPU, 1,868 passing QC, in 90 s, of which 45 s was segmentation
    including model loading;
  - in untreated wells PML–SUMO2/3 Pearson's r was 0.77 and 71 % of PML bodies overlapped SUMO2/3 foci,
    consistent with SUMO's known enrichment in PML bodies.

## Performance

Measurement runs in parallel over fields (one process per CPU core, minus two). Cellpose runs on the GPU in the
main process: about 1 s per 4096 × 3008 field on an RTX 5070 Ti, versus ~10 s on the CPU.
