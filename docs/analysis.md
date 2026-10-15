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

**Nucleus QC.** A nucleus is flagged (`qc_pass = False`) when its mean DAPI intensity is < 0.5 × or > 1.6 ×
the median of its field, or when > 1 % of its pixels are saturated in any channel. The first case catches
out-of-focus nuclei and debris; the second catches condensed mitotic or apoptotic nuclei. Flagged nuclei stay
in `nuclei.csv` and are excluded from summaries and statistics. The 1.6 default was set on RPE1 images, where
mitotic and post-mitotic nuclei measured 1.88–1.94 × and interphase nuclei ≤ 1.36 ×.

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
4. **Size window** on area-equivalent diameter (default 0.25–3 µm).

Each spot is assigned to its parent nucleus, the equivalent of CellProfiler's RelateObjects.

## 3. Measurements

**Per spot** (`spots.csv`): area, equivalent diameter, eccentricity, solidity, mean / max / integrated
intensity, integrated top-hat intensity, centroid.

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
- the fraction of A spots overlapping any B spot, and vice versa;
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

With fewer than two replicates per treatment no test is performed and only effect sizes are reported.
**Figures** (`superplots/`), per staining panel, with all cell types in the same graph (treatments on the x
axis, cell types side by side, one color per cell type):
- `<panel>_spots_per_nucleus.png`: for each channel, spots per nucleus as a SuperPlot, and the fraction of
  nuclei in each count class (0, then classes scaled to the 95th percentile of counts, e.g. 1-5, 6-10, 11-20,
  21-30, ≥ 31), pooled over replicates (`spot_count_distribution.csv`);
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
| `spot_count_distribution.csv` | fraction of nuclei per spot-count class, per channel, cell type and treatment |
| `replicate_summary.csv`, `condition_summary.csv` | medians per replicate; mean ± SD of replicate medians per condition |
| `stats_omnibus.csv`, `stats_pairwise.csv` | Welch ANOVA (+ BH q); Games-Howell, Hedges' g, Cliff's delta (treatments within cell type, cell types within treatment) |
| `superplots/` | spots-per-nucleus figure, SuperPlot overview and per-metric SuperPlots, per staining panel |
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
