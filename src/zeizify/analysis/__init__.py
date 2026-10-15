"""
Zeizify Analyze - per-nucleus quantification of nuclear bodies and colocalization, run on
the per-field TIFFs Zeizify extracts.

    nuclei      DAPI segmentation (Cellpose, or classical threshold + watershed)
    spots       speckles / bodies per signal channel inside each nucleus
    measure     per-nucleus counts, sizes, shapes, intensities, threshold-free dispersal
    coloc       pixel-based (Pearson, Spearman, Manders + Costes, Li ICQ) and object-based
    stats       replicate-level statistics (Welch ANOVA, Games-Howell) and SuperPlots

See docs/analysis.md for the methods and the sources each step follows.
"""

from .params import AnalysisParams

__all__ = ["AnalysisParams"]
