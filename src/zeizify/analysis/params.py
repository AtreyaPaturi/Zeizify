"""All tunable analysis parameters in one place, with physical units, so a run is fully described
by this dataclass (it is written verbatim to analysis_manifest.json)."""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field


@dataclass
class AnalysisParams:
    # --- channels -------------------------------------------------------------------------------
    nuclear_channel: str = ""  # "" = auto (DAPI / Hoechst)
    spot_channels: list[str] = field(default_factory=list)  # [] = every non-nuclear channel

    # --- nuclei -----------------------------------------------------------------------------------
    nucleus_method: str = "auto"  # "auto" (Cellpose if installed, else classical), "cellpose", "classical"
    nucleus_diameter_um: float = 14.0  # typical nuclear diameter; sets scale for both methods
    nucleus_min_area_frac: float = 0.25  # size window, as a fraction of the typical nuclear area ...
    nucleus_max_area_frac: float = 4.0  # ... (rejects debris and clumps)
    exclude_border_nuclei: bool = True  # nuclei cut by the image edge are incomplete
    use_gpu: bool = True  # Cellpose only

    # --- nucleus quality control (flagged nuclei are kept in tables, excluded from statistics) -----
    qc_dapi_low: float = 0.5  # mean DAPI below this x field median -> flagged (out of focus / debris)
    qc_dapi_high: float = 1.6  # above this x field median -> flagged (mitotic / apoptotic, condensed)
    qc_max_saturated_frac: float = 0.01  # > 1 % saturated pixels in any channel -> flagged

    # --- spots (speckles / nuclear bodies) ----------------------------------------------------------
    spot_min_diameter_um: float = 0.25
    spot_max_diameter_um: float = 3.0
    spot_sensitivity: float = 6.0  # noise floor: median + k * MAD of the top-hat nucleoplasm (lower = more spots)
    spot_declump: bool = True  # split touching spots at intensity saddles (watershed)

    # --- colocalization --------------------------------------------------------------------------
    coloc_randomizations: int = 0  # Costes block-scrambling significance test; 0 = off (e.g. 100)
    psf_um: float = 0.25  # block size for the randomization (about the optical resolution)

    # --- run ----------------------------------------------------------------------------------------
    workers: int = max(1, (os.cpu_count() or 2) - 2)
    qc_overlays: bool = True
    replicate_column: str = ""  # "" = plate-map "Replicate" column if present, else the well
    seed: int = 0

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "AnalysisParams":
        known = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
        return cls(**known)

    def nucleus_area_um2(self) -> float:
        import math

        return math.pi * (self.nucleus_diameter_um / 2) ** 2
