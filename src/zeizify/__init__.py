"""Zeizify - per-well extraction of Zeiss CZI multi-well plate acquisitions.

Map every scene of a multi-well-plate CZI to its well, max-project each field's
Z-stack per channel, write analysis-ready TIFFs grouped by antibody, and choose
reference contrast ranges from control wells - from a desktop GUI or the CLI.
"""

__version__ = "1.4.0"
