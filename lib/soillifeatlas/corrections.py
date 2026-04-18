"""
Modular correction functions for quantification pipeline.

Each function takes a feature × sample intensity matrix and returns a corrected
matrix of the same shape. Functions are designed to be composed via pipeline.py.

All corrections preserve NaN/zero handling: missing features stay missing.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import numpy as np
import pandas as pd


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
@dataclass
class CorrectionConfig:
    """Toggleable correction layers.

    Order of application (when multiple enabled):
      1. biomass_normalization  -> per-sample scale by extracted tissue amount
      2. IS_normalization       -> per-sample scale by IS intensity
      3. RIE_correction         -> per-feature scale by class ionization efficiency
      4. external_calibration   -> per-feature scale by external calibration curve
    """
    biomass_normalization: bool = False
    IS_normalization: bool = False
    RIE_correction: bool = False
    external_calibration: bool = False

    # Parameters
    IS_reference_compound: str = "LPE_18d7"  # which IS to anchor on
    IS_spiked_pmol: float = 100.0            # pmol of IS added per sample
    RIE_reference_class: str = "LPE"         # RIE values are relative to this
    fallback_RIE: float = 1.0                # For classes not in S10 table
    RIE_floor: float = 0.05                   # Min RIE allowed (caps amplification at 1/0.05 = 20×)
    RIE_ceiling: float = 100.0                # Max RIE allowed (caps suppression at 100×)

    # Reference-quality filters (applied during reference profile construction)
    restrict_archaea_to_archlips: bool = False  # Use only ArchLips-validated features for Archaea phyla
    aggregate_archaea: bool = False              # Pool all 4 archaeal phyla into one "Archaea" reference
    restrict_to_annotated_only: bool = False     # Keep only features with LS or ArchLips annotation

    def label(self) -> str:
        parts = []
        if self.biomass_normalization:  parts.append("biomass")
        if self.IS_normalization:        parts.append("IS")
        if self.RIE_correction:          parts.append("RIE")
        if self.external_calibration:    parts.append("EXT")
        if self.restrict_archaea_to_archlips: parts.append("archlipsArch")
        if self.aggregate_archaea:       parts.append("aggArch")
        if self.restrict_to_annotated_only: parts.append("annotOnly")
        return "+".join(parts) if parts else "raw"


# -----------------------------------------------------------------------------
# Layer 1: Biomass normalization
# -----------------------------------------------------------------------------
def apply_biomass_normalization(
    intensity: pd.DataFrame,
    biomass_mg: pd.Series,
) -> pd.DataFrame:
    """Divide each sample column by biomass in mg. Result units: intensity / mg."""
    # Align biomass to sample columns
    aligned = biomass_mg.reindex(intensity.columns).fillna(biomass_mg.median())
    out = intensity.div(aligned.values, axis=1)
    return out


# -----------------------------------------------------------------------------
# Layer 2: IS scaling
# -----------------------------------------------------------------------------
def apply_IS_normalization(
    intensity: pd.DataFrame,
    IS_intensity_per_sample: pd.Series,  # IS signal in each sample
    IS_spiked_pmol: float = 100.0,
    global_median: float | None = None,
) -> pd.DataFrame:
    """Scale each sample so its IS signal equals the global median (relative)
    and then convert to pmol-equivalent units using IS_spiked_pmol.

    corrected = raw / IS_intensity × IS_spiked_pmol

    If IS_intensity is 0 or NaN for a sample, fall back to global median scaling.
    """
    if global_median is None:
        pos = IS_intensity_per_sample[IS_intensity_per_sample > 0]
        global_median = float(pos.median()) if len(pos) else 1.0

    # Factor per sample
    factor = IS_intensity_per_sample.copy()
    factor = factor.where(factor > 0, global_median)  # fill missing IS with median
    factor = IS_spiked_pmol / factor  # pmol per intensity-count

    aligned = factor.reindex(intensity.columns).fillna(IS_spiked_pmol / global_median)
    out = intensity.mul(aligned.values, axis=1)
    return out


# -----------------------------------------------------------------------------
# Layer 3: RIE correction
# -----------------------------------------------------------------------------
def apply_RIE_correction(
    intensity: pd.DataFrame,
    feature_class: pd.Series,        # class per feature (e.g., "PC", "PE", "TG")
    feature_adduct: pd.Series,       # adduct per feature (e.g., "[M+H]+")
    rie_lookup: pd.DataFrame,        # columns: class, adduct, RIE_LPE
    fallback_RIE: float = 1.0,
    rie_floor: float = 0.05,
    rie_ceiling: float = 100.0,
) -> pd.DataFrame:
    """Divide each feature row by its class-specific RIE relative to LPE.

    corrected = intensity / RIE(class, adduct)

    Features with no class assignment get fallback_RIE (= no correction).
    """
    # Build lookup dict — (class_upper, adduct_stripped) → mean RIE (averaged across species)
    rie_lookup_clean = rie_lookup.dropna(subset=["RIE_LPE"]).copy()
    rie_lookup_clean["class_u"] = rie_lookup_clean["class"].astype(str).str.upper()
    rie_lookup_clean["adduct_clean"] = (
        rie_lookup_clean["adduct"].astype(str)
        .str.replace("[", "", regex=False)
        .str.replace("]", "", regex=False)
        .str.strip()
    )
    # Average RIE across species within each (class, adduct)
    rie_key = rie_lookup_clean.groupby(["class_u", "adduct_clean"])["RIE_LPE"].mean().to_dict()

    def rie_for(cls, add):
        if pd.isna(cls) or pd.isna(add) or cls is None:
            return fallback_RIE
        c = str(cls).upper()
        a = str(add).replace("[", "").replace("]", "").strip()
        # Strip charge trailing symbols like +, -
        if a.endswith("+") or a.endswith("-"):
            a = a[:-1]
        v = rie_key.get((c, a))
        if v is not None and v > 0 and not pd.isna(v):
            return float(v)
        return fallback_RIE

    rie_per_feature = pd.Series(
        [rie_for(c, a) for c, a in zip(feature_class, feature_adduct)],
        index=intensity.index,
    )
    # Apply floor/ceiling to prevent extreme amplification or suppression
    rie_clipped = rie_per_feature.clip(lower=rie_floor, upper=rie_ceiling)
    out = intensity.div(rie_clipped.values, axis=0)
    return out, rie_clipped


# -----------------------------------------------------------------------------
# Layer 4: External calibration (placeholder for future)
# -----------------------------------------------------------------------------
def apply_external_calibration(
    intensity: pd.DataFrame,
    calibration_curves: pd.DataFrame,  # columns: class, slope, intercept
) -> pd.DataFrame:
    """Apply per-class linear calibration: conc = (intensity - intercept) / slope.
    Placeholder until user compiles external cal data."""
    raise NotImplementedError("External calibration not yet compiled. Skip for now.")


# -----------------------------------------------------------------------------
# Pipeline orchestrator
# -----------------------------------------------------------------------------
def correct_intensity(
    raw_intensity: pd.DataFrame,
    config: CorrectionConfig,
    IS_intensity_per_sample: pd.Series | None = None,
    biomass_mg: pd.Series | None = None,
    feature_class: pd.Series | None = None,
    feature_adduct: pd.Series | None = None,
    rie_lookup: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Apply configured correction layers in order. Returns corrected matrix + diagnostics."""
    out = raw_intensity.copy()
    diag = {"config_label": config.label(), "layers_applied": [], "shape": out.shape}

    if config.biomass_normalization:
        if biomass_mg is None:
            raise ValueError("biomass_mg required when biomass_normalization=True")
        out = apply_biomass_normalization(out, biomass_mg)
        diag["layers_applied"].append("biomass")

    if config.IS_normalization:
        if IS_intensity_per_sample is None:
            raise ValueError("IS_intensity_per_sample required when IS_normalization=True")
        out = apply_IS_normalization(out, IS_intensity_per_sample, config.IS_spiked_pmol)
        diag["layers_applied"].append("IS")
        diag["IS_fraction_samples_with_signal"] = float((IS_intensity_per_sample > 0).mean())

    if config.RIE_correction:
        if feature_class is None or feature_adduct is None or rie_lookup is None:
            raise ValueError("feature_class, feature_adduct, rie_lookup required for RIE")
        out, rie_vec = apply_RIE_correction(out, feature_class, feature_adduct, rie_lookup, config.fallback_RIE)
        diag["layers_applied"].append("RIE")
        diag["RIE_fraction_features_corrected"] = float((rie_vec != config.fallback_RIE).mean())

    if config.external_calibration:
        raise NotImplementedError("External calibration layer not yet implemented.")

    return out, diag
