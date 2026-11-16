"""
Source decomposition — multi-method suite.

Methods:
  1. nnls_intensity    - scipy NNLS on L1-normalized intensity
  2. standard_bc       - Bray-Curtis similarity (1 - BC), normalized to proportion
  3. enriched_only_bc  - BC using only SIMPER-enriched features per phylum
  4. fc_weighted_bc    - BC weighted by SIMPER fold-change per phylum
  5. bayesian_log      - NumPyro log-space Dirichlet mixing (optional, slow)

All methods take numpy arrays (ref: n_phyla × n_features, soil: n_samples × n_features)
and return a composition matrix (n_samples × n_phyla) summing to 1 per row.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import nnls as scipy_nnls
from scipy.spatial.distance import braycurtis


PHYLUM_KINGDOM = {
    'Actinomycetota': 'Bacteria', 'Bacillota': 'Bacteria',
    'Cyanobacteriota': 'Bacteria', 'Pseudomonadota': 'Bacteria',
    'Euryarchaeota': 'Archaea', 'Methanobacteriota': 'Archaea',
    'Crenarchaeota': 'Archaea', 'Thermoproteota': 'Archaea',
    'Ascomycota': 'Fungi', 'Basidiomycota': 'Fungi', 'Mucoromycota': 'Fungi',
    'Arthropoda': 'Animalia', 'Nematoda': 'Animalia', 'Mollusca': 'Animalia',
    'Bryophyta': 'Plantae', 'Chlorophyta': 'Plantae', 'Marchantiophyta': 'Plantae',
    'Trachaeophyta': 'Plantae', 'Charophyta': 'Plantae', 'Magnoliophyta': 'Plantae',
    'Amoebozoa': 'Protozoa',
}
EXCLUDE_PHYLA = {'Virus', 'Mixed', 'Rootnodules', ''}


# -----------------------------------------------------------------------------
# Build phylum reference (uses mean of atlas samples per phylum)
# -----------------------------------------------------------------------------
def build_phylum_reference_array(
    atlas_intensity_df: pd.DataFrame,   # features × samples
    sample_phylum: dict,
    min_samples: int = 2,
) -> tuple[np.ndarray, list[str], list[str]]:
    """Return (ref_array[n_phyla × n_features], phyla, feature_ids)."""
    mapped = {c: p for c, p in sample_phylum.items()
              if c in atlas_intensity_df.columns and p not in EXCLUDE_PHYLA}
    counts = defaultdict(int)
    for p in mapped.values(): counts[p] += 1
    phyla = sorted([p for p, n in counts.items() if n >= min_samples])

    feature_ids = list(atlas_intensity_df.index)
    ref = np.zeros((len(phyla), len(feature_ids)))
    for i, p in enumerate(phyla):
        cols = [c for c, pp in mapped.items() if pp == p]
        ref[i] = atlas_intensity_df[cols].values.astype(float).mean(axis=1)
    return ref, phyla, feature_ids


def build_soil_matrix_array(soil_intensity_df: pd.DataFrame, feature_ids: list[str]) -> tuple[np.ndarray, list[str]]:
    """Align soil data to feature_ids order. Returns (soil_array[n_samples × n_features], sample_cols)."""
    common = [f for f in feature_ids if f in soil_intensity_df.index]
    aligned = soil_intensity_df.reindex(feature_ids).fillna(0)
    return aligned.values.T, list(aligned.columns)


# -----------------------------------------------------------------------------
# Method 1: NNLS (intensity-based)
# -----------------------------------------------------------------------------
def decompose_nnls(ref: np.ndarray, soil: np.ndarray) -> np.ndarray:
    """ref: n_phyla × n_features; soil: n_samples × n_features.
    Normalize both per-sample/per-phylum to L1, then NNLS."""
    # L1 normalize phylum refs (each phylum sums to 1 across features)
    ref_sum = ref.sum(axis=1, keepdims=True)
    ref_sum[ref_sum == 0] = 1.0
    R = (ref / ref_sum).T  # n_features × n_phyla

    n_samples, n_phyla = soil.shape[0], ref.shape[0]
    comp = np.zeros((n_samples, n_phyla))
    for i in range(n_samples):
        s = soil[i]
        total = s.sum()
        if total == 0: continue
        s_norm = s / total
        w, _ = scipy_nnls(R, s_norm)
        if w.sum() > 0:
            comp[i] = w / w.sum()
    return comp


# -----------------------------------------------------------------------------
# Method 2: Standard BC
# -----------------------------------------------------------------------------
def decompose_standard_bc(ref: np.ndarray, soil: np.ndarray) -> np.ndarray:
    n_samples, n_phyla = soil.shape[0], ref.shape[0]
    comp = np.zeros((n_samples, n_phyla))
    for i in range(n_samples):
        if soil[i].sum() == 0: continue
        sims = np.zeros(n_phyla)
        for p in range(n_phyla):
            if ref[p].sum() == 0: continue
            sims[p] = max(0, 1 - braycurtis(soil[i], ref[p]))
        if sims.sum() > 0:
            comp[i] = sims / sims.sum()
    return comp


# -----------------------------------------------------------------------------
# Build direction masks + FC weights from SIMPER atlas
# -----------------------------------------------------------------------------
def build_direction_masks(
    simper_full: pd.DataFrame,
    phyla: list[str],
    feature_ids: list[str],
) -> tuple[dict, dict]:
    """Return enriched_masks[phylum] = bool array, fc_weights[phylum] = float array."""
    fid_to_idx = {f: i for i, f in enumerate(feature_ids)}
    n_features = len(feature_ids)

    enriched_masks = {}
    fc_weights = {}

    for p in phyla:
        mask = np.zeros(n_features, dtype=bool)
        weights = np.ones(n_features) * 0.01
        p_simper = simper_full[simper_full['phylum'] == p]
        for _, row in p_simper.iterrows():
            fid = row['feature_id']
            if fid not in fid_to_idx: continue
            idx = fid_to_idx[fid]
            if row.get('direction') == 'enriched':
                mask[idx] = True
                fc = row.get('fold_change', 1.0)
                if np.isfinite(fc) and fc > 0:
                    weights[idx] = min(fc, 100)
                else:
                    weights[idx] = 100
            else:
                weights[idx] = 0.01
        enriched_masks[p] = mask
        fc_weights[p] = weights
    return enriched_masks, fc_weights


# -----------------------------------------------------------------------------
# Method 3: Enriched-only BC
# -----------------------------------------------------------------------------
def decompose_enriched_only_bc(ref, soil, enriched_masks, phyla):
    n_samples, n_phyla, n_features = soil.shape[0], len(phyla), ref.shape[1]
    comp = np.zeros((n_samples, n_phyla))
    for i in range(n_samples):
        if soil[i].sum() == 0: continue
        sims = np.zeros(n_phyla)
        for p_idx, p in enumerate(phyla):
            mask = enriched_masks.get(p, np.ones(n_features, dtype=bool))
            if mask.sum() < 10: continue
            s_m = soil[i][mask]
            r_m = ref[p_idx][mask]
            if s_m.sum() > 0 and r_m.sum() > 0:
                sims[p_idx] = max(0, 1 - braycurtis(s_m, r_m))
        if sims.sum() > 0:
            comp[i] = sims / sims.sum()
    return comp


# -----------------------------------------------------------------------------
# Method 4: FC-weighted BC
# -----------------------------------------------------------------------------
def _weighted_bc(a, b, w):
    wn = w / (w.sum() + 1e-10)
    num = np.sum(wn * np.abs(a - b))
    den = np.sum(wn * (a + b))
    return num / (den + 1e-10) if den > 0 else 1.0


def decompose_fc_weighted_bc(ref, soil, fc_weights, phyla):
    n_samples, n_phyla, n_features = soil.shape[0], len(phyla), ref.shape[1]
    comp = np.zeros((n_samples, n_phyla))
    for i in range(n_samples):
        if soil[i].sum() == 0: continue
        sims = np.zeros(n_phyla)
        for p_idx, p in enumerate(phyla):
            w = fc_weights.get(p, np.ones(n_features))
            sims[p_idx] = max(0, 1 - _weighted_bc(soil[i], ref[p_idx], w))
        if sims.sum() > 0:
            comp[i] = sims / sims.sum()
    return comp


# -----------------------------------------------------------------------------
# Method 5: Bayesian log-space Dirichlet mixing (NumPyro, slow)
# -----------------------------------------------------------------------------
def decompose_bayesian_log(ref, soil, num_warmup=500, num_samples=500):
    """Geometric (log-space) mixing model with Dirichlet prior on alpha.

    log(soil + 1) ~ Normal(alpha @ log(ref + 1), sigma)

    Slow (~1-5 min per soil sample on CPU). Returns mean posterior composition.
    """
    try:
        import jax
        import jax.numpy as jnp
        import numpyro
        import numpyro.distributions as dist
        from numpyro.infer import MCMC, NUTS
    except ImportError:
        print("  [bayesian_log] NumPyro not available; skipping.")
        return None

    n_samples, n_phyla, n_features = soil.shape[0], ref.shape[0], ref.shape[1]
    comp = np.zeros((n_samples, n_phyla))

    log_ref = np.log(ref + 1.0)  # n_phyla × n_features

    def model(log_y, log_R, n_phyla):
        alpha = numpyro.sample("alpha", dist.Dirichlet(jnp.ones(n_phyla)))
        sigma = numpyro.sample("sigma", dist.HalfNormal(1.0))
        mu = jnp.dot(alpha, log_R)
        numpyro.sample("obs", dist.Normal(mu, sigma), obs=log_y)

    for i in range(n_samples):
        if soil[i].sum() == 0:
            continue
        log_y = np.log(soil[i] + 1.0)
        rng_key = jax.random.PRNGKey(42 + i)
        kernel = NUTS(model)
        mcmc = MCMC(kernel, num_warmup=num_warmup, num_samples=num_samples, progress_bar=False)
        mcmc.run(rng_key, log_y=jnp.array(log_y), log_R=jnp.array(log_ref), n_phyla=n_phyla)
        samples = mcmc.get_samples()
        alpha_mean = np.asarray(samples["alpha"].mean(axis=0))
        comp[i] = alpha_mean
    return comp


# -----------------------------------------------------------------------------
# Phylum -> Kingdom rollup
# -----------------------------------------------------------------------------
def phylum_to_kingdom_array(
    comp: np.ndarray,     # n_samples × n_phyla
    phyla: list[str],
) -> tuple[np.ndarray, list[str]]:
    kingdoms = sorted({PHYLUM_KINGDOM.get(p, 'Unknown') for p in phyla})
    out = np.zeros((comp.shape[0], len(kingdoms)))
    for p_idx, p in enumerate(phyla):
        k = PHYLUM_KINGDOM.get(p, 'Unknown')
        k_idx = kingdoms.index(k)
        out[:, k_idx] += comp[:, p_idx]
    return out, kingdoms


# -----------------------------------------------------------------------------
# Unified runner — runs ALL methods on one (ref, soil) pair
# -----------------------------------------------------------------------------
def run_all_methods(
    ref: np.ndarray,
    soil: np.ndarray,
    phyla: list[str],
    simper_full: pd.DataFrame,
    feature_ids: list[str],
    include_bayesian: bool = False,
) -> dict[str, np.ndarray]:
    """Run all decomposition methods. Returns dict[method_name] = composition array."""
    results = {}
    print("    → NNLS...", end=" ", flush=True)
    results["nnls"] = decompose_nnls(ref, soil); print("done")

    print("    → Standard BC...", end=" ", flush=True)
    results["standard_bc"] = decompose_standard_bc(ref, soil); print("done")

    # Build direction masks once for enriched/fc methods
    enriched_masks, fc_weights = build_direction_masks(simper_full, phyla, feature_ids)

    print("    → Enriched-only BC...", end=" ", flush=True)
    results["enriched_only_bc"] = decompose_enriched_only_bc(ref, soil, enriched_masks, phyla); print("done")

    print("    → FC-weighted BC...", end=" ", flush=True)
    results["fc_weighted_bc"] = decompose_fc_weighted_bc(ref, soil, fc_weights, phyla); print("done")

    if include_bayesian:
        print("    → Bayesian log-space (slow)...", end=" ", flush=True)
        bayes = decompose_bayesian_log(ref, soil)
        if bayes is not None:
            results["bayesian_log"] = bayes
            print("done")
        else:
            print("skipped")

    return results
