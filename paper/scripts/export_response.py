#!/usr/bin/env python3
"""Issue #3: response tables from Samrat et al. (2025), Table S10, for Fig. 5A-B.

  response_by_class_adduct.tsv  Fig. 5A: RIE (relative to LPE) per standard, class and adduct
  ie_loco_predictions.tsv       Fig. 5B: measured vs predicted log10 RIE per standard
                                (leave-one-class-out: each class predicted by a model that never saw it)
  ie_model_stats.tsv            R2, RMSE and MAE (log10 units) per model; best model marked

The structure builder and descriptors are imported from the SLA repository's
methods/08_climgrass/scripts/rie_prediction.py (path and SHA recorded in
provenance). Unlike that script, descriptors are screened on the standards
only (no archaeal prediction set), so the statistics are an S10-only
re-evaluation, not a copy of the earlier numbers.
Positive-mode adducts only (M+H, M+NH4), as in the original model.
Run with an environment that has rdkit and scikit-learn.
"""
from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge

sys.path.insert(0, str(Path(__file__).resolve().parent))
import provenance  # noqa: E402

SEED = 20260808
MECHANISTIC = ["MolLogP", "TPSA", "NumHAcceptors", "NumHDonors", "NumRotatableBonds", "FractionCSP3",
               "MolWt", "MolMR", "NOCount", "NHOHCount", "LabuteASA", "BalabanJ"]


def load_builder(path: Path):
    spec = importlib.util.spec_from_file_location("rie_prediction", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def loco(X: np.ndarray, y: np.ndarray, classes: np.ndarray, make) -> np.ndarray:
    pred = np.full(len(y), np.nan)
    for held in sorted(set(classes)):
        m = classes == held
        model = make()
        model.fit(X[~m], y[~m])
        pred[m] = model.predict(X[m])
    return pred


def stats(pred: np.ndarray, y: np.ndarray) -> dict:
    err = pred - y
    return {"r2": float(np.corrcoef(pred, y)[0, 1] ** 2),
            "rmse_log10": float(np.sqrt(np.mean(err ** 2))),
            "mae_log10": float(np.mean(np.abs(err))),
            "fold_error": float(10 ** np.sqrt(np.mean(err ** 2)))}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rie", required=True)
    ap.add_argument("--builder", required=True, help="SLA repo methods/08_climgrass/scripts/rie_prediction.py")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rie = pd.read_csv(a.rie)

    fig_a = rie.copy()
    fig_a["n_in_group"] = fig_a.groupby(["class", "adduct"])["RIE_LPE"].transform("size")
    fig_a["group_mean_RIE_LPE"] = fig_a.groupby(["class", "adduct"])["RIE_LPE"].transform("mean")
    fig_a.sort_values(["adduct", "group_mean_RIE_LPE"]).to_csv(out / "response_by_class_adduct.tsv",
                                                               sep="\t", index=False)

    b = load_builder(Path(a.builder))
    tr = rie[rie["adduct"].isin(b.POS_ADDUCTS)].copy()
    tr["chains"] = [b.parse_species(c, s) for c, s in zip(tr["class"], tr["species"])]
    tr["smiles"] = [b.build_smiles(c, ch) if ch else None for c, ch in zip(tr["class"], tr["chains"])]
    skipped = tr[tr["smiles"].isna()]
    tr = tr[tr["smiles"].notna()].reset_index(drop=True)
    desc = b.descriptor_frame(tr["smiles"].tolist())
    desc = desc.loc[:, desc.notna().all() & (desc.std() > 0) & (desc.abs().max() < 1e8)]
    for ad in b.POS_ADDUCTS:
        desc[f"adduct_{ad}"] = tr["adduct"].eq(ad).astype(float).to_numpy()
    desc["quaternary_N"] = tr["smiles"].str.count(r"\[N\+\]").astype(float).to_numpy()
    mech = [c for c in MECHANISTIC if c in desc.columns] + [c for c in desc.columns
                                                            if c.startswith("adduct_") or c == "quaternary_N"]

    def z(d):
        return ((d - d.mean()) / d.std().replace(0, 1.0)).to_numpy()

    X, Xm = z(desc), z(desc[mech])
    y, cls = tr["logRIE_LPE"].astype(float).to_numpy(), tr["class"].to_numpy()
    models = {
        "random_forest": (X, lambda: RandomForestRegressor(n_estimators=500, min_samples_leaf=2, random_state=SEED)),
        "ridge": (X, lambda: Ridge(alpha=10.0)),
        "rf_mechanistic": (Xm, lambda: RandomForestRegressor(n_estimators=500, min_samples_leaf=2, random_state=SEED)),
        "ridge_mechanistic": (Xm, lambda: Ridge(alpha=1.0)),
    }
    rows, preds = [], {}
    for name, (XX, make) in models.items():
        preds[name] = loco(XX, y, cls, make)
        rows.append({"model": name, "n_standards": len(y), "n_classes": len(set(cls)), **stats(preds[name], y)})
    st = pd.DataFrame(rows).sort_values("rmse_log10")
    best = st.iloc[0]["model"]
    st["best"] = st["model"] == best
    st.round(4).to_csv(out / "ie_model_stats.tsv", sep="\t", index=False)
    pr = tr[["class", "species", "adduct", "logRIE_LPE"]].copy()
    for name in models:
        pr[f"pred_{name}"] = preds[name]
    pr["pred_best"] = preds[best]
    pr.to_csv(out / "ie_loco_predictions.tsv", sep="\t", index=False)
    if len(skipped):
        skipped[["class", "species", "adduct"]].to_csv(out / "ie_skipped_no_structure.tsv", sep="\t", index=False)

    sla_sha = subprocess.run(["git", "-C", str(Path(a.builder).parent), "rev-parse", "HEAD"],
                             capture_output=True, text=True).stdout.strip()
    provenance.write(out, {"rie_table_s10": a.rie, "rie_prediction_builder": a.builder},
                     {"cv": "leave-one-class-out", "adducts": list(b.POS_ADDUCTS), "seed": SEED,
                      "builder_repo_sha": sla_sha, "best_model": best},
                     notes="Matrix effect not computed: no post-extraction spikes in the ClimGrass sequence.")


if __name__ == "__main__":
    main()
