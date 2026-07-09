"""
Reconstruct `Raw_feature_tables.xlsx` from the compiled Griffin profile CSVs.

Background
----------
The notebook `PE_source_code.ipynb` loads a workbook `Raw_feature_tables.xlsx`
with three sheets that are merged with the metadata on `sample_ID`:

    Griffin_scATAC_features = 'Raw feature Table 1'   # scATAC-derived tissue sites
    Griffin_DNase_features  = 'Raw feature Table 2'   # DNase (DHS) tissue sites
    Griffin_TFBS_features   = 'Raw feature Table 3'   # transcription-factor sites

That workbook is missing.  Per private communication the underlying data are the
four compiled Griffin runs in `data/`:

    Nucleosome_TFBS_profiles_Griffin_run_compiled_v1.csv.gz       (fragments 120-180 bp)
    Nucleosome_Tissue_profiles_Griffin_run_compiled_v1.csv.gz     (fragments 120-180 bp)
    Subnucleosome_TFBS_profiles_Griffin_run_compiled_v1.csv.gz    (fragments  35-80  bp)
    Subnucleosome_Tissue_profiles_Griffin_run_compiled_v1.csv.gz  (fragments  35-80  bp)

Each row is one Griffin coverage profile: 132 columns of GC-corrected normalized
coverage sampled every 15 bp from -990 .. +975 bp around the site, followed by
Griffin metadata columns (..., mean_coverage, central_coverage, ..., site_name,
..., sample).

Feature definitions (Adil et al. 2025, Nature Medicine, Fig. 1c & Methods)
--------------------------------------------------------------------------
    MCV     mean coverage around the site, +/- 1000 bp
            -> mean of all 132 position columns (equals Griffin's `mean_coverage`)
    NDR-30  normalized coverage at the nucleosome-depleted region, +/- 30 bp
            -> mean of position columns with |pos| <= 30  (bins -30,-15,0,15,30)
    NDR-75  normalized coverage at the nucleosome-depleted region, +/- 75 bp
            -> mean of position columns with |pos| <= 75  (11 bins)

Feature naming used by the notebook
-----------------------------------
    <site>_<fragment>_<metric>   e.g. Placental_120_180_MCV, GRHL2_120_180_NDR-75,
                                       Endothelial_35_80_NDR-30
For TFBS the Griffin `site_name` looks like `GRHL2.hg38.10000_120_180`; the
`.hg38.10000` token is stripped so it becomes `GRHL2_120_180`.

The tissue sheets are split into scATAC vs DNase using the bed-file directories
under `data/Bed_files_hg38/`.
"""

import os
import numpy as np
import pandas as pd

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
BED_DIR = os.path.join(DATA_DIR, "Bed_files_hg38")
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Raw_feature_tables.xlsx")

NUC_TFBS = os.path.join(DATA_DIR, "Nucleosome_TFBS_profiles_Griffin_run_compiled_v1.csv.gz")
NUC_TISSUE = os.path.join(DATA_DIR, "Nucleosome_Tissue_profiles_Griffin_run_compiled_v1.csv.gz")
SUB_TFBS = os.path.join(DATA_DIR, "Subnucleosome_TFBS_profiles_Griffin_run_compiled_v1.csv.gz")
SUB_TISSUE = os.path.join(DATA_DIR, "Subnucleosome_Tissue_profiles_Griffin_run_compiled_v1.csv.gz")


def _clean_site_name(site_name: str) -> str:
    """`GRHL2.hg38.10000_120_180` -> `GRHL2_120_180`; tissue names are unchanged."""
    return site_name.replace(".hg38.10000", "")


def _compute_features(csv_path: str) -> pd.DataFrame:
    """Read one compiled Griffin CSV and return a wide (sample_ID x feature) frame.

    Feature columns are `<clean_site>_<MCV|NDR-30|NDR-75>`.
    """
    df = pd.read_csv(csv_path)

    position_cols = [c for c in df.columns if c.lstrip("-").isdigit()]
    w30 = [c for c in position_cols if abs(int(c)) <= 30]
    w75 = [c for c in position_cols if abs(int(c)) <= 75]

    prof = df[position_cols].astype(np.float64)
    out = pd.DataFrame({"sample_ID": df["sample"].values})
    out["site"] = df["site_name"].map(_clean_site_name).values
    out["MCV"] = prof.mean(axis=1).values
    out["NDR-30"] = prof[w30].mean(axis=1).values
    out["NDR-75"] = prof[w75].mean(axis=1).values

    long = out.melt(id_vars=["sample_ID", "site"],
                    value_vars=["MCV", "NDR-30", "NDR-75"],
                    var_name="metric", value_name="value")
    long["feature"] = long["site"] + "_" + long["metric"]
    wide = long.pivot_table(index="sample_ID", columns="feature", values="value")
    wide.columns.name = None
    return wide.reset_index()


def _bed_site_names(subdir: str) -> set:
    """Site base names (without fragment suffix) taken from a bed-file directory."""
    path = os.path.join(BED_DIR, subdir)
    return {os.path.splitext(f)[0] for f in os.listdir(path) if f.endswith(".bed")}


def build_tables():
    print("Computing nucleosome (120-180 bp) tissue features ...")
    nuc_tissue = _compute_features(NUC_TISSUE)
    print("Computing sub-nucleosome (35-80 bp) tissue features ...")
    sub_tissue = _compute_features(SUB_TISSUE)
    print("Computing nucleosome (120-180 bp) TFBS features ...")
    nuc_tfbs = _compute_features(NUC_TFBS)
    print("Computing sub-nucleosome (35-80 bp) TFBS features ...")
    sub_tfbs = _compute_features(SUB_TFBS)

    tissue = pd.merge(nuc_tissue, sub_tissue, on="sample_ID", how="inner")
    tfbs = pd.merge(nuc_tfbs, sub_tfbs, on="sample_ID", how="inner")

    scatac_sites = _bed_site_names("scATAC_sites")
    dnase_sites = _bed_site_names("DHS_sites")

    def _sites_in(colname: str) -> str:
        # feature col -> base tissue name, e.g. `Placental_trophoblast_120_180_MCV`
        # -> `Placental_trophoblast`. Fragment tokens are always `120_180` or `35_80`.
        base = colname
        for frag in ("_120_180", "_35_80"):
            idx = base.find(frag)
            if idx != -1:
                return base[:idx]
        return base

    scatac_cols = ["sample_ID"] + [c for c in tissue.columns
                                   if c != "sample_ID" and _sites_in(c) in scatac_sites]
    dnase_cols = ["sample_ID"] + [c for c in tissue.columns
                                  if c != "sample_ID" and _sites_in(c) in dnase_sites]

    unmatched = [c for c in tissue.columns
                 if c != "sample_ID" and _sites_in(c) not in scatac_sites
                 and _sites_in(c) not in dnase_sites]
    if unmatched:
        print("WARNING: tissue columns not matched to scATAC or DHS bed dirs:")
        print(sorted(set(_sites_in(c) for c in unmatched)))

    scatac = tissue[scatac_cols]
    dnase = tissue[dnase_cols]
    return scatac, dnase, tfbs


def main():
    scatac, dnase, tfbs = build_tables()

    print("\nReconstructed table shapes (rows = samples, cols include sample_ID):")
    print("  Raw feature Table 1 (scATAC):", scatac.shape)
    print("  Raw feature Table 2 (DNase) :", dnase.shape)
    print("  Raw feature Table 3 (TFBS)  :", tfbs.shape)

    with pd.ExcelWriter(OUT_PATH, engine="openpyxl") as xw:
        scatac.to_excel(xw, sheet_name="Raw feature Table 1", index=False)
        dnase.to_excel(xw, sheet_name="Raw feature Table 2", index=False)
        tfbs.to_excel(xw, sheet_name="Raw feature Table 3", index=False)
    print("\nWrote", OUT_PATH)


if __name__ == "__main__":
    main()
