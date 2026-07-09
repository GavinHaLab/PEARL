"""
Reconstruct `Supplementary_Tables_v2.xlsx` (the `Meta_data` input of the notebook)
from the published Nature Medicine supplementary workbook.

The notebook loads metadata with:

    Meta_data = pd.read_excel("Supplementary_Tables_v2.xlsx",
                              sheet_name='Supplementary Table 2')

That internal "v2" workbook is missing.  The published supplementary file
`41591_2025_3509_MOESM2_ESM.xlsx` contains the same per-sample metadata, but:

  * it lives on sheet 'Supplementary Table 1' (not 'Supplementary Table 2'),
  * several columns were renamed (NIPS_* -> PDNAS_*, Griffin_* -> PEARL_*,
    Non-Pregnant_Sample -> Non-Pregnant_Female, ...),
  * an `Excluded` flag column was added (178 samples excluded from analysis),
  * cohort / subtype labels were renamed for publication.

Filtering to `Excluded == 0` and applying the label maps below reproduces the
cohort sizes reported in the paper exactly (PE-Training n=450 with
EPE=38 / LPE-PB=35 / LPE-TB=62 / NP=315; PDNAS validation n=831;
FF-training n=395), matching Extended Data Fig. 7.

This script writes a fresh `Supplementary_Tables_v2.xlsx` whose
'Supplementary Table 2' sheet has the column/label conventions the notebook
expects, so `PE_source_code.ipynb` runs unchanged.
"""

import os
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(ROOT, "41591_2025_3509_MOESM2_ESM.xlsx")
OUT = os.path.join(ROOT, "Supplementary_Tables_v2.xlsx")

# published column name -> name expected by the notebook
COLUMN_MAP = {
    "PDNAS_gestational_age": "NIPS_gestational_age",
    "PEARL_Fetal_Fraction": "Griffin_Fetal_Fraction",
    "Non-Pregnant_Female": "Non-Pregnant_Sample",
    "PDNAS_Maternal_age": "Maternal_age",
    "PDNAS_Maternal_bmi": "NIPS_Maternal_bmi",
    "PDNAS_sbp": "NIPS_sbp",
    "PDNAS_dbp": "NIPS_dbp",
    "Race": "race",
    "Ethnicity": "ethnicity",
}

# published cohort label -> notebook label
COHORT_MAP = {"PDNAS Validation Cohort": "Screening Cohort"}

# published subtype label -> notebook label
#   EPE    (early PE)                 -> EOPE
#   LPE    (late PE, term birth)      -> LOPE-TB
#   LPE-PB (late PE, preterm birth)   -> LOPE-PB
SUBTYPE_MAP = {"EPE": "EOPE", "LPE": "LOPE-TB", "LPE-PB": "LOPE-PB"}

# column order used by the notebook's Meta_data
COLUMN_ORDER = [
    "sample_ID", "Cohort", "Outcome", "Subtype", "NIPS_gestational_age", "Year",
    "mean_sequencing_coverage", "ChrY_Fetal_Fraction", "Griffin_Fetal_Fraction",
    "Neonatal_Sex", "Non-Pregnant_Sample", "Maternal_age", "NIPS_Maternal_bmi",
    "NIPS_sbp", "NIPS_dbp", "Onset_of_hypertension_Gestational_age",
    "Delivery_Gestational_age", "Birthweight", "Chronic_hypertension",
    "Pre_gestational_diabetes", "race", "ethnicity",
]


def build_meta_data() -> pd.DataFrame:
    df = pd.read_excel(SRC, sheet_name="Supplementary Table 1")

    df = df[df["Excluded"] == 0].copy()
    df = df.drop(columns=["Excluded"])

    df = df.rename(columns=COLUMN_MAP)
    df["Cohort"] = df["Cohort"].replace(COHORT_MAP)
    df["Subtype"] = df["Subtype"].replace(SUBTYPE_MAP)

    df = df[COLUMN_ORDER]
    return df


def main():
    meta = build_meta_data()

    print("Meta_data shape:", meta.shape)
    print("\nCohort counts:")
    print(meta["Cohort"].value_counts(dropna=False))
    print("\nPE-Training subtype counts:")
    print(meta[meta["Cohort"] == "PE-Training Cohort"]["Subtype"].value_counts(dropna=False))

    with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
        meta.to_excel(xw, sheet_name="Supplementary Table 2", index=False)
    print("\nWrote", OUT)


if __name__ == "__main__":
    main()
