
import pandas as pd
from scipy.stats import chisquare
from tqdm import tqdm

# --------------------------------------------------
# 1. Load simulated patient data
# --------------------------------------------------

FILE_PATH = "simulated_patients_formatted.jsonl"
OUTPUT_PATH = "simulated_age_chi_squared_results.csv"

df = pd.read_json(FILE_PATH, lines=True)

# Filter for diseases of interest
diseases_of_interest = [3318, 90695, 586]

df = df[
    df["disease_id"].isin(diseases_of_interest)
].copy()

# --------------------------------------------------
# 2. Convert age categories to numerical values
# --------------------------------------------------

mapping = {
    "Onset_Infant": 0,
    "Onset_Child": 1,
    "Onset_Adolescent": 2,
    "Onset_Adult": 3,
    "Onset_Elderly": 4
}

df["age_numeric"] = df["age"].map(mapping)

# Remove patients with missing or unrecognized ages
df = df.dropna(subset=["age_numeric"]).copy()

df["age_numeric"] = df["age_numeric"].astype(int)

# --------------------------------------------------
# 3. Perform chi-square goodness-of-fit tests
# --------------------------------------------------

results = []

age_categories = [0, 1, 2, 3, 4]

for orpha_code, group in tqdm(
    df.groupby("disease_id"),
    total=df["disease_id"].nunique(),
    desc="Performing chi-square tests"
):

    # Count patients in each age category
    observed = (
        group["age_numeric"]
        .value_counts()
        .reindex(age_categories, fill_value=0)
    )

    n_patients = len(group)

    # Expected counts assuming uniform distribution
    expected_per_category = n_patients / 5
    expected = [expected_per_category] * 5

    # Chi-square goodness-of-fit test
    chi_squared, p_value = chisquare(
        f_obs=observed.values,
        f_exp=expected
    )

    # Save results
    results.append({
        "ORPHA_CODE": orpha_code,
        "n_patients": n_patients,
        "chi_squared": chi_squared,
        "p_value": p_value,
        "significant_0.05": p_value < 0.05,
        "expected_per_category": expected_per_category,
        "n_infant": observed[0],
        "n_child": observed[1],
        "n_adolescent": observed[2],
        "n_adult": observed[3],
        "n_elderly": observed[4]
    })

# --------------------------------------------------
# 4. Save results to CSV
# --------------------------------------------------

results_df = pd.DataFrame(results)

results_df = results_df.sort_values(
    by="n_patients",
    ascending=False
)

results_df.to_csv(OUTPUT_PATH, index=False)

print(f"\nResults saved to {OUTPUT_PATH}")

# --------------------------------------------------
# 5. Display summary
# --------------------------------------------------

print("\nChi-square goodness-of-fit results:")
print(results_df.to_string(index=False))

print("\nSummary:")
print(f"Total diseases tested: {len(results_df)}")

print(
    "Diseases with p < 0.05:",
    (results_df["p_value"] < 0.05).sum()
)

print(
    "Diseases with p >= 0.05:",
    (results_df["p_value"] >= 0.05).sum()
)

# Check chi-square expected-count assumption
low_expected = results_df[
    results_df["expected_per_category"] < 5
]

if not low_expected.empty:
    print(
        "\nWarning:",
        len(low_expected),
        "diseases have expected counts below 5."
    )
    print(
        "Their chi-square p-values may be unreliable."
    )
