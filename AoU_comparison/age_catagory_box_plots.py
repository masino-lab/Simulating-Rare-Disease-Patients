# creating box plots of the disease age distribution for each disease

# read in the simulated_patients_formatted.jsonol file
# group by disease 
# map ages cats as numbers
# for each disease
#   convert age cat to number
#   make box plot of distribution
#   send output to another file

# in output file
#   have the orphacode at the top
#   then follow up with the box plot

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm
import io
import base64

# Age category mapping
mapping = {
    'Onset_Infant': 0,
    'Onset_Child': 1,
    'Onset_Adolescent': 2,
    'Onset_Adult': 3,
    'Onset_Elderly': 4
}

# Read simulated patients
df = pd.read_json("simulated_patients_formatted.jsonl", lines=True)

# Convert age categories to numerical values
df["age_numeric"] = df["age"].map(mapping)

# Remove patients with missing or unrecognized categories
df_plot = df.dropna(subset=["age_numeric"])

# Generate HTML with one box plot per disease
with open("disease_age_distributions.html", "w") as f:

    f.write("<html><body>\n")
    f.write("<h1>Simulated Disease Age Distributions</h1>\n")

    for orpha_code, group in tqdm(
        df_plot.groupby("disease_id"),
        total=df_plot["disease_id"].nunique(),
        desc="Generating box plots"
    ):

        fig, ax = plt.subplots(figsize=(8, 5))

        ax.boxplot(group["age_numeric"])

        ax.set_title(f"Orphacode: {orpha_code}")
        ax.set_ylabel("Onset Age Category")

        ax.set_yticks(range(5))
        ax.set_yticklabels([
            "Infant",
            "Child",
            "Adolescent",
            "Adult",
            "Elderly"
        ])

        ax.set_ylim(-0.5, 4.5)
        ax.set_xticks([])

        plt.tight_layout()

        # Convert plot to image in memory
        buffer = io.BytesIO()
        fig.savefig(buffer, format="png", dpi=100)
        buffer.seek(0)

        # Embed image into HTML
        image_base64 = base64.b64encode(buffer.read()).decode("utf-8")

        f.write(f'<img src="data:image/png;base64,{image_base64}">\n')
        f.write("<hr>\n")

        buffer.close()
        plt.close(fig)

    f.write("</body></html>")

print("Plots saved to disease_age_distributions.html")