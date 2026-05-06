import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
import sys


def plot(final_tsv, color_contigs_tsv, reference_colors_tsv, plot_png, title):
    # Read table data
    df = pd.read_csv(final_tsv, sep="\t", index_col=0)

    feature_names = df.columns

    # PCA down to 2-D
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(df.values)
    df_pca = pd.DataFrame(X_pca, index=df.index, columns=["PC1", "PC2"])

    # Filter out outlier points in plot
    center = df_pca.mean().values.reshape(1, -1)
    distances = cdist(df_pca, center)
    threshold = np.percentile(distances, 95)  # keep 95% closest points
    df_pca = df_pca[distances.flatten() < threshold]

    # Colors
    color_df = pd.read_csv(color_contigs_tsv, sep="\t")
    color_map = dict(zip(color_df.iloc[:, 0], color_df.iloc[:, 1]))

    point_colors = df_pca.index.map(color_map).fillna("black").tolist()

    # Loadings
    loadings = pca.components_.T

    # Create plot
    fig, ax = plt.subplots(figsize=(14, 10), dpi=150)
    ax.scatter(
        df_pca["PC1"],
        df_pca["PC2"],
        c=point_colors,
        alpha=0.7,
        s=30,
        edgecolor="k",
        linewidth=0.2,
    )

    # Feature arrows
    scale = 0.75
    for i, feature in enumerate(feature_names):
        x, y = loadings[i] * scale
        ax.arrow(
            0,
            0,
            x,
            y,
            color="gray",
            alpha=0.5,
            head_width=0.03,
            length_includes_head=True,
        )
        ax.text(
            x * 1.1,
            y * 1.1,
            feature,
            fontsize=8,
        )

    # Axis labels
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")

    # Legend
    ref_df = pd.read_csv(reference_colors_tsv, sep="\t")
    legend_handles = [
        mpatches.Patch(color=row.color, label=row.species)
        for row in ref_df.itertuples(index=False)
    ]
    ax.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc="upper left")

    # Origin lines
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)

    # Plot styling
    ax.set_title(title)
    plt.tight_layout(rect=[0, 0, 0.8, 1])
    plt.savefig(plot_png, bbox_inches="tight", dpi=300)


if __name__ == "__main__":
    final_tsv = sys.argv[1]
    color_contigs_tsv = sys.argv[2]
    reference_colors_tsv = sys.argv[3]
    plot_png = sys.argv[4]
    title = sys.argv[5].upper().replace("_", " ")

    plot(final_tsv, color_contigs_tsv, reference_colors_tsv, plot_png, title)
