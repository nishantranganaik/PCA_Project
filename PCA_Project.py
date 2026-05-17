

import matplotlib
matplotlib.use("Agg")  
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



expression = pd.read_csv("data/filtered.tsv.gz", sep="\t")


labels = pd.read_csv("data/class.csv", header=None)[0].values

print(f"Expression matrix shape: {expression.shape}")
print(f"Number of patients: {len(labels)}")
print(f"ER+ patients: {sum(labels == 1)}, ER- patients: {sum(labels == 0)}")


xbp1 = expression[" 4404"].values   # XBP1 expression
gata3 = expression[" 4359"].values  # GATA3 expression

print(f"\nXBP1 range: [{xbp1.min():.2f}, {xbp1.max():.2f}]")
print(f"GATA3 range: [{gata3.min():.2f}, {gata3.max():.2f}]")



fig1, ax1 = plt.subplots(figsize=(6, 5))


er_pos = labels == 1
er_neg = labels == 0

ax1.scatter(gata3[er_neg], xbp1[er_neg], c="black", s=30, marker="s", label="ER-")
ax1.scatter(gata3[er_pos], xbp1[er_pos], c="red", s=30, marker="s", label="ER+")

ax1.set_xlabel("GATA3", fontsize=13, fontstyle="italic")
ax1.set_ylabel("XBP1", fontsize=13, fontstyle="italic")
ax1.legend(loc="lower right", fontsize=10)


ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig("figure_1a.png", dpi=150, bbox_inches="tight")

data_2d = np.column_stack([gata3, xbp1])  # shape: (105, 2)


mean_vec = data_2d.mean(axis=0)
centered = data_2d - mean_vec


cov_matrix = np.cov(centered, rowvar=False)
print(f"\nCovariance matrix:\n{cov_matrix}")

eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)


sort_idx = np.argsort(eigenvalues)[::-1]
eigenvalues = eigenvalues[sort_idx]
eigenvectors = eigenvectors[:, sort_idx]

pc1_direction = eigenvectors[:, 0]  # first principal component
print(f"PC1 direction: {pc1_direction}")
print(f"Variance explained by PC1: {eigenvalues[0]/eigenvalues.sum()*100:.1f}%")

# Project all points onto PC1
projections = centered @ pc1_direction  # dot product with PC1

# Flip sign if needed so that ER+ patients have higher PC1 values (matching paper)
if projections[er_pos].mean() < projections[er_neg].mean():
    projections = -projections



fig1b, ax1b = plt.subplots(figsize=(6, 5))

# Same scatter plot as Figure 1a
ax1b.scatter(gata3[er_neg], xbp1[er_neg], c="black", s=30, marker="s", label="ER-")
ax1b.scatter(gata3[er_pos], xbp1[er_pos], c="red", s=30, marker="s", label="ER+")

# Get the PC directions (in GATA3, XBP1 space)
pc1_vec = eigenvectors[:, 0]  # PC1 direction
pc2_vec = eigenvectors[:, 1]  # PC2 direction


line_length = 5  

# PC1 line
pc1_start = mean_vec - line_length * pc1_vec
pc1_end = mean_vec + line_length * pc1_vec
ax1b.plot([pc1_start[0], pc1_end[0]], [pc1_start[1], pc1_end[1]],
          "k-", linewidth=1.2)

# PC2 line
pc2_start = mean_vec - line_length * pc2_vec
pc2_end = mean_vec + line_length * pc2_vec
ax1b.plot([pc2_start[0], pc2_end[0]], [pc2_start[1], pc2_end[1]],
          "k-", linewidth=1.2)


if pc1_vec[0] < 0:
    pc1_vec = -pc1_vec  # flip so it points toward upper-right
if pc2_vec[0] > 0 and pc2_vec[1] < 0:
    pc2_vec = -pc2_vec  # flip so it points toward upper-left

arrow_len = 2.0
# PC1 arrow
ax1b.annotate("", xy=(mean_vec[0] + arrow_len * pc1_vec[0],
                       mean_vec[1] + arrow_len * pc1_vec[1]),
              xytext=(mean_vec[0], mean_vec[1]),
              arrowprops=dict(arrowstyle="->", lw=2))
ax1b.text(mean_vec[0] + (arrow_len + 0.2) * pc1_vec[0],
          mean_vec[1] + (arrow_len + 0.2) * pc1_vec[1],
          "PC1", fontsize=12, fontweight="bold")

# PC2 arrow
ax1b.annotate("", xy=(mean_vec[0] + arrow_len * pc2_vec[0],
                       mean_vec[1] + arrow_len * pc2_vec[1]),
              xytext=(mean_vec[0], mean_vec[1]),
              arrowprops=dict(arrowstyle="->", lw=2))
ax1b.text(mean_vec[0] + (arrow_len + 0.2) * pc2_vec[0] - 0.5,
          mean_vec[1] + (arrow_len + 0.2) * pc2_vec[1],
          "PC2", fontsize=12, fontweight="bold")

ax1b.set_xlabel("GATA3", fontsize=13, fontstyle="italic")
ax1b.set_ylabel("XBP1", fontsize=13, fontstyle="italic")
ax1b.spines["top"].set_visible(False)
ax1b.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig("figure_1b.png", dpi=150, bbox_inches="tight")


fig2, ax2 = plt.subplots(figsize=(7, 4))

# We need three horizontal strips: "All", "ER-", "ER+"
y_positions = {"ER+": 1, "ER-": 2, "All": 3}

# Plot "All" row — both classes together
ax2.scatter(projections[er_neg], np.full(er_neg.sum(), y_positions["All"]),
            c="black", s=20, marker="s")
ax2.scatter(projections[er_pos], np.full(er_pos.sum(), y_positions["All"]),
            c="red", s=20, marker="s")

# Plot "ER-" row
ax2.scatter(projections[er_neg], np.full(er_neg.sum(), y_positions["ER-"]),
            c="black", s=20, marker="s")
ax2.scatter(projections[er_pos], np.full(er_pos.sum(), y_positions["ER-"]),
            c="red", s=20, marker="s")

# Plot "ER+" row
ax2.scatter(projections[er_pos], np.full(er_pos.sum(), y_positions["ER+"]),
            c="red", s=20, marker="s")
ax2.scatter(projections[er_neg], np.full(er_neg.sum(), y_positions["ER+"]),
            c="black", s=20, marker="s")

ax2.axhline(y=1.5, color="gray", linewidth=0.5)
ax2.axhline(y=2.5, color="gray", linewidth=0.5)


ax2.set_yticks([1, 2, 3])
ax2.set_yticklabels(["ER+", "ER−", "All"], fontsize=12)
ax2.set_xlabel("Projection onto PC1", fontsize=12)


ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.spines["left"].set_visible(False)
ax2.tick_params(left=False)
ax2.set_ylim(0.4, 3.6)

plt.tight_layout()
plt.savefig("figure_1c.png", dpi=150, bbox_inches="tight")