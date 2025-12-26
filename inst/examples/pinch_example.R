# Pinch Force Data Analysis with LFPCA
#
# This example demonstrates how to use the LFPCA package to analyze
# the pinch force data from the fda package.

# Load required packages
library(LFPCA)
library(fda)
library(refund)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# ============================================================
# Load and Smooth the "Pinch" Data
# ============================================================

data("pinch")
X        <- t(pinch)        # 301 time points × 24 curves
fp       <- fpca.face(X)    # default PVE cutoff
X.smooth <- fp$Yhat         # fitted curves

# ============================================================
# Decompose and Extract Contiguous Segments
# ============================================================

decomp           <- bdsvd::bdsvd(X.smooth, anp = "1")
numeric_segments <- convert_segments_to_numeric(decomp)
nSeg             <- length(numeric_segments)

# Build lookup: Time → Subprocess index as character
seg_df <- bind_rows(
  lapply(seq_along(numeric_segments), function(i) {
    data.frame(
      Time       = numeric_segments[[i]],
      Subprocess = as.character(i),
      stringsAsFactors = FALSE
    )
  })
)

# ============================================================
# Pivot into Long Form and Join Subprocess
# ============================================================

time_idx <- 1:151
long_df  <- as.data.frame(t(X.smooth)) %>%
  mutate(Time = time_idx) %>%
  pivot_longer(
    cols      = -Time,
    names_to  = "Curve",
    values_to = "Value"
  ) %>%
  left_join(seg_df, by = "Time") %>%
  filter(!is.na(Subprocess))  # keep only inside-segment points

# ============================================================
# Define Palettes and Labels
# ============================================================

pal_subp <- setNames(
  brewer.pal(max(3, nSeg), "Accent")[1:nSeg],
  as.character(seq_len(nSeg))
)

# ============================================================
# Run LFPCA and FPCA
# ============================================================

res.lfpca <- LFPCA(
  data  = X.smooth,
  Tmax  = length(time_idx),
  range = range(time_idx),
  Nsel  = 6,
  anp   = "1"
)

phi.df    <- res.lfpca$phi.df        # columns: t, phi, nr, Subprocess
lambdas   <- res.lfpca$eigenvalues   # length-6 eigenvalues
subprocs  <- res.lfpca$subprocess_numbers  # length-6 subprocess IDs

res.fpca <- fpca.face(X.smooth, npc = 6)
phi.fpca <- res.fpca$efunctions

# Number of PCs selected
N <- length(lambdas)

# Facet labels: % variance per PC
explvar_values <- lambdas / sum(lambdas)
custom_labels_pc <- setNames(
  paste0("PC", 1:N, ": Expl. Var: ", round(100 * explvar_values, 2), "%"),
  as.character(1:N)
)

# Legend labels: % variance per subprocess
pve_subp <- sapply(seq_len(nSeg), function(i) {
  sum(lambdas[subprocs == i]) / sum(lambdas)
})
custom_labels_subp <- setNames(
  paste0(seq_len(nSeg), ": PVE: ", round(100 * pve_subp, 2), "%"),
  as.character(seq_len(nSeg))
)

# ============================================================
# Custom Theme
# ============================================================

custom_theme <- theme(
  plot.title   = element_text(size = 24, face = "bold"),
  axis.title.x = element_text(size = 22),
  axis.title.y = element_text(size = 22),
  axis.text.x  = element_text(size = 20),
  axis.text.y  = element_text(size = 20),
  legend.text  = element_text(size = 20),
  legend.title = element_text(size = 22),
  strip.text   = element_text(size = 22)
)

# ============================================================
# Plot All Curves Colored by Subprocess
# ============================================================

p.curves <- ggplot(long_df,
                   aes(x = (Time - 1) / (151 - 1) * 0.30,
                       y = Value,
                       group = Curve,
                       colour = Subprocess)) +
  geom_line() +
  scale_colour_manual(
    name   = "Subprocess",
    values = pal_subp,
    labels = custom_labels_subp,
    guide  = guide_legend(nrow = 2)
  ) +
  scale_x_continuous(
    limits = c(0, 0.30),
    breaks = seq(0, 0.30, by = 0.05),
    labels = scales::number_format(accuracy = 0.05)
  ) +
  labs(
    x = "Time (s)",
    y = "Force (N)"
  ) +
  theme_minimal() +
  custom_theme +
  theme(legend.position = "none")

# ============================================================
# Sign-correct and Match FPCA Eigenfunctions
# ============================================================

# Scaled time grid
t.scaled <- seq(0, 0.3, length.out = length(time_idx))

# Build LFPCA matrix: rows=time, cols=1..N
lfpca_mat <- sapply(1:N, function(k) {
  phi.df %>%
    filter(nr == k) %>%
    arrange(t) %>%
    pull(phi)
})

# FPCA matrix (take first N columns)
fpca_mat <- phi.fpca[, 1:N]

# Correlation matrix
corr_mat <- cor(lfpca_mat, fpca_mat)

# Greedy one-to-one matching by max |corr|
match_idx      <- integer(N)
available_fpca <- 1:N
for (i in order(-apply(abs(corr_mat), 1, max))) {
  j_best <- available_fpca[which.max(abs(corr_mat[i, available_fpca]))]
  match_idx[i] <- j_best
  available_fpca <- setdiff(available_fpca, j_best)
}

# Flip sign of FPCA curves where corr < 0
for (i in seq_len(N)) {
  if (corr_mat[i, match_idx[i]] < 0) {
    fpca_mat[, match_idx[i]] <- -fpca_mat[, match_idx[i]]
  }
}

# Build a long data.frame for all N matched FPCA curves
phi_fpca_df <- data.frame(
  t.scaled = rep(t.scaled, times = N),
  phi      = as.vector(fpca_mat[, match_idx]),
  nr       = factor(rep(1:N, each = length(t.scaled)))
)

# ============================================================
# Plot First-3 LFPCA + Matched FPCA Overlays
# ============================================================

phi3 <- phi.df %>% filter(nr %in% 1:3)
phi3$Subprocess <- factor(phi3$Subprocess,
                          levels = as.character(seq_len(nSeg)))
phi3$t.scaled <- rep(t.scaled, 3)

p.eig3 <- ggplot() +
  # LFPCA curves
  geom_line(
    data   = phi3,
    aes(x = t.scaled, y = phi, colour = Subprocess, group = nr),
    linewidth = 1.2
  ) +
  # Matched FPCA curves, dotted black
  geom_line(
    data      = filter(phi_fpca_df, nr %in% factor(1:3)),
    aes(x = t.scaled, y = phi, group = nr),
    colour    = "black",
    linetype  = "dotted",
    linewidth = 0.8
  ) +
  facet_wrap(~nr, nrow = 1, labeller = labeller(nr = custom_labels_pc)) +
  theme_minimal() +
  scale_colour_manual(
    name   = "Subprocess",
    values = pal_subp,
    labels = custom_labels_subp,
    guide  = guide_legend(nrow = 2)
  ) +
  scale_x_continuous(
    limits = c(0, 0.3),
    breaks = seq(0, 0.30, by = 0.05),
    labels = scales::number_format(accuracy = 0.05)
  ) +
  labs(
    x = "Time (s)",
    y = expression(psi(t))
  ) +
  custom_theme +
  theme(legend.position = "bottom")

# ============================================================
# Combine and Display
# ============================================================

# To save as PDF:
# pdf("pinch_curves_and_eigen3.pdf", width = 16.5, height = 8)
# grid.arrange(p.curves, p.eig3, ncol = 1, heights = c(1, 1))
# dev.off()

# Display plots
print(p.curves)
print(p.eig3)
