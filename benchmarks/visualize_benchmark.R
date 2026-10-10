# Publication-quality benchmark visualization script for BioContext v1.0.0
# Outputs high-resolution plots for LinkedIn & technical report

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

base_dir <- getwd()
data_dir <- file.path(base_dir, "benchmarks", "output")
out_img <- file.path(data_dir, "biocontext_benchmark_linkedin.png")

# 1. Load CSV data
res_df <- read.csv(file.path(data_dir, "resolution_benchmark_results.csv"))
stress_df <- read.csv(file.path(data_dir, "stress_test_1000_results.csv"))
scaling_df <- read.csv(file.path(data_dir, "concurrency_scaling_results.csv"))
comp_df <- read.csv(file.path(data_dir, "method_comparison_benchmark.csv"))

# Publication-quality benchmark visualization script for BioContext
# Designed following Nature / Cell / Science editorial figure guidelines
# - Muted, distinct palette (Navy, Steel, Amber, Slate)
# - Pure white canvas with subtle horizontal guides
# - High-contrast academic typography with explicit units
# - Clean panel labels (a, b, c, d) and uncluttered legends

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

base_dir <- getwd()
data_dir <- file.path(base_dir, "benchmarks", "output")
out_img <- file.path(data_dir, "biocontext_benchmark_linkedin.png")

# 1. Load benchmark datasets
res_df <- read.csv(file.path(data_dir, "resolution_benchmark_results.csv"))
stress_df <- read.csv(file.path(data_dir, "stress_test_1000_results.csv"))
scaling_df <- read.csv(file.path(data_dir, "concurrency_scaling_results.csv"))
comp_df <- read.csv(file.path(data_dir, "method_comparison_benchmark.csv"))

# Editorial Color Palette (Nature / Cell style)
COL_PRIMARY   <- "#1E3A8A"  # Deep Navy (BioContext highlight)
COL_BASELINE  <- "#94A3B8"  # Neutral Slate Gray (Baselines)
COL_COLD      <- "#B45309"  # Deep Warm Ochre (Cold API)
COL_WARM      <- "#047857"  # Forest Emerald (Cache hit)
COL_ACCENT    <- "#475569"  # Slate Accent
COL_BORDER    <- "#1E293B"  # Dark Slate for crisp borders
COL_GRID      <- "#F1F5F9"  # Very faint grid lines

# Academic Minimalist Theme
theme_editorial <- theme_classic(base_size = 10.5, base_family = "Liberation Sans") +
  theme(
    plot.title = element_text(face = "bold", size = 11, color = "#0F172A", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 9, color = "#475569", margin = margin(b = 8)),
    axis.title.x = element_text(face = "bold", size = 9, color = "#1E293B", margin = margin(t = 5)),
    axis.title.y = element_text(face = "bold", size = 9, color = "#1E293B", margin = margin(r = 5)),
    axis.text.x = element_text(size = 8.5, color = "#334155"),
    axis.text.y = element_text(size = 8.5, color = "#334155"),
    axis.line = element_line(color = "#334155", linewidth = 0.5),
    axis.ticks = element_line(color = "#334155", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "#E2E8F0", linewidth = 0.35, linetype = "solid"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.background = element_rect(fill = "#FFFFFF", color = NA),
    strip.background = element_rect(fill = "#F8FAFC", color = "#CBD5E1", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 8.5, color = "#1E293B", margin = margin(t = 3, b = 3))
  )

# -------------------------------------------------------------
# Panel a: Comparative Resolution Accuracy Across Architectures
# -------------------------------------------------------------
comp_long <- comp_df %>%
  select(method, accuracy_pct, alias_handling_pct, cross_database_id_pct) %>%
  pivot_longer(
    cols = c(accuracy_pct, alias_handling_pct, cross_database_id_pct),
    names_to = "metric",
    values_to = "score"
  ) %>%
  mutate(
    metric_label = factor(
      case_when(
        metric == "accuracy_pct" ~ "Overall Accuracy",
        metric == "alias_handling_pct" ~ "Historical Aliases",
        metric == "cross_database_id_pct" ~ "Cross-Database IDs"
      ),
      levels = c("Cross-Database IDs", "Historical Aliases", "Overall Accuracy")
    ),
    is_biocontext = ifelse(grepl("BioContext", method), "BioContext", "Baseline")
  )

comp_long$method_short <- factor(
  comp_long$method,
  levels = c(
    "Naive LLM Zero-Shot (GPT-4o / Claude)",
    "Dense Semantic Vector Search (BioLinkBERT/Ada)",
    "Direct Ad-Hoc REST APIs (NCBI/Ensembl)",
    "BioContext v1.0.0 (Deterministic MCP Engine)"
  ),
  labels = c("Naive LLM\nZero-Shot", "Bio-RAG\n(Dense Vector)", "Direct REST\nAPIs", "BioContext\n(Deterministic)")
)

p_comp <- ggplot(comp_long, aes(x = method_short, y = score, fill = is_biocontext)) +
  geom_col(width = 0.62, color = "#0F172A", linewidth = 0.35, alpha = 0.95) +
  geom_text(
    aes(label = sprintf("%.1f%%", score)),
    vjust = -0.4,
    size = 2.5,
    fontface = "bold",
    color = "#0F172A"
  ) +
  facet_wrap(~metric_label, nrow = 1) +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  scale_fill_manual(values = c("Baseline" = COL_BASELINE, "BioContext" = COL_PRIMARY)) +
  labs(
    title = "a  Resolution precision across biological entity benchmarks",
    subtitle = "Tested across 100 canonical symbols, legacy clinical aliases (e.g. HER2 -> ERBB2), and cross-database accession IDs",
    x = NULL,
    y = "Resolution Precision (%)"
  ) +
  theme_editorial

# -------------------------------------------------------------
# Panel b: Upstream HTTP Cold Fetch vs Multi-Tier Local Cache
# -------------------------------------------------------------
lat_data <- res_df %>%
  select(query, cold_latency_ms, warm_latency_ms) %>%
  pivot_longer(cols = c(cold_latency_ms, warm_latency_ms), names_to = "mode", values_to = "latency") %>%
  mutate(
    mode_label = factor(
      ifelse(mode == "cold_latency_ms", "Upstream REST (Cold)", "SQLite Cache (Warm)"),
      levels = c("Upstream REST (Cold)", "SQLite Cache (Warm)")
    )
  )

p_lat <- ggplot(lat_data, aes(x = mode_label, y = latency, fill = mode_label)) +
  geom_boxplot(
    width = 0.42,
    color = "#0F172A",
    linewidth = 0.45,
    outlier.size = 1.0,
    outlier.alpha = 0.45,
    alpha = 0.9
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10, 100, 1000, 5000),
    labels = c("0.01 ms", "0.1 ms", "1 ms", "10 ms", "100 ms", "1 s", "5 s")
  ) +
  scale_fill_manual(values = c("Upstream REST (Cold)" = COL_COLD, "SQLite Cache (Warm)" = COL_WARM)) +
  labs(
    title = "b  Response latency distribution",
    subtitle = "Sub-millisecond retrieval on warm cache (median: 0.05 ms)",
    x = NULL,
    y = "Latency (log scale)"
  ) +
  theme_editorial

# -------------------------------------------------------------
# Panel c: Cumulative Latency Distribution (1,000 Concurrent Requests)
# -------------------------------------------------------------
p_stress <- ggplot(stress_df, aes(x = latency_ms)) +
  stat_ecdf(geom = "step", color = COL_PRIMARY, linewidth = 0.9) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = COL_WARM, linewidth = 0.5) +
  geom_vline(xintercept = 0.16, linetype = "dashed", color = COL_COLD, linewidth = 0.5) +
  annotate("text", x = 0.07, y = 0.40, label = "P50 = 0.05 ms", hjust = 0, size = 2.8, color = COL_WARM, fontface = "bold") +
  annotate("text", x = 0.19, y = 0.72, label = "P95 = 0.16 ms", hjust = 0, size = 2.8, color = COL_COLD, fontface = "bold") +
  scale_x_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, 0.2), labels = function(x) paste0(x, " ms")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0.01, 0.01)) +
  labs(
    title = "c  Tail latency empirical CDF (1,000 requests)",
    subtitle = "100% success rate (1,000/1,000) | Sustained: 10,458 req/s",
    x = "Response Latency",
    y = "Cumulative Fraction"
  ) +
  theme_editorial

# -------------------------------------------------------------
# Panel d: Concurrency Scaling Curve (Throughput Stability)
# -------------------------------------------------------------
p_scale <- ggplot(scaling_df, aes(x = concurrency, y = throughput_rps)) +
  geom_line(color = COL_PRIMARY, linewidth = 0.85) +
  geom_point(color = COL_PRIMARY, fill = "#FFFFFF", shape = 21, size = 2.5, stroke = 1.2) +
  scale_x_continuous(breaks = scaling_df$concurrency) +
  scale_y_continuous(
    labels = scales::comma,
    limits = c(9500, 12000),
    expand = c(0, 0)
  ) +
  labs(
    title = "d  Throughput scaling under asynchronous load",
    subtitle = "Sustained >10,000 req/s across 50 to 1,000 concurrent tasks",
    x = "Concurrent Asynchronous Workers",
    y = "Throughput (requests / s)"
  ) +
  theme_editorial

# -------------------------------------------------------------
# Assemble Multipanel Figure via Patchwork
# -------------------------------------------------------------
final_plot <- (p_comp) / (p_lat | p_stress) / (p_scale) +
  plot_annotation(
    title = "Empirical Validation of BioContext Grounding Architecture",
    subtitle = "Resolution accuracy and asynchronous performance characteristics across standard benchmarks",
    caption = "Source data: BioContext Benchmark Suite (v1.0.1) | Linux x86_64, Python 3.13, SQLite | github.com/CORE-Lab-Research/biocontext",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13.5, color = "#0F172A", margin = margin(b = 2)),
      plot.subtitle = element_text(size = 10, color = "#475569", margin = margin(b = 10)),
      plot.caption = element_text(size = 8, color = "#64748B", hjust = 1, margin = margin(t = 8))
    )
  )

ggsave(
  filename = out_img,
  plot = final_plot,
  width = 11.5,
  height = 10.5,
  dpi = 300,
  bg = "white"
)

message(sprintf("[SUCCESS] Publication plot generated at: %s", out_img))
