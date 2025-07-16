# Load required libraries
library(ggplot2)
library(circlize)
library(dplyr)

# Set working directory
setwd("Y:/Bigdata/computing/Shima/circos")

# Helper to convert comma-quoted numbers to numeric
comma_to_numeric <- function(x) {
  as.numeric(gsub(",", "", gsub('"', "", x)))
}

# 1. Read and clean chromosome lengths
chr_len <- read.table("Chifuu_V4.1_lenth.txt", header = FALSE, 
                      col.names = c("chr", "length")) %>%
  mutate(length = comma_to_numeric(length),
         chr = toupper(chr),
         length = length * 1.10)  # Increased padding to 10%

# 2. Read and clean SyRI data with strict coordinate validation
syri_data <- read.table("Exact_twotetraploid_pA7TA7_minimapsyri.out.txt", 
                        header = FALSE, stringsAsFactors = FALSE, sep = "\t", quote = "") %>%
  transmute(
    chrA = toupper(V1),
    startA = comma_to_numeric(V2),
    endA = comma_to_numeric(V3),
    chrB = toupper(V7),
    startB = comma_to_numeric(V8),
    endB = comma_to_numeric(V9),
    type = V13
  ) %>%
  # Filter for valid SV types and chromosomes
  filter(type %in% c("INV", "TRANS", "INVTR", "DUP", "INVDP"),
         chrA %in% chr_len$chr, 
         chrB %in% chr_len$chr) %>%
  # Ensure start <= end coordinates
  mutate(
    startA = pmin(startA, endA),
    endA = pmax(startA, endA),
    startB = pmin(startB, endB),
    endB = pmax(startB, endB)
  ) %>%
  # Join with chromosome lengths and validate coordinates
  left_join(chr_len, by = c("chrA" = "chr")) %>%
  rename(lengthA = length) %>%
  left_join(chr_len, by = c("chrB" = "chr")) %>%
  rename(lengthB = length) %>%
  filter(
    startA >= 1, endA <= lengthA,  # Validate against chromosome lengths
    startB >= 1, endB <= lengthB
  ) %>%
  select(-lengthA, -lengthB) %>%
  # Factorize chromosomes for proper ordering
  mutate(
    chrA = factor(chrA, levels = chr_len$chr),
    chrB = factor(chrB, levels = chr_len$chr)
  )

message("✔ SyRI structural variants loaded: ", nrow(syri_data), " valid SVs")

# 3. Enhanced Linear Plot (shows only connected SVs)
png("linear_genome_plot_connected.png", width = 2000, height = 1200, res = 300)
ggplot() +
  # Chromosome backbones
  geom_segment(data = chr_len, 
               aes(x = 0, xend = length, y = chr, yend = chr), 
               linewidth = 3, color = "gray80") +
  # SV connections (curves) - drawn first
  geom_curve(data = syri_data,
             aes(x = startA, xend = startB, y = chrA, yend = chrB),
             curvature = 0.2, alpha = 0.8, linewidth = 0.8, color = "purple4") +
  # Reference genome segments (red)
  geom_segment(data = syri_data,
               aes(x = startA, xend = endA, y = chrA, yend = chrA),
               color = "red3", linewidth = 1.5) +
  # Query genome segments (blue)
  geom_segment(data = syri_data,
               aes(x = startB, xend = endB, y = chrB, yend = chrB),
               color = "blue3", linewidth = 1.5) +
  labs(title = "Structural Variants (Validated Connections)", 
       subtitle = "Red: Reference Genome | Blue: Query Genome | Purple: Links",
       x = "Position (bp)", 
       y = "Chromosome") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor = element_blank())
dev.off()

# 4. Enhanced Circos Plot with guaranteed connections
png("circos_plot_connected.png", width = 2000, height = 2000, res = 300)

# Initialize circos with tighter margins
circos.clear()
circos.par(
  "track.height" = 0.1,
  "gap.degree" = 3,
  "start.degree" = 90,
  "canvas.xlim" = c(-1.2, 1.2),
  "canvas.ylim" = c(-1.2, 1.2),
  "points.overflow.warning" = FALSE
)

# Initialize chromosomes
circos.initialize(
  factors = chr_len$chr,
  xlim = matrix(c(rep(0, nrow(chr_len)), chr_len$length), ncol = 2)
)

# Custom color palette for SV types
sv_colors <- c(
  "INV" = "#FF0000",    # Red
  "TRANS" = "#0000FF",  # Blue
  "DUP" = "#00AA00",    # Green
  "INVDP" = "#AA00AA",  # Purple
  "INVTR" = "#FFA500"   # Orange
)

# Draw chromosome tracks
circos.track(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    circos.rect(xlim[1], 0, xlim[2], 1, 
                col = rand_color(1, luminosity = "light"),
                border = NA)
    circos.text(mean(xlim), 1.2, chr, cex = 0.8, 
                facing = "downward", niceFacing = TRUE)
  },
  bg.border = NA
)

# Draw all valid links
for (i in 1:nrow(syri_data)) {
  sv_type <- syri_data$type[i]
  link_color <- sv_colors[sv_type]
  
  circos.link(
    sector.index1 = as.character(syri_data$chrA[i]),
    point1 = c(syri_data$startA[i], syri_data$endA[i]),
    sector.index2 = as.character(syri_data$chrB[i]),
    point2 = c(syri_data$startB[i], syri_data$endB[i]),
    col = paste0(link_color, "40"),  # 40% transparency
    border = link_color,             # Solid border matching SV type
    lwd = 3.5,                       # Thicker lines
    h.ratio = 0.5
  )
}

# Add legend
legend("bottomright", 
       legend = c("Inversion", "Translocation", "Duplication", "Inverted Dup", "Inverted Trans"),
       fill = c("#FF0000", "#0000FF", "#00AA00", "#AA00AA", "#FFA500"),
       border = "black",
       bty = "n",
       cex = 1.1,
       title = "Variant Type")

title("Circos Plot of Validated Structural Variants", cex.main = 1.5)
dev.off()

message("✅ Plots created successfully:")
message("• linear_genome_plot_connected.png")
message("• circos_plot_connected.png")