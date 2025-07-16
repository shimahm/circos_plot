# circos_plot
# Structural Variant Visualization in R

This repository contains an R script for visualizing structural variants (SVs) between two genomes using both linear and circos-style plots. The workflow is designed for comparative genomics, particularly for visualizing results from SyRI or similar SV callers.

## Features

- **Automated Data Cleaning:** Handles chromosome length and SyRI output files, parsing and validating coordinates.
- **Strict SV Filtering:** Only well-mapped, valid SVs (e.g., INV, TRANS, DUP, etc.) between known chromosomes are plotted.
- **Linear Plot:** Visualizes SVs as arcs and colored segments on a chromosome-by-chromosome basis.
- **Circos Plot:** Generates a circular ideogram with color-coded links for each variant type.
- **Customizable Colors and Layout:** Distinct colors for SV types, chromosome padding, and high-resolution outputs.

---

## Requirements

- **R version 4.0+**
- **R packages:**  
  - `ggplot2`
  - `circlize`
  - `dplyr`

Install missing packages with:
```r
install.packages(c("ggplot2", "dplyr"))
if (!requireNamespace("circlize", quietly = TRUE)) BiocManager::install("circlize")
```

---

## Input Files

- `Chifuu_V4.1_lenth.txt`  
  Tab-delimited file with chromosome names and lengths.  
  **Example:**
  ```
  chr1    150000000
  chr2    130000000
  ...
  ```

- `Exact_twotetraploid_pA7TA7_minimapsyri.out.txt`  
  Tab-delimited SyRI or SV caller output with at least 13 columns, including chromosome names, coordinates, and event type.

---

## How It Works

1. **Set Working Directory**  
   The script starts by setting the working directory to the location of your files.

2. **Data Loading & Cleaning**  
   - Reads chromosome lengths, adds a 10% padding, and standardizes chromosome names.
   - Parses SV data, standardizes coordinates, filters for allowed types, and ensures coordinates are valid.

3. **Plotting**
   - **Linear Plot:**  
     Each SV is visualized as a curve (arc) connecting its position on the reference and query chromosomes, with colored segments showing variant spans.
     Output: `linear_genome_plot_connected.png`
   - **Circos Plot:**  
     Chromosomes are arranged in a circle. Each SV is a colored link based on its type, with a legend for quick interpretation.
     Output: `circos_plot_connected.png`

---

## Usage

1. Save the script as `plot_structural_variants.R`.
2. Place it in the same directory as your data files.
3. Adjust the working directory in the script if needed.
4. Run the script in R:
    ```r
    source("plot_structural_variants.R")
    ```
5. Check the output PNG files in your working directory.

---

## Output

- `linear_genome_plot_connected.png`  
  High-resolution plot showing chromosomes and SV connections.

- `circos_plot_connected.png`  
  Circos-style plot with color-coded links for SVs.

---

## Customization

- **SV Types:**  
  Edit the `sv_colors` mapping in the script to add or change colors for event types.
- **Input Paths:**  
  Change the file names at the top of the script if your files are named differently.
- **Plot Appearance:**  
  Adjust parameters in `geom_segment`, `geom_curve`, or `circos.par` for different styles.

---

## License

MIT License

---

## Acknowledgments

- [SyRI](https://schneebergerlab.github.io/syri/) for structural variant detection.
- [`circlize`](https://jokergoo.github.io/circlize_book/book/) and [`ggplot2`](https://ggplot2.tidyverse.org/) for visualization tools.

---

## Example

![Linear Genome Plot Example](linear_genome_plot_connected.png)
![Circos Plot Example](circos_plot_connected.png)
