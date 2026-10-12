# Container for the R steps. Build and push once, then set process.container in nextflow.config.
#   docker build -t ghcr.io/OWNER/scrna-nf-seurat:0.1.0 .
#   docker push ghcr.io/OWNER/scrna-nf-seurat:0.1.0
# NOTE: not yet built or tested; expect to adjust system libraries on the first build.

FROM rocker/r-ver:4.4.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libhdf5-dev libglpk-dev \
    libpng-dev libjpeg-dev libtiff5-dev libcairo2-dev libfontconfig1-dev \
    libfreetype6-dev libharfbuzz-dev libfribidi-dev procps \
    && rm -rf /var/lib/apt/lists/*

# rocker/r-ver pins CRAN to a dated snapshot, so package versions are reproducible
RUN R -e 'install.packages(c("BiocManager", "optparse", "Seurat", "harmony", "ggplot2"))'

RUN R -e 'BiocManager::install(c("scDblFinder", "SingleR", "celldex", "DESeq2"), ask = FALSE, update = FALSE)'
