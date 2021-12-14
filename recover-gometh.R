# recover-gometh.R
# by Umair Khan

# Try to recover the performance of GOmeth in the B-cell vs. NK contrast.

#########################
# 1: Basic housekeeping #
#########################

# Load libraries
library(minfi)
library(limma)
library(DMRcate)
library(missMethyl)
library(dplyr)
library(tibble)
library(ggplot2)

# Define dictionary for figures
legendLabels <- c("goregion-gometh" = "GOregion",
                  "gometh-probe-top" = "GOmeth (5000)",
                  "gometh-probe-fdr-0.05" = "GOmeth (FDR < 0.05)",
                  "gometh-probe-fdr-0.04" = "GOmeth (FDR < 0.04)",
                  "gometh-probe-fdr-0.03" = "GOmeth (FDR < 0.03)",
                  "gometh-probe-fdr-0.02" = "GOmeth (FDR < 0.02)",
                  "gometh-probe-fdr-0.01" = "GOmeth (FDR < 0.01)",
                  "gometh-probe-fdr-0.001" = "GOmeth (FDR < 0.001)",
                  "gometh-probe-fdr-1e-04" = "GOmeth (FDR < 0.0001)")

# Read in pre-computed data files
load("data/GSE110554-data.RData")
load("data/annEPIC.RData")
dmrList <- readRDS("data/dmrcate-results.rds")
dmrGO <- readRDS("data/dmrcate-go.rds")

# Rename existing p < 0.05 category
dmrGO$method[dmrGO$method == "gometh-probe-fdr"] <- "gometh-probe-fdr-0.05"

# Read in truth set
immuneGO <- unique(read.csv("data/GO-immune-system-process.txt",
                            stringsAsFactors = FALSE, header = FALSE,
                            col.names = "GOID"))

# Perform statistical analysis (same as paper)
mVals <- getM(fltGr)
bVals <- getBeta(fltGr)
design <- model.matrix(~0+targets$CellType)
colnames(design) <- levels(factor(targets$CellType))
fit <- lmFit(mVals, design)
cont.matrix <- makeContrasts(CD4vCD8 = CD4T-CD8T, MonovNeu = Mono-Neu,
                             BcellvNK = Bcell-NK, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
tfit <- eBayes(fit2, robust = TRUE, trend = TRUE)
tfit <- treat(tfit, lfc = 0.5)

#########################
# 2: try different FDRs #
#########################

# Specify p-values to run through
pvals <- c(0.04, 0.03, 0.02, 0.01, 0.001, 0.0001)

# Loop over p-values
for (pval in pvals) {

    # Loop over dmrList
    for(i in 1:length(dmrList)){

        # Perform GSEA
        tmp <- topGSA(gometh(rownames(topTreat(tfit, coef = i, num = Inf, p.value = pval)),
                             anno = ann, array.type = "EPIC"), number = Inf)
        tmp <- rownames_to_column(tmp, var = "GO")[, c("GO", "P.DE")]
        tmp$method <- paste("gometh-probe-fdr-", pval, sep = "")
        tmp$contrast <- colnames(cont.matrix)[i]

        # Add to data frame
        dmrGO <- bind_rows(dmrGO, tmp)

        # Print progress
        print(paste("Finished running contrast ", i, " with FDR ", pval, ".", sep = ""))

    }

}

# Data frame manipulation
dmrGO %>% dplyr::filter(method %in% c("goregion-gometh", "gometh-probe-top",
                                      "gometh-probe-fdr-0.05",
                                      "gometh-probe-fdr-0.04",
                                      "gometh-probe-fdr-0.03",
                                      "gometh-probe-fdr-0.02",
                                      "gometh-probe-fdr-0.01",
                                      "gometh-probe-fdr-0.001",
                                      "gometh-probe-fdr-1e-04")) %>%
          dplyr::filter(contrast %in% c("BcellvNK")) %>%
          mutate(method = unname((legendLabels[method]))) %>%
          arrange(contrast, method, P.DE) %>%
          group_by(contrast, method) %>%
          mutate(csum = cumsum(GO %in% immuneGO$GOID)) %>%
          mutate(rank = 1:n()) %>%
          dplyr::filter(rank <= 100) -> sub

# Plot new curves
p1 <- ggplot(sub, aes(x = rank, y = csum, colour = method)) +
      geom_line() +
      facet_wrap(vars(contrast), ncol = 9) +
      geom_vline(xintercept = 10, linetype = "dotted") +
      labs(colour = "Method", x = "Rank", y = "Cumulative no. in truth set") +
      theme(legend.position = "bottom")

##########################
# 3: random CpG sampling #
##########################

# Define numbers of CpGs to use, FDR thresholds, and replicates
cpg_nums <- c(40000, 35000, 30000, 25000, 20000, 15000, 10000, 5000)
pvals <- c(0.0001, 0.05)
repls <- 1:10

# Loop over CpG numbers
for (cpg_num in cpg_nums) {

    # For each CpG number, loop through thresholds
    for (pval in pvals) {

        # Get list of CpGs at threshold (only B-cell vs NK contrast)
        cpg_list <- rownames(topTreat(tfit, coef = 3, num = Inf, p.value = pval))

        # For each threshold, loop through replicates
        for (repl in repls) {

            # Make run name and sample CpGs to use
            run_name <- paste("gometh-probe-fdr", pval, cpg_num, repl, sep = "-")
            cpg_sample <- sample(cpg_list, cpg_num)

            # Run through GOmeth
            tmp <- topGSA(gometh(cpg_sample, anno = ann, array.type = "EPIC"),
                          number = Inf)
            tmp <- rownames_to_column(tmp, var = "GO")[, c("GO", "P.DE")]
            tmp$method <- run_name
            tmp$contrast <- colnames(cont.matrix)[3]

            # Add to overall data frame
            dmrGO <- bind_rows(dmrGO, tmp)

            # Print progress
            print(paste("Finished running replicate", repl, "for FDR", pval,
                        "with", cpg_num, "CpGs.", sep = " "))

        }

    }

}

# Define new data frame to store results
sampling_data <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(sampling_data) <- c("cpgs", "method", "rep", "csum10", "csum100")
sampling_data$method <- as.character(sampling_data$method)

# Same loops as before, extracting sums this time
# (this could have all been done at once, but separating it makes
# it easier to use a pre-computed data file as input)
for (cpg_num in cpg_nums) {
    for (pval in pvals) {
        for (repl in repls) {

            # Make run name and method name
            run_name <- paste("gometh-probe-fdr", pval, cpg_num, repl, sep = "-")
            method_name <- paste("gometh-probe-fdr", pval, sep = "-")

            # Rank terms and compare to truth set
            dmrGO %>% dplyr::filter(method %in% c(run_name)) %>%
                      dplyr::filter(contrast %in% c("BcellvNK")) %>%
                      arrange(contrast, method, P.DE) %>%
                      group_by(contrast, method) %>%
                      mutate(csum = cumsum(GO %in% immuneGO$GOID)) %>%
                      mutate(rank = 1:n()) %>%
                      dplyr::filter(rank <= 100) -> ranked

            # Get the cumulative sums
            sum10 <- ranked$csum[ranked$rank == 10][1]
            sum100 <- ranked$csum[ranked$rank == 100][1]

            # Add to new data frame
            sampling_data %>% add_row(cpgs = cpg_num, method = method_name, rep = repl,
                                      csum10 = sum10, csum100 = sum100) -> sampling_data

            # Print progress
            print(paste("Extracted sums for replicate", repl, "for FDR", pval,
                        "with", cpg_num, "CpGs.", sep = " "))

        }
    }
}

# Specify order of x-axis since we won't be using it as a number
x_order <- as.character(cpg_nums)

# We are comparing to GOregion and GOmeth (5000)
csum100_targets <- data.frame(Targets = c("GOregion", "GOmeth (5000)"), target = c(48, 43))
csum10_targets <- data.frame(Targets = c("GOregion", "GOmeth (5000)"), target = c(8, 9))

# Rename methods using dictionary
sampling_data %>% mutate(method = unname((legendLabels[method]))) -> sampling_data

# Plot results for top 100 terms
p2 <- ggplot(sampling_data, aes(x = factor(cpgs, level = x_order), y = csum100, fill = method)) +
      geom_boxplot() +
      geom_hline(data = csum100_targets, aes(yintercept = target, linetype = Targets)) +
      labs(fill = "Dataset", x = "num. of CpGs sampled",
           y = "num. of top 100 terms in truth set") +
      theme(legend.position = "bottom")

# Plot results for top 10 terms
p3 <- ggplot(sampling_data, aes(x = factor(cpgs, level = x_order), y = csum10, fill = method)) +
      geom_boxplot() +
      geom_hline(data = csum10_targets, aes(yintercept = target, linetype = Targets)) +
      labs(fill = "Dataset", x = "num. of CpGs sampled",
           y = "num. of top 10 terms in truth set") +
      theme(legend.position = "bottom")
