# Literature-informed BRD4-associated gene panel for exploratory validation.
# This is NOT a strict gold-standard set of direct BRD4 target genes.
# Genes are selected based on published BRD4/BET inhibitor sensitivity, BRD4 occupancy,
# BRD4-loaded enhancer/super-enhancer association, or BRD4-associated oncogenic transcriptional dependency.

brd4_core_panel <- c("MYC", "BCL2", "CDK6", "IL7R", "MYCN", "MYCL", "AR", "ERG", "CD274", "TP63", "MET", "FOSL1")

brd4_context_panel <- c("BCL6", "IRF8", "PAX5", "POU2AF1", "IRF4", "CRLF2", "ESR1", "TFF1", "GREB1", "PGR", "KLK3", "TMPRSS2", "FKBP5", "NKX3-1")

brd4_associated_panel <- unique(c(brd4_core_panel, brd4_context_panel))
known_targets <- brd4_associated_panel  # backward compatibility with existing scripts
