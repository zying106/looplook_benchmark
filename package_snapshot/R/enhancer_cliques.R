# Functions for clustering enhancer-like anchors into spatial
# enhancer cliques (E-E family) or promoter hubs (P-P family).
#
# The module consumes the output of annotate_peaks_and_loops() and,
# optionally, its expression- and chromatin-aware refinement results.
# It does NOT re-run the refinements: reclassified anchor types and
# chromatin evidence grades are inherited from the supplied object.

utils::globalVariables(c(
  "clique_id", "community_id", "effective_type", "promotion_reason",
  "weakest_evidence", "n_chromatin_validated", "is_trans_chromosomal",
  "n_internal_loops", "n_genes", "anchor_genes",
  "linked_promoter_genes", "gene_role", "n_paths", "n_anchors",
  "chromatin_evidence", "anchor_gene", "anchor_type", "chr", "start", "end",
  "prom_ok", "ev_ok", "clique", "community", "drop_category", "dropped",
  "intra_clique_degree", "mean_intra_degree", "max_intra_degree"
))

#' Internal: Build the clique anchor map
#'
#' Derives one row per anchor ID. Coordinates come from the canonical
#' anchor registry when available (merged-anchor hulls, matching
#' \code{validate_epeG_by_chromatin()} semantics); anchor types and genes
#' are overridden by the loop-table values so that expression/chromatin
#' refinements are inherited. Anchors missing from the registry are
#' recovered from the loop rows.
#' @keywords internal
#' @noRd
.clique_anchor_map <- function(annotation_res, sel) {
  reg <- tryCatch(.get_anchor_registry(annotation_res), error = function(e) NULL)
  if (!is.null(reg)) {
    m <- S4Vectors::mcols(reg)
    am <- data.frame(
      anchor_id = as.character(m$anchor_id),
      chr = as.character(GenomicRanges::seqnames(reg)),
      start = as.integer(GenomicRanges::start(reg)),
      end = as.integer(GenomicRanges::end(reg)),
      anchor_type = as.character(m$type_code),
      anchor_gene = as.character(m$SYMBOL),
      stringsAsFactors = FALSE
    )
  } else {
    am <- data.frame(
      anchor_id = character(), chr = character(),
      start = integer(), end = integer(), anchor_type = character(),
      anchor_gene = character(), stringsAsFactors = FALSE
    )
  }
  needed <- unique(c(as.character(sel$a1_id), as.character(sel$a2_id)))
  needed <- needed[!is.na(needed) & needed != ""]
  missing_ids <- setdiff(needed, am$anchor_id)
  if (length(missing_ids) > 0) {
    side1 <- sel[sel$a1_id %in% missing_ids, c("a1_id", "chr1", "start1", "end1", "anchor1_type", "anchor1_gene")]
    side2 <- sel[sel$a2_id %in% missing_ids, c("a2_id", "chr2", "start2", "end2", "anchor2_type", "anchor2_gene")]
    colnames(side1) <- c("anchor_id", "chr", "start", "end", "anchor_type", "anchor_gene")
    colnames(side2) <- c("anchor_id", "chr", "start", "end", "anchor_type", "anchor_gene")
    sides <- rbind(side1, side2)
    sides <- unique(sides)
    sides <- sides[!is.na(sides$anchor_id), , drop = FALSE]
    am <- rbind(am, sides)
  }
  am <- am[am$anchor_id %in% needed, , drop = FALSE]
  am <- am[!duplicated(am$anchor_id), , drop = FALSE]

  # Authoritative types/genes come from the loop table so that
  # expression (P -> eP) and chromatin reclassifications are inherited.
  side_types <- rbind(
    data.frame(
      anchor_id = as.character(sel$a1_id),
      anchor_type = as.character(sel$anchor1_type),
      anchor_gene = as.character(sel$anchor1_gene),
      stringsAsFactors = FALSE
    ),
    data.frame(
      anchor_id = as.character(sel$a2_id),
      anchor_type = as.character(sel$anchor2_type),
      anchor_gene = as.character(sel$anchor2_gene),
      stringsAsFactors = FALSE
    )
  )
  side_types <- side_types[!is.na(side_types$anchor_id) & !duplicated(side_types$anchor_id), , drop = FALSE]
  idx_t <- match(am$anchor_id, side_types$anchor_id)
  hit <- !is.na(idx_t)
  am$anchor_type[hit] <- side_types$anchor_type[idx_t[hit]]
  am$anchor_gene[hit] <- side_types$anchor_gene[idx_t[hit]]
  am
}

#' Internal: Merge chromatin evidence into the anchor map
#'
#' @return A list with \code{am} (enriched anchor map) and
#'   \code{evidence_active} (logical: was a validation table supplied).
#' @keywords internal
#' @noRd
.clique_merge_evidence <- function(am, validation_res, require_chromatin_evidence,
                                   quiet) {
  evidence_active <- !is.null(validation_res)
  am$chromatin_evidence <- NA_character_
  am$H3K4me3 <- NA
  if (!evidence_active) {
    if (isTRUE(require_chromatin_evidence) && !quiet) {
      warning(
        "require_chromatin_evidence = TRUE but no validation_res was ",
        "supplied (see validate_epeG_by_chromatin). Chromatin evidence ",
        "filtering is inactive: anchors are clustered by type alone. ",
        "Anchors absent from a supplied validation_res are treated as ",
        "lacking evidence.",
        call. = FALSE
      )
    }
    return(list(am = am, evidence_active = FALSE))
  }
  if (!is.data.frame(validation_res) ||
    !"anchor_id" %in% colnames(validation_res) ||
    !"enhancer_evidence" %in% colnames(validation_res)) {
    stop("`validation_res` must be the output of ",
      "validate_epeG_by_chromatin(): a data frame with 'anchor_id' and ",
      "'enhancer_evidence' columns.",
      call. = FALSE
    )
  }
  ev <- validation_res[, c("anchor_id", "enhancer_evidence"), drop = FALSE]
  ev$anchor_id <- as.character(ev$anchor_id)
  ev$chromatin_evidence <- as.character(ev$enhancer_evidence)
  ev$H3K4me3 <- if ("H3K4me3" %in% colnames(validation_res)) {
    as.logical(validation_res$H3K4me3)
  } else {
    rep(NA, nrow(validation_res))
  }
  idx_e <- match(am$anchor_id, ev$anchor_id)
  hit_e <- !is.na(idx_e)
  am$chromatin_evidence[hit_e] <- ev$chromatin_evidence[idx_e[hit_e]]
  am$H3K4me3[hit_e] <- ev$H3K4me3[idx_e[hit_e]]
  list(am = am, evidence_active = TRUE)
}

#' Internal: Classify anchors into effective network roles
#'
#' Implements the boundary rules for enhancer-mode and promoter-mode
#' clustering:
#' \itemize{
#'   \item E/dual anchors enter the enhancer network; P/dual anchors enter
#'     the promoter network.
#'   \item eP/eG anchors are promoted by chromatin evidence: enhancer marks
#'     (canonical/strong/supported) promote them into the enhancer network;
#'     an H3K4me3-positive signal promotes them into the promoter network.
#'   \item Without evidence, eP/eG anchors are dropped by default
#'     (conservative) unless \code{include_unresolved = TRUE}.
#' }
#' @return A list with the enriched anchor map (columns
#'   \code{effective_type}, \code{keep}, \code{promotion_reason}) and
#'   boundary counts.
#' @keywords internal
#' @noRd
.clique_classify_anchors <- function(am, enhancer_mode, promoter_mode,
                                     require_chromatin_evidence,
                                     promote_evidenced_eP,
                                     include_unresolved,
                                     evidence_active) {
  n <- nrow(am)
  fam_E <- am$anchor_type %in% c("E", "dual")
  fam_P <- am$anchor_type %in% c("P", "dual")
  fam_e <- am$anchor_type %in% c("eP", "eG")

  ev_ok <- if (evidence_active) {
    am$chromatin_evidence %in% c("canonical", "strong", "supported")
  } else {
    rep(FALSE, n)
  }
  prom_ok <- if (evidence_active) {
    !is.na(am$H3K4me3) & am$H3K4me3
  } else {
    rep(FALSE, n)
  }

  keep <- rep(FALSE, n)
  effective <- rep(NA_character_, n)
  reason <- rep(NA_character_, n)
  drop_category <- rep(NA_character_, n)

  # --- enhancer-mode: base-eligible anchors (E / dual) ---
  if (enhancer_mode) {
    keep[fam_E] <- TRUE
    effective[fam_E] <- "E"
    if (evidence_active && require_chromatin_evidence) {
      drop_e <- fam_E & !ev_ok
      keep[drop_e] <- FALSE
      reason[drop_e] <- "insufficient_chromatin_evidence"
      drop_category[drop_e] <- "insufficient_chromatin_evidence"
    }
    # eP/eG branch
    if (evidence_active) {
      promo <- fam_e & isTRUE(promote_evidenced_eP) & ev_ok
      keep[promo] <- TRUE
      effective[promo] <- "E"
      reason[promo] <- "promoted_by_enhancer_marks"
      if (!require_chromatin_evidence) {
        rest <- fam_e & !promo
        keep[rest] <- TRUE
        effective[rest] <- am$anchor_type[rest]
        reason[rest] <- "kept_without_evidence_requirement"
      } else {
        rest <- fam_e & !promo
        cat_rest <- ifelse(prom_ok[rest] & ev_ok[rest],
          "insufficient_chromatin_evidence",
          ifelse(prom_ok[rest],
            "promoter_signal_not_enhancer",
            "insufficient_chromatin_evidence"
          )
        )
        reason[rest] <- cat_rest
        drop_category[rest] <- cat_rest
      }
    } else if (include_unresolved) {
      keep[fam_e] <- TRUE
      effective[fam_e] <- am$anchor_type[fam_e]
      reason[fam_e] <- "unresolved_no_chromatin_evidence"
    }
  }

  # --- promoter-mode: base-eligible anchors (P / dual) ---
  if (promoter_mode) {
    keep[fam_P] <- TRUE
    effective[fam_P] <- "P"
    # eP/eG branch
    if (evidence_active) {
      promo <- fam_e & prom_ok
      keep[promo] <- TRUE
      effective[promo] <- "P"
      reason[promo] <- "promoted_by_promoter_marks"
      reroute <- fam_e & !prom_ok & ev_ok
      reason[reroute] <- "rerouted_to_enhancer"
      drop_category[reroute] <- "rerouted_to_enhancer"
      if (!require_chromatin_evidence) {
        rest <- fam_e & !promo & !reroute
        keep[rest] <- TRUE
        effective[rest] <- am$anchor_type[rest]
        reason[rest] <- "kept_without_evidence_requirement"
      }
    } else if (include_unresolved) {
      keep[fam_e] <- TRUE
      effective[fam_e] <- am$anchor_type[fam_e]
      reason[fam_e] <- "unresolved_no_chromatin_evidence"
    }
  }

  # Any eligible-but-excluded anchor that was not set above (e.g. an anchor
  # whose type is outside every active loop family, such as G, or an
  # eP/eG with no evidence under a strict rule) is recorded for traceability.
  unresolved_no_ev <- fam_e & !evidence_active & !include_unresolved
  drop_category[unresolved_no_ev] <- "unresolved_no_chromatin_evidence"
  reason[unresolved_no_ev] <- "unresolved_no_chromatin_evidence"
  not_family <- !keep & is.na(drop_category)
  drop_category[not_family] <- "not_in_loop_family"
  reason[not_family] <- "not_in_loop_family"

  # --- boundary counters (derived from the drop categories, so counts
  #     always agree with the traceability table) ---
  drop_cats <- factor(drop_category[!keep],
    levels = c(
      "insufficient_chromatin_evidence", "promoter_signal_not_enhancer",
      "rerouted_to_enhancer", "unresolved_no_chromatin_evidence",
      "not_in_loop_family"
    )
  )
  drop_tab <- tabulate(drop_cats, nbins = 5L)
  counts <- list(
    n_dropped_insufficient_evidence = drop_tab[1L],
    n_promoter_signal_not_enhancer = drop_tab[2L],
    n_rerouted_to_enhancer = drop_tab[3L],
    n_unresolved_dropped = drop_tab[4L],
    n_not_in_family = drop_tab[5L],
    n_promoted_to_enhancer = if (enhancer_mode) {
      sum(fam_e & isTRUE(promote_evidenced_eP) & ev_ok)
    } else {
      0L
    },
    n_promoted_to_promoter = if (promoter_mode) {
      sum(fam_e & prom_ok)
    } else {
      0L
    }
  )

  am$keep <- keep
  am$effective_type <- effective
  am$promotion_reason <- reason
  am$drop_category <- drop_category
  list(am = am, counts = counts)
}

#' Internal: Build the clique graph and membership
#'
#' @return A list with the edge data frame (deduplicated), the igraph
#'   object, clique assignments, community assignments, per-clique
#'   internal edge counts, and the Leiden partition modularity (NA when
#'   \code{community_method = "components"} or when the partition fails).
#' @keywords internal
#' @noRd
.clique_build_graph <- function(am, sel, min_clique_size, community_method,
                                resolution_parameter) {
  am_keep <- am[am$keep, , drop = FALSE]
  keep_lookup <- setNames(am$keep, am$anchor_id)

  pairs <- unique(data.frame(
    a1 = as.character(sel$a1_id),
    a2 = as.character(sel$a2_id),
    stringsAsFactors = FALSE
  ))
  pairs <- pairs[!is.na(pairs$a1) & !is.na(pairs$a2), , drop = FALSE]
  pairs <- pairs[pairs$a1 != pairs$a2, , drop = FALSE]
  ok <- keep_lookup[pairs$a1] & keep_lookup[pairs$a2]
  ok[is.na(ok)] <- FALSE
  edges <- pairs[ok, , drop = FALSE]
  colnames(edges) <- c("from", "to")

  n_vertices <- nrow(am_keep)
  if (n_vertices == 0L || nrow(edges) == 0L) {
    g <- igraph::make_empty_graph(n = n_vertices, directed = FALSE)
    if (n_vertices > 0L) {
      igraph::V(g)$name <- am_keep$anchor_id
    }
    return(list(
      edges = edges, g = g,
      clique_of = setNames(rep(NA_integer_, n_vertices), am_keep$anchor_id),
      community_of = setNames(rep(NA_integer_, n_vertices), am_keep$anchor_id),
      modularity = NA_real_
    ))
  }

  g <- igraph::graph_from_data_frame(edges, directed = FALSE,
    vertices = am_keep$anchor_id
  )
  comp <- igraph::components(g)
  clique_of <- comp$membership
  sizes <- tabulate(clique_of)
  small <- which(sizes < min_clique_size)
  if (length(small) > 0L) {
    clique_of[clique_of %in% small] <- NA_integer_
  }

  community_of <- rep(NA_integer_, igraph::vcount(g))
  names(community_of) <- igraph::V(g)$name
  modularity <- NA_real_
  if (community_method == "leiden" && igraph::ecount(g) > 0L) {
    # Modularity-objective Leiden; igraph >= 2.1 uses `resolution` while
    # older versions use `resolution_parameter` -- try both, silently.
    comm <- tryCatch(
      igraph::cluster_leiden(igraph::simplify(g),
        objective_function = "modularity",
        resolution = resolution_parameter
      ),
      error = function(e) NULL
    )
    if (is.null(comm)) {
      comm <- tryCatch(
        igraph::cluster_leiden(igraph::simplify(g),
          objective_function = "modularity",
          resolution_parameter = resolution_parameter
        ),
        error = function(e) NULL
      )
    }
    if (!is.null(comm)) {
      # membership is ordered by vertex ID; assign positionally.
      community_of <- as.integer(comm$membership)
      names(community_of) <- igraph::V(g)$name
      modularity <- tryCatch(
        igraph::modularity(igraph::simplify(g), as.integer(comm$membership)),
        error = function(e) NA_real_
      )
    }
  }

  list(
    edges = edges, g = g,
    clique_of = clique_of,
    community_of = community_of,
    modularity = modularity
  )
}

#' Internal: Assemble the traceability table of anchors excluded from the network
#'
#' Records every anchor that was not clustered, with its type, chromatin
#' evidence grade, and the reason (drop category) for exclusion. This keeps
#' the boundary decisions fully auditable even when anchors are removed.
#' @keywords internal
#' @noRd
.clique_dropped_table <- function(am) {
  empty <- data.frame(
    anchor_id = character(), chr = character(),
    start = integer(), end = integer(),
    anchor_type = character(), chromatin_evidence = character(),
    drop_category = character(),
    stringsAsFactors = FALSE
  )
  if (is.null(am) || is.null(am$keep) || is.null(am$drop_category)) {
    return(empty)
  }
  bad <- am[!am$keep & !is.na(am$drop_category), , drop = FALSE]
  if (nrow(bad) == 0L) {
    return(empty)
  }
  data.frame(
    anchor_id = bad$anchor_id,
    chr = bad$chr,
    start = bad$start,
    end = bad$end,
    anchor_type = bad$anchor_type,
    chromatin_evidence = bad$chromatin_evidence,
    drop_category = bad$drop_category,
    stringsAsFactors = FALSE
  )
}

#' Internal: Assemble the membership table
#' @keywords internal
#' @noRd
.clique_membership_table <- function(am_keep, graph_res) {
  ids <- am_keep$anchor_id
  clique_rank <- graph_res$clique_of[ids]
  clique_id <- ifelse(is.na(clique_rank),
    NA_character_,
    paste0("Clique_", sprintf("%02d", clique_rank))
  )
  data.frame(
    anchor_id = ids,
    chr = am_keep$chr,
    start = am_keep$start,
    end = am_keep$end,
    anchor_type = am_keep$anchor_type,
    anchor_gene = am_keep$anchor_gene,
    effective_type = am_keep$effective_type,
    chromatin_evidence = am_keep$chromatin_evidence,
    promotion_reason = am_keep$promotion_reason,
    clique_id = clique_id,
    community_id = as.integer(graph_res$community_of[ids]),
    intra_clique_degree = .clique_intra_degree(ids, clique_rank, graph_res$edges),
    stringsAsFactors = FALSE
  )
}

#' Internal: Per-anchor degree within its own clique (local connectivity)
#'
#' Counts, for each kept anchor, the number of clique-graph edges to other
#' anchors of the same clique. This is an anchor-local, parameter-free measure
#' of how connected an anchor is within the module it belongs to -- as opposed
#' to the annotation-layer, whole-network labels. Anchors in cliques dropped
#' by \code{min_clique_size} (clique_rank NA) get degree \code{NA}.
#' @param anchor_ids Character vector.
#' @param clique_rank Named integer vector (anchor -> clique rank, NA for
#'   clique-lost anchors).
#' @param edges Data frame with columns \code{from}, \code{to}.
#' @return Integer vector parallel to \code{anchor_ids}.
#' @keywords internal
#' @noRd
.clique_intra_degree <- function(anchor_ids, clique_rank, edges) {
  out <- rep(NA_integer_, length(anchor_ids))
  names(out) <- anchor_ids
  if (nrow(edges) == 0L) {
    return(out)
  }
  keep_edges <- !is.na(clique_rank[edges$from]) &
    !is.na(clique_rank[edges$to]) &
    clique_rank[edges$from] == clique_rank[edges$to]
  e <- edges[keep_edges, , drop = FALSE]
  if (nrow(e) == 0L) {
    return(out)
  }
  deg <- tabulate(c(match(e$from, anchor_ids), match(e$to, anchor_ids)),
    nbins = length(anchor_ids)
  )
  out[] <- deg
  out
}

#' Internal: Build gene-link rows with recycled role/evidence labels
#' @keywords internal
#' @noRd
.clique_gene_rows <- function(loop_ID, gene, anchor_role, evidence) {
  n <- length(loop_ID)
  data.frame(
    loop_ID = loop_ID, gene = gene,
    anchor_role = rep(anchor_role, n),
    evidence = rep(evidence, n),
    stringsAsFactors = FALSE
  )
}

#' Internal: Compute per-clique gene links
#'
#' @return A data frame with columns \code{clique_id}, \code{gene},
#'   \code{gene_role}, \code{loop_ID}, \code{anchor_role},
#'   \code{evidence}, \code{n_paths}.
#' @keywords internal
#' @noRd
.clique_genes_table <- function(membership, loop_df, neighbor_hop,
                                compute_anchors, compute_linked) {
  out <- data.frame(
    clique_id = character(), gene = character(), gene_role = character(),
    loop_ID = character(), anchor_role = character(),
    evidence = character(), n_paths = integer(),
    stringsAsFactors = FALSE
  )
  kept <- membership[!is.na(membership$clique_id), , drop = FALSE]
  if (nrow(kept) == 0L) {
    return(out)
  }

  # --- anchor-gene rows (promoter hubs) ---
  if (compute_anchors) {
    ag <- kept[!is.na(kept$anchor_gene) & kept$anchor_gene != "", , drop = FALSE]
    if (nrow(ag) > 0L) {
      gene_vec <- unlist(strsplit(as.character(ag$anchor_gene), ";"))
      genes <- trimws(gene_vec)
      clique_vec <- rep(ag$clique_id, lengths(strsplit(as.character(ag$anchor_gene), ";")))
      valid <- !is.na(genes) & genes != ""
      anchor_rows <- data.frame(
        clique_id = clique_vec[valid],
        gene = genes[valid],
        gene_role = "anchor_gene",
        loop_ID = NA_character_,
        anchor_role = "clique_anchor",
        evidence = "anchor_gene",
        n_paths = 1L,
        stringsAsFactors = FALSE
      )
      out <- rbind(out, anchor_rows)
    }
  }

  # --- linked-promoter rows (enhancer cliques) ---
  if (compute_linked) {
    clique_ids <- unique(kept$clique_id)

    # Precompute, once, the promoter-side link rows and per-anchor row
    # indexes. Keeps the per-clique query O(clique rows) instead of
    # O(total rows). Index selection via sort(unlist(...)) reproduces the
    # original loop_df row order exactly, so output is byte-identical to a
    # full-table scan.
    a2prom <- .is_target_promoter_like(loop_df$anchor2_type) &
      !is.na(loop_df$anchor2_gene) & loop_df$anchor2_gene != ""
    a1prom <- .is_target_promoter_like(loop_df$anchor1_type) &
      !is.na(loop_df$anchor1_gene) & loop_df$anchor1_gene != ""
    idx_by_a1 <- split(which(a2prom), as.character(loop_df$a1_id)[a2prom])
    idx_by_a2 <- split(which(a1prom), as.character(loop_df$a2_id)[a1prom])

    .clique_rows_at <- function(anchor_ids, idx_map) {
      present <- anchor_ids[anchor_ids %in% names(idx_map)]
      if (length(present) == 0L) {
        return(NULL)
      }
      sort(unlist(idx_map[present], use.names = FALSE))
    }

    adjacency <- unique(data.frame(
      a1 = as.character(loop_df$a1_id),
      a2 = as.character(loop_df$a2_id),
      stringsAsFactors = FALSE
    ))
    adjacency <- adjacency[!is.na(adjacency$a1) & !is.na(adjacency$a2) &
      adjacency$a1 != adjacency$a2, , drop = FALSE]

    expanded_map <- list()
    if (neighbor_hop >= 1L && nrow(adjacency) > 0L) {
      g_full <- igraph::graph_from_data_frame(adjacency, directed = FALSE)
      for (cid in clique_ids) {
        anchors_c <- kept$anchor_id[kept$clique_id == cid]
        ego_hits <- igraph::ego(g_full, order = 1L, nodes = anchors_c, mode = "all")
        nb <- unique(unlist(lapply(ego_hits, function(x) {
          igraph::V(g_full)$name[as.integer(x)]
        })))
        expanded_map[[cid]] <- setdiff(nb, anchors_c)
      }
    }

    link_rows <- lapply(clique_ids, function(cid) {
      anchors_c <- kept$anchor_id[kept$clique_id == cid]
      expanded <- if (neighbor_hop >= 1L) expanded_map[[cid]] else character(0)

      i1 <- .clique_rows_at(anchors_c, idx_by_a1)
      i2 <- .clique_rows_at(expanded, idx_by_a1)
      i3 <- .clique_rows_at(anchors_c, idx_by_a2)
      i4 <- .clique_rows_at(expanded, idx_by_a2)

      cand1 <- loop_df[i1, , drop = FALSE]
      exp1 <- loop_df[i2, , drop = FALSE]
      cand2 <- loop_df[i3, , drop = FALSE]
      exp2 <- loop_df[i4, , drop = FALSE]

      direct <- rbind(
        .clique_gene_rows(cand1$loop_ID, cand1$anchor2_gene,
          "clique_anchor", "loop_anchor"
        ),
        .clique_gene_rows(cand2$loop_ID, cand2$anchor1_gene,
          "clique_anchor", "loop_anchor"
        )
      )
      expanded_rows <- rbind(
        .clique_gene_rows(exp1$loop_ID, exp1$anchor2_gene,
          "expanded_anchor", "expanded"
        ),
        .clique_gene_rows(exp2$loop_ID, exp2$anchor1_gene,
          "expanded_anchor", "expanded"
        )
      )
      rows <- rbind(direct, expanded_rows)
      if (is.null(rows) || nrow(rows) == 0L) {
        return(NULL)
      }
      gene_vec <- unlist(strsplit(as.character(rows$gene), ";"))
      genes <- trimws(gene_vec)
      row_idx <- rep(seq_len(nrow(rows)), lengths(strsplit(as.character(rows$gene), ";")))
      valid <- !is.na(genes) & genes != ""
      data.frame(
        clique_id = cid,
        gene = genes[valid],
        gene_role = "linked_promoter",
        loop_ID = rows$loop_ID[row_idx[valid]],
        anchor_role = rows$anchor_role[row_idx[valid]],
        evidence = rows$evidence[row_idx[valid]],
        n_paths = 1L,
        stringsAsFactors = FALSE
      )
    })
    link_rows <- link_rows[!vapply(link_rows, is.null, logical(1))]
    if (length(link_rows) > 0L) {
      linked <- do.call(rbind, link_rows)
      out <- rbind(out, linked)
    }
  }

  if (nrow(out) == 0L) {
    return(out)
  }
  # Collapse loop IDs per (clique, gene, role, evidence)
  key <- paste(out$clique_id, out$gene, out$gene_role, out$anchor_role, out$evidence, sep = "\r")
  loop_ids <- tapply(out$loop_ID, key, function(x) {
    x <- unique(stats::na.omit(x))
    paste(x, collapse = ";")
  })
  n_paths <- tapply(out$n_paths, key, sum)
  first <- match(unique(key), key)
  data.frame(
    clique_id = out$clique_id[first],
    gene = out$gene[first],
    gene_role = out$gene_role[first],
    loop_ID = unname(loop_ids[unique(key)]),
    anchor_role = out$anchor_role[first],
    evidence = out$evidence[first],
    n_paths = unname(n_paths[unique(key)]),
    stringsAsFactors = FALSE
  )
}

#' Internal: Overlay a BED of genomic features (peaks/variants) onto clique anchors
#'
#' Returns a long, traceable table of which peptide peaks overlap which clique
#' anchor. This is descriptive (no per-clique statistics) -- see
#' \code{.clique_peak_enrichment} for the aggregate matched-background test.
#' @return A data frame with \code{clique_id}, \code{peak_id},
#'   \code{anchor_id}, \code{overlap_bp}; empty schema when no peaks given.
#' @keywords internal
#' @noRd
.clique_peaks_table <- function(kept, peak_gr) {
  empty <- data.frame(
    clique_id = character(), peak_id = character(),
    anchor_id = character(), overlap_bp = integer(),
    stringsAsFactors = FALSE
  )
  if (is.null(peak_gr) || length(peak_gr) == 0L || nrow(kept) == 0L) {
    return(empty)
  }
  anchors_gr <- GenomicRanges::makeGRangesFromDataFrame(
    data.frame(
      seqnames = kept$chr, start = kept$start, end = kept$end,
      anchor_id = kept$anchor_id, clique_id = kept$clique_id,
      stringsAsFactors = FALSE
    ),
    keep.extra.columns = TRUE
  )
  hits <- GenomicRanges::findOverlaps(anchors_gr, peak_gr)
  if (length(hits) == 0L) {
    return(empty)
  }
  q <- anchors_gr[S4Vectors::queryHits(hits)]
  s <- peak_gr[S4Vectors::subjectHits(hits)]
  data.frame(
    clique_id = q$clique_id,
    peak_id = paste0("Peak_", S4Vectors::subjectHits(hits)),
    anchor_id = q$anchor_id,
    overlap_bp = GenomicRanges::width(GenomicRanges::pintersect(q, s)),
    stringsAsFactors = FALSE
  )
}

#' Internal: Network-level clustering significance against a degree-preserving null
#'
#' Tests whether the observed E-E/E-P module structure is more modular than a
#' random network with the same degree sequence (degree-preserving rewiring,
#' \code{igraph::rewire}). The observed Leiden modularity (modularity objective)
#' is compared against modularity values of the same partition applied to
#' rewired graphs; the empirical p-value is the fraction of null modularities
#' at least as large as observed (with a +1 pseudocount). This answers the
#' structural anti-noise question only -- per-clique significance is NOT
#' reported (multiple-testing of 1000+ cliques and size confounding make it
#' indefensible); clique relevance is supported by the functional layers
#' (peak enrichment, perturbation response).
#' @return A one-row data frame, or \code{NULL} when the graph is too small.
#' @keywords internal
#' @noRd
.clique_clustering_pvalue <- function(g, n_perm = 200L, seed = NULL) {
  sg <- igraph::simplify(g)
  if (igraph::vcount(sg) < 3L || igraph::ecount(sg) < 2L) {
    return(NULL)
  }
  if (!is.null(seed)) {
    withr::local_seed(seed)
  }
  comm <- tryCatch(
    igraph::cluster_leiden(sg, objective_function = "modularity"),
    error = function(e) NULL
  )
  if (is.null(comm)) {
    return(NULL)
  }
  mem <- as.integer(comm$membership)
  q_obs <- igraph::modularity(sg, mem)
  niter <- max(igraph::ecount(g), 100L)
  q_null <- vapply(seq_len(n_perm), function(i) {
    rw <- tryCatch(
      igraph::rewire(g, with = igraph::keeping_degseq(niter = niter)),
      error = function(e) NULL
    )
    if (is.null(rw)) {
      return(NA_real_)
    }
    tryCatch(igraph::modularity(igraph::simplify(rw), mem),
      error = function(e) NA_real_
    )
  }, numeric(1))
  q_null <- q_null[!is.na(q_null)]
  if (length(q_null) < 10L) {
    return(NULL)
  }
  p <- (sum(q_null >= q_obs) + 1) / (length(q_null) + 1)
  data.frame(
    statistic = "leiden_modularity",
    observed_modularity = q_obs,
    null_mean_modularity = mean(q_null),
    null_sd_modularity = stats::sd(q_null),
    p_value = p,
    n_permutations = length(q_null),
    n_nodes = igraph::vcount(sg),
    n_edges = igraph::ecount(sg),
    stringsAsFactors = FALSE
  )
}

#' Internal: Matched-background permutation enrichment of clique anchors
#'   against a genomic feature (peak/variant) set
#'
#' Collapses the enhancer-clique anchor intervals into a query region set and
#' runs a permutation test (\code{regioneR::overlapPermTest}) against the
#' feature set, randomising the query across the genome (optionally masked by
#' a mappability/blacklist BED). An optional 1D comparison region set
#' (e.g. super-enhancers) is benchmarked identically so 3D-clique enrichment
#' can be contrasted against a 1D baseline. Fails soft when \pkg{regioneR} is
#' unavailable. This is an exploratory, background-matched enrichment, not a
#' per-clique significance test.
#' @return A data frame with one row per region set, or \code{NULL}.
#' @keywords internal
#' @noRd
.clique_peak_enrichment <- function(kept, peak_bed, comparison_region_bed,
                                    peak_ntimes, peak_mappability_bed,
                                    annotation_res) {
  if (!requireNamespace("regioneR", quietly = TRUE)) {
    warning("regioneR is required for peak enrichment; skipping enrichment.",
      call. = FALSE
    )
    return(NULL)
  }
  if (is.null(kept) || nrow(kept) == 0L) {
    return(NULL)
  }
  peak_gr <- read_simple_bed(peak_bed, quiet = TRUE)
  if (length(peak_gr) == 0L) {
    warning("peak_bed contains no features; skipping enrichment.", call. = FALSE)
    return(NULL)
  }
  mask_gr <- if (!is.null(peak_mappability_bed)) {
    read_simple_bed(peak_mappability_bed, quiet = TRUE)
  } else {
    NULL
  }

  genome_df <- .clique_enrichgenome(annotation_res, kept)
  if (is.null(genome_df)) {
    return(NULL)
  }

  .run_one <- function(label, kept_sub) {
    if (is.null(kept_sub) || nrow(kept_sub) == 0L) {
      return(NULL)
    }
    query_gr <- unique(GenomicRanges::makeGRangesFromDataFrame(
      data.frame(
        seqnames = kept_sub$chr, start = kept_sub$start, end = kept_sub$end,
        stringsAsFactors = FALSE
      )
    ))
    per <- .with_known_upstream_noise_suppressed(
      regioneR::overlapPermTest(
        A = query_gr, B = peak_gr,
        ntimes = peak_ntimes,
        genome = genome_df,
        mask = mask_gr,
        alternative = "greater",
        count.once = TRUE,
        force.parallel = FALSE
      )
    )
    res_t <- per$numOverlaps
    .scalar_or_na <- function(x) {
      if (length(x) == 0L) NA_real_ else as.numeric(x[1L])
    }
    obs <- .scalar_or_na(res_t$observed)
    perm_vec <- as.numeric(res_t$permuted)
    exp <- if (length(perm_vec) > 0L) stats::median(perm_vec) else NA_real_
    data.frame(
      region_set = label,
      n_regions = length(query_gr),
      n_observed_overlaps = obs,
      mean_expected = exp,
      fold_enrichment = if (!is.na(exp) && exp > 0) obs / exp else NA_real_,
      z_score = .scalar_or_na(res_t$zscore),
      p_value = .scalar_or_na(res_t$pval),
      n_permutations = peak_ntimes,
      stringsAsFactors = FALSE
    )
  }

  rows <- list()
  rows[[1L]] <- .run_one("enhancer_clique_anchors", kept)
  if (!is.null(comparison_region_bed)) {
    ref_gr <- read_simple_bed(comparison_region_bed, quiet = TRUE)
    if (length(ref_gr) > 0L) {
      ref_df <- data.frame(
        chr = as.character(GenomicRanges::seqnames(ref_gr)),
        start = GenomicRanges::start(ref_gr),
        end = GenomicRanges::end(ref_gr),
        stringsAsFactors = FALSE
      )
      rows[[2L]] <- .run_one(basename(comparison_region_bed), ref_df)
    }
  }
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    return(NULL)
  }
  do.call(rbind, rows)
}

#' Internal: Build a genome data frame for regioneR from the anchor registry
#' @keywords internal
#' @noRd
.clique_enrichgenome <- function(annotation_res, kept) {
  reg <- tryCatch(.get_anchor_registry(annotation_res), error = function(e) NULL)
  if (!is.null(reg)) {
    sl <- GenomeInfoDb::seqlengths(reg)
    sl <- sl[!is.na(sl) & sl > 0]
    if (length(sl) > 0L) {
      return(data.frame(
        chr = names(sl), start = rep(1L, length(sl)), end = as.integer(sl),
        stringsAsFactors = FALSE
      ))
    }
  }
  # fallback: infer chromosome sizes from observed anchor extents
  chrs <- unique(kept$chr)
  if (length(chrs) == 0L) {
    return(NULL)
  }
  mx <- vapply(chrs, function(c) {
    max(c(
      kept$end[kept$chr == c],
      rep(1L, length(c))
    ), na.rm = TRUE)
  }, numeric(1))
  data.frame(chr = chrs, start = rep(1L, length(chrs)),
    end = as.integer(mx), stringsAsFactors = FALSE)
}


#' Internal: Assemble the cliques summary table
#' @keywords internal
#' @noRd
.clique_summary_table <- function(membership, edges, clique_genes,
                                  evidence_active, clique_peaks = NULL) {
  kept <- membership[!is.na(membership$clique_id), , drop = FALSE]
  empty_peaks <- is.null(clique_peaks) || nrow(clique_peaks) == 0L
  n_peaks_by_clique <- if (empty_peaks) {
    setNames(numeric(0), character(0))
  } else {
    stats::setNames(
      table(clique_peaks$clique_id),
      as.character(unique(clique_peaks$clique_id))
    )
  }
  if (nrow(kept) == 0L) {
    return(data.frame(
      clique_id = character(), n_anchors = integer(),
      n_internal_loops = integer(),
      chr = character(), start = integer(), end = integer(),
      is_trans_chromosomal = logical(),
      n_chromatin_validated = integer(), weakest_evidence = character(),
      mean_intra_degree = numeric(), max_intra_degree = integer(),
      anchor_genes = character(), linked_promoter_genes = character(),
      n_genes = integer(), n_target_peaks = integer(),
      stringsAsFactors = FALSE
    ))
  }

  edge_clique <- if (nrow(edges) > 0L) {
    cl1 <- membership$clique_id[match(edges$from, membership$anchor_id)]
    cl2 <- membership$clique_id[match(edges$to, membership$anchor_id)]
    both <- !is.na(cl1) & !is.na(cl2) & cl1 == cl2
    data.frame(clique_id = cl1[both], stringsAsFactors = FALSE)
  } else {
    data.frame(clique_id = character(), stringsAsFactors = FALSE)
  }

  evidence_rank <- c(
    canonical = 1L, strong = 2L, supported = 3L,
    limited = 4L, uncertain = 5L
  )

  clique_ids <- unique(kept$clique_id)
  per <- lapply(clique_ids, function(cid) {
    m <- kept[kept$clique_id == cid, , drop = FALSE]
    chrs <- unique(m$chr)
    trans <- length(chrs) > 1L
    n_valid <- if (evidence_active) {
      sum(m$chromatin_evidence %in% c("canonical", "strong", "supported"),
        na.rm = TRUE
      )
    } else {
      NA_integer_
    }
    ranks <- evidence_rank[m$chromatin_evidence]
    weakest <- if (length(stats::na.omit(ranks)) > 0L) {
      names(evidence_rank)[max(ranks, na.rm = TRUE)]
    } else {
      NA_character_
    }
    genes <- clique_genes$gene[clique_genes$clique_id == cid &
      clique_genes$gene_role == "anchor_gene"]
    linked <- clique_genes$gene[clique_genes$clique_id == cid &
      clique_genes$gene_role == "linked_promoter"]
    data.frame(
      clique_id = cid,
      n_anchors = nrow(m),
      n_internal_loops = sum(edge_clique$clique_id == cid),
      chr = paste(sort(chrs), collapse = ";"),
      start = if (trans) NA_integer_ else min(m$start),
      end = if (trans) NA_integer_ else max(m$end),
      is_trans_chromosomal = trans,
      n_chromatin_validated = n_valid,
      weakest_evidence = weakest,
      mean_intra_degree = mean(m$intra_clique_degree[!is.na(m$intra_clique_degree)]),
      max_intra_degree = if (any(!is.na(m$intra_clique_degree))) {
        max(m$intra_clique_degree, na.rm = TRUE)
      } else {
        NA_integer_
      },
      anchor_genes = if (length(genes) > 0L) {
        paste(sort(unique(genes)), collapse = ";")
      } else {
        NA_character_
      },
      linked_promoter_genes = if (length(linked) > 0L) {
        paste(sort(unique(linked)), collapse = ";")
      } else {
        NA_character_
      },
      n_genes = length(unique(c(genes, linked))),
      n_target_peaks = if (empty_peaks) {
        NA_integer_
      } else {
        unname(n_peaks_by_clique[cid])
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, per)
}

#' Internal: Empty clique result with stable schema
#' @keywords internal
#' @noRd
.empty_clique_result <- function(qc_status, qc_reason, loop_types,
                                 community_method, resolution_parameter,
                                 gene_mode,
                                 peak_bed = NULL, run_peak_enrichment = FALSE,
                                 clustering_significance = FALSE,
                                 clustering_ntimes = 200L) {
  membership <- data.frame(
    anchor_id = character(), chr = character(),
    start = integer(), end = integer(),
    anchor_type = character(), anchor_gene = character(),
    effective_type = character(), chromatin_evidence = character(),
    promotion_reason = character(), clique_id = character(),
    community_id = integer(), intra_clique_degree = integer(),
    stringsAsFactors = FALSE
  )
  cliques <- data.frame(
    clique_id = character(), n_anchors = integer(),
    n_internal_loops = integer(),
    chr = character(), start = integer(), end = integer(),
    is_trans_chromosomal = logical(),
    n_chromatin_validated = integer(), weakest_evidence = character(),
    mean_intra_degree = numeric(), max_intra_degree = integer(),
    anchor_genes = character(), linked_promoter_genes = character(),
    n_genes = integer(),
    stringsAsFactors = FALSE
  )
  clique_genes <- data.frame(
    clique_id = character(), gene = character(), gene_role = character(),
    loop_ID = character(), anchor_role = character(),
    evidence = character(), n_paths = integer(),
    stringsAsFactors = FALSE
  )
  dropped <- data.frame(
    anchor_id = character(), chr = character(),
    start = integer(), end = integer(),
    anchor_type = character(), chromatin_evidence = character(),
    drop_category = character(),
    stringsAsFactors = FALSE
  )
  clique_peaks <- data.frame(
    clique_id = character(), peak_id = character(),
    anchor_id = character(), overlap_bp = integer(),
    stringsAsFactors = FALSE
  )
  list(
    cliques = cliques,
    membership = membership,
    clique_genes = clique_genes,
    dropped = dropped,
    clique_peaks = clique_peaks,
    peak_enrichment = NULL,
    clustering_test = NULL,
    network = igraph::make_empty_graph(n = 0L, directed = FALSE),
    summary = list(
      n_input_loops = 0L, n_edges = 0L, n_cliques = 0L,
      clique_size_median = NA_real_, clique_size_max = 0L,
      n_anchors_eligible = 0L,
      n_dropped = 0L,
      n_dropped_insufficient_evidence = 0L,
      n_promoter_signal_not_enhancer = 0L,
      n_rerouted_to_enhancer = 0L,
      n_promoted_to_enhancer = 0L,
      n_promoted_to_promoter = 0L,
      n_unresolved_dropped = 0L,
      n_not_in_family = 0L,
      n_peak_overlaps = 0L,
      clustering_p = NA_real_,
      modularity = NA_real_
    ),
    metadata = .build_looplook_metadata(
      fun = "cluster_enhancer_cliques",
      params = list(
        loop_types = loop_types,
        community_method = community_method,
        resolution_parameter = resolution_parameter,
        gene_mode = gene_mode,
        peak_bed = peak_bed,
        run_peak_enrichment = run_peak_enrichment,
        clustering_significance = clustering_significance,
        clustering_ntimes = clustering_ntimes,
        qc_status = qc_status,
        qc_reason = qc_reason
      )
    )
  )
}

#' Cluster enhancer-like anchors into spatial cliques or promoter hubs
#'
#' Extracts enhancer-like (E/eP family) or promoter (P family) loops from
#' an annotated (optionally refined) looplook result, builds the anchor
#' network, and identifies enhancer cliques (dense E-E modules) or
#' promoter hubs (P-P modules) via connected components, with an optional
#' Leiden community refinement.
#'
#' @details
#' \strong{Refinement inheritance:} the function consumes
#' \code{annotate_peaks_and_loops()} output directly, or its
#' expression-/chromatin-refined variants. Reclassified anchor types
#' (e.g., P -> eP from expression refinement) are inherited from the
#' supplied object; the function never re-runs refinements. The
#' recommended chain is annotation -> \code{refine_loop_anchors_by_expression()}
#' -> \code{refine_loop_anchors_by_chromatin()} (or
#' \code{validate_epeG_by_chromatin()}) -> clustering.
#'
#' \strong{Which mode?} One call runs one mode, selected by
#' \code{loop_types}:
#' \itemize{
#'   \item \strong{Enhancer cliques} (default): \code{loop_types =
#'     c("E-E","E-eP","eP-eP")}. Asks \emph{which distal enhancers cluster in
#'     3D, and which promoters they loop to}. \code{clique_genes} is populated
#'     with \code{gene_role = "linked_promoter"} and \code{cliques$anchor_genes}
#'     is \code{NA}.
#'   \item \strong{Promoter hubs}: \code{loop_types = c("P-P","eP-P")}. Asks
#'     \emph{which promoters cluster in space}. \code{clique_genes} is populated
#'     with \code{gene_role = "anchor_gene"} (the anchor's own genes).
#' }
#' The two modes are independent calls and are exported as separate workbooks;
#' run the function twice to compare.
#'
#' \strong{Boundary rules:}
#' \itemize{
#'   \item E/dual anchors enter the enhancer network; P/dual anchors enter
#'     the promoter network.
#'   \item eP/eG anchors are dual-identity: with enhancer marks
#'     (canonical/strong/supported) they are promoted into the enhancer
#'     network; with an H3K4me3-positive signal they are promoted into the
#'     promoter network.
#'   \item Anchors absent from \code{validation_res} are treated as lacking
#'     evidence. Without evidence, eP/eG anchors are dropped by default
#'     (conservative); set \code{include_unresolved = TRUE} for structural
#'     exploration.
#'   \item When \code{require_chromatin_evidence = TRUE}, E anchors without
#'     sufficient enhancer marks are dropped from the enhancer network.
#' }
#'
#' @param annotation_res List. Output of
#'   \code{\link{annotate_peaks_and_loops}} or a refined variant.
#' @param validation_res Data frame or \code{NULL}. Output of
#'   \code{\link{validate_epeG_by_chromatin}}, providing per-anchor
#'   \code{enhancer_evidence} grades and \code{H3K4me3} overlap status.
#'   Default \code{NULL} (evidence filtering inactive with a warning when
#'   \code{require_chromatin_evidence = TRUE}).
#' @param loop_types Character vector. Loop types to cluster. Default
#'   \code{c("E-E", "E-eP", "eP-eP")} for enhancer cliques; use
#'   \code{c("P-P", "eP-P")} for promoter hubs.
#' @param require_chromatin_evidence Logical. If \code{TRUE}, only anchors
#'   with enhancer marks (canonical/strong/supported) enter the enhancer
#'   network. Default \code{TRUE}.
#' @param promote_evidenced_eP Logical. Promote eP/eG anchors with enhancer
#'   marks into the enhancer network. Default \code{TRUE}.
#' @param community_method Character. \code{"components"} (default;
#'   connected components, consistent with the package's cluster
#'   semantics) or \code{"leiden"} (additional modularity-based
#'   refinement recorded in \code{community_id} and
#'   \code{summary$modularity}). The Leiden refinement is an exploratory
#'   partition of an unweighted binary interaction graph -- it is not a
#'   significance test and may be sensitive to \code{resolution_parameter}
#'   on sparse chain-like networks; components remain the primary,
#'   interpretation-grade grouping in \code{clique_id}.
#' @param resolution_parameter Numeric. Leiden resolution parameter,
#'   used only when \code{community_method = "leiden"}. Values > 1 favour
#'   finer partitions, values in (0, 1) coarser ones. Default \code{1}.
#'   This is an exploratory granularity control, not a statistical
#'   significance threshold.
#' @param min_clique_size Integer. Minimum number of anchors per clique
#'   (component). Default \code{2}.
#' @param neighbor_hop Integer. \code{0} (default): clique genes are
#'   promoters directly looped to clique anchors. \code{1}: additionally
#'   includes promoters reached through one extra anchor hop.
#' @param include_unresolved Logical. If \code{TRUE}, eP/eG anchors without
#'   chromatin evidence stay in the network under their original type
#'   (structural exploration). Default \code{FALSE}.
#' @param gene_mode Character. \code{"auto"} (default): anchor genes for
#'   promoter-mode types, linked promoter genes for enhancer-mode types
#'   (both when loop types mix families). \code{"anchors"} /
#'   \code{"linked"} force one behaviour.
#' @param peak_bed Character or \code{NULL}. Optional BED of genomic features
#'   (GWAS/eQTL variants, TF ChIP peaks, ATAC regions, etc.) to overlay onto
#'   clique anchors. When supplied, produces the descriptive
#'   \code{clique_peaks} table and the \code{n_target_peaks} column. This is
#'   an annotation layer only -- it does not affect anchor inclusion or gene
#'   assignment. Default \code{NULL}.
#' @param run_peak_enrichment Logical. If \code{TRUE} (requires
#'   \code{peak_bed}), an aggregate matched-background permutation test
#'   (\pkg{regioneR} \code{overlapPermTest}) is run on the merged clique-anchor
#'   regions vs the feature set, reporting fold-enrichment and p-value. This is
#'   not a per-clique significance test; per-clique confidence is descriptive
#'   only. Fails soft when \pkg{regioneR} is unavailable or the genome cannot
#'   be resolved. Default \code{FALSE}.
#' @param peak_ntimes Integer. Number of permutations for
#'   \code{run_peak_enrichment}. Default \code{1000L}.
#' @param peak_mappability_bed Character or \code{NULL}. Optional BED mask
#'   (e.g. mappability/blacklist) used to constrain randomisation of the
#'   enrichment background. Default \code{NULL}.
#' @param comparison_region_bed Character or \code{NULL}. Optional BED of a 1D
#'   baseline region set (e.g. super-enhancers) run through the identical
#'   enrichment procedure, so 3D-clique enrichment can be contrasted with a 1D
#'   baseline. Default \code{NULL}.
#' @param clustering_significance Logical. If \code{TRUE}, a network-level
#'   anti-noise test is run: the observed Leiden modularity (modularity
#'   objective) of the clique graph is compared with a degree-preserving
#'   rewiring null (\code{igraph::rewire(keeping_degseq)}). The resulting
#'   empirical p-value (\code{summary$clustering_p}, \code{clustering_test})
#'   answers whether the E-E/E-P module structure is more modular than random
#'   wiring with the same degrees. It is a structural check only: per-clique
#'   p-values are deliberately NOT computed (multiple testing over hundreds to
#'   thousands of cliques and size confounding), and clique relevance is
#'   supported by the functional layers (peak enrichment and perturbation
#'   response). Default \code{FALSE}.
#' @param clustering_ntimes Integer. Number of rewiring permutations for
#'   \code{clustering_significance}. Default \code{200L}.
#' @param clustering_seed Integer or \code{NULL}. Optional seed for the
#'   rewiring permutation, applied with \code{withr::local_seed} so the
#'   caller's RNG state is restored on exit (no leakage). \code{NULL}
#'   (default) uses the current RNG state, matching the package convention in
#'   \code{profile_target_genes()}.
#' @param write_output Logical. If \code{TRUE}, write an Excel workbook
#'   with the Cliques, Membership and CliqueGenes sheets to
#'   \code{out_dir}. Default \code{FALSE}.
#' @param out_dir Character. Output directory (used when
#'   \code{write_output = TRUE}). Default \code{"./results"}.
#' @param project_name Character. File-name prefix for exported workbooks.
#'   Default \code{"EnhancerCliques"}.
#' @param quiet Logical. If \code{TRUE}, suppress progress messages.
#'   Default \code{FALSE}.
#'
#' @return A named list:
#' \itemize{
#'   \item \code{cliques} -- one row per clique: size, span, trans-chromosomal
#'     status, chromatin evidence composition, internal connectivity
#'     (\code{mean_intra_degree}, \code{max_intra_degree}), and genes.
#'   \item \code{membership} -- one row per eligible anchor: coordinates,
#'     original and effective types, chromatin evidence, promotion reason,
#'     clique/community assignment, and \code{intra_clique_degree} (number of
#'     edges to other anchors of the same clique; anchor-local connectivity,
#'     not a whole-network label).
#'   \item \code{clique_genes} -- long-format clique-to-gene provenance
#'     (\code{gene_role}, \code{anchor_role}, \code{evidence},
#'     \code{n_paths}).
#'   \item \code{dropped} -- traceability table of anchors excluded from the
#'     network (\code{anchor_id}, coordinates, \code{anchor_type},
#'     \code{chromatin_evidence}, \code{drop_category}: one of
#'     \code{insufficient_chromatin_evidence},
#'     \code{promoter_signal_not_enhancer}, \code{rerouted_to_enhancer},
#'     \code{unresolved_no_chromatin_evidence}, \code{not_in_loop_family}).
#'     Kept in the output so every boundary decision remains auditable.
#'   \item \code{clique_peaks} -- descriptive overlay of \code{peak_bed}
#'     features onto clique anchors (\code{clique_id}, \code{peak_id},
#'     \code{anchor_id}, \code{overlap_bp}); \code{NULL} when \code{peak_bed}
#'     is not supplied.
#'   \item \code{peak_enrichment} -- aggregate matched-background permutation
#'     enrichment of merged clique-anchor regions against \code{peak_bed},
#'     optionally alongside a \code{comparison_region_bed} baseline
#'     (\code{region_set}, \code{n_regions}, \code{n_observed_overlaps},
#'     \code{mean_expected}, \code{fold_enrichment}, \code{z_score},
#'     \code{p_value}); \code{NULL} unless \code{run_peak_enrichment = TRUE}.
#'   \item \code{clustering_test} -- network-level degree-preserving rewiring
#'     test of Leiden modularity (\code{observed_modularity},
#'     \code{null_mean_modularity}, \code{null_sd_modularity},
#'     \code{p_value}, \code{n_permutations}); \code{NULL} unless
#'     \code{clustering_significance = TRUE}. Structural anti-noise check only;
#'     not a per-clique test.
#'   \item \code{network} -- igraph object with vertex attributes
#'     \code{anchor_id}, \code{anchor_type}, \code{effective_type},
#'     \code{chromatin_evidence}, \code{clique}, \code{community}.
#'   \item \code{summary} -- clustering and boundary-rule counts
#'     (\code{n_dropped}, \code{n_dropped_insufficient_evidence},
#'     \code{n_promoter_signal_not_enhancer}, \code{n_rerouted_to_enhancer},
#'     \code{n_unresolved_dropped}, \code{n_not_in_family}, promotions), plus
#'     \code{modularity} (NA unless \code{community_method = "leiden"};
#'     reported for partition-quality assessment, not as a significance
#'     test).
#'   \item \code{metadata} -- parameters, QC status and boundary rules.
#' }
#'
#' @importFrom igraph make_empty_graph components cluster_leiden
#'   graph_from_data_frame V ego
#' @export
#'
#' @examples
#' loop_df <- data.frame(
#'   loop_ID = c("L1", "L2", "L3", "L4"),
#'   chr1 = c("chr1", "chr1", "chr2", "chr1"),
#'   start1 = c(100L, 200L, 50L, 400L),
#'   end1 = c(200L, 300L, 150L, 500L),
#'   chr2 = c("chr1", "chr1", "chr2", "chr1"),
#'   start2 = c(210L, 310L, 160L, 510L),
#'   end2 = c(310L, 410L, 260L, 610L),
#'   anchor1_type = c("E", "E", "E", "eP"),
#'   anchor2_type = c("E", "E", "E", "E"),
#'   anchor1_gene = NA_character_,
#'   anchor2_gene = NA_character_,
#'   loop_type = c("E-E", "E-E", "E-E", "E-eP"),
#'   a1_id = c("A1", "A2", "A3", "A4"),
#'   a2_id = c("A2", "A3", "A4", "A5"),
#'   stringsAsFactors = FALSE
#' )
#' res <- list(loop_annotation = loop_df)
#' cl <- cluster_enhancer_cliques(res,
#'   require_chromatin_evidence = FALSE, quiet = TRUE
#' )
#' head(cl$cliques)
#'
#' # Promoter hubs: switch loop_types (and this time show anchor genes)
#' pp_df <- data.frame(
#'   loop_ID = c("L1", "L2"),
#'   chr1 = c("chr1", "chr1"),
#'   start1 = c(100L, 210L), end1 = c(200L, 310L),
#'   chr2 = c("chr1", "chr1"),
#'   start2 = c(210L, 320L), end2 = c(310L, 420L),
#'   anchor1_type = c("P", "P"), anchor2_type = c("P", "P"),
#'   anchor1_gene = c("GENE1", "GENE2"),
#'   anchor2_gene = c("GENE2", "GENE3"),
#'   loop_type = "P-P",
#'   a1_id = c("B1", "B2"), a2_id = c("B2", "B3"),
#'   stringsAsFactors = FALSE
#' )
#' ph <- cluster_enhancer_cliques(list(loop_annotation = pp_df),
#'   loop_types = c("P-P", "eP-P"), quiet = TRUE
#' )
#' head(ph$cliques)
cluster_enhancer_cliques <- function(
  annotation_res,
  validation_res = NULL,
  loop_types = c("E-E", "E-eP", "eP-eP"),
  require_chromatin_evidence = TRUE,
  promote_evidenced_eP = TRUE,
  community_method = c("components", "leiden"),
  resolution_parameter = 1,
  min_clique_size = 2,
  neighbor_hop = 0,
  include_unresolved = FALSE,
  gene_mode = c("auto", "anchors", "linked"),
  peak_bed = NULL,
  run_peak_enrichment = FALSE,
  peak_ntimes = 1000L,
  peak_mappability_bed = NULL,
  comparison_region_bed = NULL,
  clustering_significance = FALSE,
  clustering_ntimes = 200L,
  clustering_seed = NULL,
  write_output = FALSE,
  out_dir = "./results",
  project_name = "EnhancerCliques",
  quiet = FALSE
) {
  community_method <- match.arg(community_method)
  gene_mode <- match.arg(gene_mode)
  log_message <- .make_log_message(quiet)
  .assert_scalar_number(resolution_parameter, "resolution_parameter", min = 0)
  .assert_scalar_count(min_clique_size, "min_clique_size", min = 1)
  .assert_scalar_count(neighbor_hop, "neighbor_hop")
  if (!is.null(peak_bed)) {
    .assert_file_exists(peak_bed, "peak_bed")
  }
  if (!is.null(comparison_region_bed)) {
    .assert_file_exists(comparison_region_bed, "comparison_region_bed")
  }
  if (!is.null(peak_mappability_bed)) {
    .assert_file_exists(peak_mappability_bed, "peak_mappability_bed")
  }
  if (isTRUE(run_peak_enrichment) && is.null(peak_bed)) {
    warning("run_peak_enrichment = TRUE but no peak_bed was supplied; ",
      "enrichment is skipped.", call. = FALSE)
    run_peak_enrichment <- FALSE
  }
  if (isTRUE(clustering_significance)) {
    .assert_scalar_count(clustering_ntimes, "clustering_ntimes", min = 10)
    if (!is.null(clustering_seed)) {
      .assert_scalar_count(clustering_seed, "clustering_seed", min = 1)
    }
  }
  if (neighbor_hop > 1L) {
    stop("`neighbor_hop` must be 0 or 1.", call. = FALSE)
  }
  if (!is.character(loop_types) || length(loop_types) == 0L ||
    anyNA(loop_types) || !all(nzchar(trimws(loop_types)))) {
    stop("`loop_types` must be a non-empty character vector of loop-type codes.",
      call. = FALSE
    )
  }

  base_reason <- paste(
    "Clustering requested for loop_types:",
    paste(loop_types, collapse = ", ")
  )

  if (.is_empty_annotation_result(annotation_res)) {
    return(.empty_clique_result("empty_input", base_reason,
      loop_types = loop_types, community_method = community_method,
      resolution_parameter = resolution_parameter,
      gene_mode = gene_mode,
      peak_bed = peak_bed, run_peak_enrichment = run_peak_enrichment,
      clustering_significance = clustering_significance,
      clustering_ntimes = clustering_ntimes
    ))
  }
  loop_df <- annotation_res$loop_annotation
  if (is.null(loop_df) || !is.data.frame(loop_df)) {
    stop("'annotation_res$loop_annotation' is missing.", call. = FALSE)
  }
  required <- c("a1_id", "a2_id", "anchor1_type", "anchor2_type", "loop_type")
  missing_cols <- setdiff(required, colnames(loop_df))
  if (length(missing_cols) > 0L) {
    stop("'loop_annotation' is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  sel <- loop_df[loop_df$loop_type %in% loop_types, , drop = FALSE]
  sel <- sel[!is.na(sel$a1_id) & !is.na(sel$a2_id) &
    sel$a1_id != "" & sel$a2_id != "", , drop = FALSE]
  if (nrow(sel) == 0L) {
    return(.empty_clique_result("no_matching_loops", base_reason,
      loop_types = loop_types, community_method = community_method,
      resolution_parameter = resolution_parameter,
      gene_mode = gene_mode,
      peak_bed = peak_bed, run_peak_enrichment = run_peak_enrichment,
      clustering_significance = clustering_significance,
      clustering_ntimes = clustering_ntimes
    ))
  }

  # --- mode detection ---
  .tokens <- function(x) {
    strsplit(as.character(x), "-", fixed = TRUE)
  }
  enhancer_mode <- any(vapply(.tokens(loop_types), function(t) {
    any(t == "E")
  }, logical(1)))
  promoter_mode <- any(vapply(.tokens(loop_types), function(t) {
    any(t == "P")
  }, logical(1)))
  if (!enhancer_mode && !promoter_mode) {
    enhancer_mode <- TRUE
  }
  if (gene_mode == "auto") {
    compute_anchors <- promoter_mode
    compute_linked <- enhancer_mode
  } else {
    compute_anchors <- gene_mode == "anchors"
    compute_linked <- gene_mode == "linked"
  }

  # --- anchor map + evidence ---
  am <- .clique_anchor_map(annotation_res, sel)
  if (nrow(am) == 0L) {
    return(.empty_clique_result("no_valid_anchors", base_reason,
      loop_types = loop_types, community_method = community_method,
      resolution_parameter = resolution_parameter,
      gene_mode = gene_mode,
      peak_bed = peak_bed, run_peak_enrichment = run_peak_enrichment,
      clustering_significance = clustering_significance,
      clustering_ntimes = clustering_ntimes
    ))
  }
  ev <- .clique_merge_evidence(am, validation_res,
    require_chromatin_evidence, quiet
  )
  am <- ev$am
  evidence_active <- ev$evidence_active

  cls <- .clique_classify_anchors(
    am = am, enhancer_mode = enhancer_mode, promoter_mode = promoter_mode,
    require_chromatin_evidence = require_chromatin_evidence,
    promote_evidenced_eP = promote_evidenced_eP,
    include_unresolved = include_unresolved,
    evidence_active = evidence_active
  )
  am <- cls$am
  counts <- cls$counts
  dropped <- .clique_dropped_table(am)
  if (!any(am$keep)) {
    return(.empty_clique_result("no_eligible_anchors", base_reason,
      loop_types = loop_types, community_method = community_method,
      resolution_parameter = resolution_parameter,
      gene_mode = gene_mode,
      peak_bed = peak_bed, run_peak_enrichment = run_peak_enrichment,
      clustering_significance = clustering_significance,
      clustering_ntimes = clustering_ntimes
    ))
  }

  # --- graph + membership ---
  graph_res <- .clique_build_graph(am, sel, min_clique_size, community_method,
    resolution_parameter
  )
  am_keep <- am[am$keep, , drop = FALSE]
  membership <- .clique_membership_table(am_keep, graph_res)
  log_message(
    ">>> ", nrow(membership), " eligible anchors, ",
    nrow(graph_res$edges), " edges, ",
    sum(!is.na(unique(membership$clique_id))), " clique(s) kept."
  )

  # --- clique genes ---
  clique_genes <- .clique_genes_table(
    membership = membership, loop_df = loop_df,
    neighbor_hop = neighbor_hop,
    compute_anchors = compute_anchors,
    compute_linked = compute_linked
  )

  # --- peak overlay (optional, descriptive) ---
  peak_gr <- if (!is.null(peak_bed)) {
    read_simple_bed(peak_bed, quiet = TRUE)
  } else {
    NULL
  }
  clique_peaks <- NULL
  if (!is.null(peak_gr) && length(peak_gr) > 0L) {
    clique_peaks <- .clique_peaks_table(
      membership[!is.na(membership$clique_id), , drop = FALSE], peak_gr
    )
    log_message(
      ">>> ", nrow(clique_peaks),
      " peak-anchor overlap(s) recorded across ",
      length(unique(clique_peaks$clique_id)), " clique(s)."
    )
  }

  # --- cliques table ---
  cliques <- .clique_summary_table(
    membership = membership, edges = graph_res$edges,
    clique_genes = clique_genes, evidence_active = evidence_active,
    clique_peaks = clique_peaks
  )
  cliques <- cliques[order(as.numeric(sub("Clique_", "", cliques$clique_id))), , drop = FALSE]

  # --- network-level clustering significance (optional anti-noise test) ---
  clustering_test <- NULL
  if (isTRUE(clustering_significance)) {
    t0 <- Sys.time()
    clustering_test <- .clique_clustering_pvalue(
      g = graph_res$g, n_perm = clustering_ntimes, seed = clustering_seed
    )
    if (!is.null(clustering_test)) {
      log_message(
        ">>> Network clustering significance: modularity = ",
        round(clustering_test$observed_modularity, 4),
        " vs null ", round(clustering_test$null_mean_modularity, 4),
        " (p = ", round(clustering_test$p_value, 4), ")",
        " [", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s]"
      )
    }
  }

  # --- aggregate matched-background enrichment (optional) ---
  peak_enrichment <- NULL
  if (isTRUE(run_peak_enrichment)) {
    kept_for_enrichment <- membership[!is.na(membership$clique_id), , drop = FALSE]
    peak_enrichment <- .clique_peak_enrichment(
      kept = kept_for_enrichment,
      peak_bed = peak_bed,
      comparison_region_bed = comparison_region_bed,
      peak_ntimes = peak_ntimes,
      peak_mappability_bed = peak_mappability_bed,
      annotation_res = annotation_res
    )
  }

  # --- summary ---
  clique_sizes <- cliques$n_anchors
  summary <- list(
    n_input_loops = nrow(sel),
    n_anchors_eligible = nrow(membership),
    n_edges = nrow(graph_res$edges),
    n_cliques = nrow(cliques),
    clique_size_median = if (length(clique_sizes) > 0L) {
      stats::median(clique_sizes)
    } else {
      NA_real_
    },
    clique_size_max = if (length(clique_sizes) > 0L) max(clique_sizes) else 0L,
    n_trans_chromosomal = sum(cliques$is_trans_chromosomal),
    n_dropped = nrow(dropped),
    n_dropped_insufficient_evidence = counts$n_dropped_insufficient_evidence,
    n_promoter_signal_not_enhancer = counts$n_promoter_signal_not_enhancer,
    n_rerouted_to_enhancer = counts$n_rerouted_to_enhancer,
    n_promoted_to_enhancer = counts$n_promoted_to_enhancer,
    n_promoted_to_promoter = counts$n_promoted_to_promoter,
    n_unresolved_dropped = counts$n_unresolved_dropped,
    n_not_in_family = counts$n_not_in_family,
    n_peak_overlaps = if (!is.null(clique_peaks)) nrow(clique_peaks) else 0L,
    clustering_p = if (!is.null(clustering_test)) clustering_test$p_value else NA_real_,
    modularity = graph_res$modularity
  )

  # --- network annotation ---
  g <- graph_res$g
  if (igraph::vcount(g) > 0L) {
    ids <- igraph::V(g)$name
    igraph::V(g)$anchor_type <- am_keep$anchor_type[match(ids, am_keep$anchor_id)]
    igraph::V(g)$effective_type <- am_keep$effective_type[match(ids, am_keep$anchor_id)]
    igraph::V(g)$chromatin_evidence <- am_keep$chromatin_evidence[match(ids, am_keep$anchor_id)]
    igraph::V(g)$clique <- graph_res$clique_of[ids]
    igraph::V(g)$community <- graph_res$community_of[ids]
  }

  # --- metadata ---
  input_refinement <- if (is.list(annotation_res$metadata) &&
    is.character(annotation_res$metadata$function_name)) {
    annotation_res$metadata$function_name
  } else {
    "unknown"
  }
  boundary_rules <- paste0(
    "enhancer_mode=", enhancer_mode,
    "; promoter_mode=", promoter_mode,
    "; require_chromatin_evidence=", require_chromatin_evidence,
    "; promote_evidenced_eP=", promote_evidenced_eP,
    "; include_unresolved=", include_unresolved,
    "; validation_provided=", evidence_active,
    "; absent_from_validation_treated_as_no_evidence=TRUE"
  )
  metadata <- .build_looplook_metadata(
    fun = "cluster_enhancer_cliques",
    params = list(
      loop_types = loop_types,
      require_chromatin_evidence = require_chromatin_evidence,
      promote_evidenced_eP = promote_evidenced_eP,
      community_method = community_method,
      resolution_parameter = resolution_parameter,
      min_clique_size = min_clique_size,
      neighbor_hop = neighbor_hop,
      include_unresolved = include_unresolved,
      gene_mode = gene_mode,
      peak_bed = peak_bed,
      run_peak_enrichment = run_peak_enrichment,
      peak_ntimes = peak_ntimes,
      comparison_region_bed = comparison_region_bed,
      clustering_significance = clustering_significance,
      clustering_ntimes = clustering_ntimes,
      clustering_seed = clustering_seed,
      input_refinement = input_refinement,
      boundary_rules = boundary_rules,
      qc_status = "ok"
    ),
    diagnostics = summary
  )

  out <- list(
    cliques = cliques,
    membership = membership,
    clique_genes = clique_genes,
    dropped = dropped,
    clique_peaks = clique_peaks,
    peak_enrichment = peak_enrichment,
    clustering_test = clustering_test,
    network = g,
    summary = summary,
    metadata = metadata
  )

  # --- export ---
  if (isTRUE(write_output)) {
    .ensure_out_dir(TRUE, out_dir)
    wb <- openxlsx::createWorkbook()
    .add_sheet(wb, "Cliques", cliques)
    .add_sheet(wb, "Membership", membership)
    .add_sheet(wb, "CliqueGenes", clique_genes)
    .add_sheet(wb, "Dropped", dropped)
    .add_sheet(wb, "CliquePeaks", clique_peaks)
    .add_sheet(wb, "PeakEnrichment", peak_enrichment)
    .save_workbook(wb, out_dir, project_name,
      suffix = "_Cliques.xlsx",
      fail_prefix = "Failed to save clique workbook: "
    )
  }

  out
}

# --- Gene co-expression validation for clique gene sets ---

#' Internal: Resolve gene coordinates for distance-matched null sampling
#'
#' Accepts a data frame with columns \code{gene_id}, \code{chr},
#' \code{start}, \code{end} (gene IDs matching the expression matrix row
#' names / clique gene symbols, case-insensitive), or a TxDb object
#' (coordinates extracted from \code{GenomicFeatures::genes()}; note the
#' resulting gene_id column may use a different ID space from gene symbols,
#' so a data frame is preferred).
#'
#' @param gene_coords Data frame or TxDb.
#' @return A data frame with columns \code{gene_id} (uppercase), \code{chr},
#'   \code{mid} (gene mid-point), de-duplicated by gene_id.
#' @keywords internal
#' @noRd
.clique_gene_coords <- function(gene_coords) {
  if (inherits(gene_coords, "TxDb")) {
    gr <- .with_known_upstream_noise_suppressed(
      GenomicFeatures::genes(gene_coords)
    )
    df <- data.frame(
      gene_id = as.character(S4Vectors::mcols(gr)$gene_id),
      chr = as.character(GenomicRanges::seqnames(gr)),
      start = GenomicRanges::start(gr),
      end = GenomicRanges::end(gr),
      stringsAsFactors = FALSE
    )
  } else if (is.data.frame(gene_coords)) {
    df <- gene_coords
  } else {
    stop("`gene_coords` must be a data.frame with columns ",
      "gene_id, chr, start, end, or a TxDb object.",
      call. = FALSE
    )
  }
  req <- c("gene_id", "chr", "start", "end")
  missing_cols <- setdiff(req, colnames(df))
  if (length(missing_cols) > 0L) {
    stop("`gene_coords` is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  df$gene_id <- trimws(toupper(as.character(df$gene_id)))
  df$chr <- trimws(as.character(df$chr))
  df$mid <- suppressWarnings(
    as.numeric(df$start) + (as.numeric(df$end) - as.numeric(df$start)) / 2
  )
  df <- df[!is.na(df$gene_id) & df$gene_id != "" &
    !is.na(df$chr) & df$chr != "" & !is.na(df$mid), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("No valid gene coordinates remain after validation.",
      call. = FALSE
    )
  }
  df <- df[!duplicated(df$gene_id), , drop = FALSE]
  df[, c("gene_id", "chr", "mid"), drop = FALSE]
}

#' Internal: Assign a genomic distance to a distance bin
#'
#' @param d Numeric distance (same chromosome), or \code{Inf} for
#'   trans-chromosomal pairs.
#' @return A single bin label.
#' @keywords internal
#' @noRd
.clique_distance_bin <- function(d) {
  if (is.infinite(d)) return("trans")
  if (d < 1e4) return("lt10kb")
  if (d < 1e5) return("10kb_100kb")
  if (d < 1e6) return("100kb_1Mb")
  if (d < 1e7) return("1Mb_10Mb")
  "gt10Mb"
}

#' Internal: Bin-wise lower/upper distance bounds
#' @return A list with \code{lo} and \code{hi} named by bin.
#' @keywords internal
#' @noRd
.clique_bin_bounds <- function() {
  list(
    lt10kb = c(lo = 0, hi = 1e4),
    `10kb_100kb` = c(lo = 1e4, hi = 1e5),
    `100kb_1Mb` = c(lo = 1e5, hi = 1e6),
    `1Mb_10Mb` = c(lo = 1e6, hi = 1e7),
    gt10Mb = c(lo = 1e7, hi = Inf)
  )
}

#' Internal: Vectorised per-row correlation between two expression matrices
#'
#' Computes the correlation of each row pair \code{(A[k, ], B[k, ])}.  Uses
#' the Pearson closed form when possible (fast); falls back to per-row
#' \code{stats::cor()} for Spearman or when missing values are present.
#'
#' @param A Numeric matrix (pairs x samples).
#' @param B Numeric matrix (pairs x samples).
#' @param cor_method Character. \code{"pearson"} or \code{"spearman"}.
#' @return Numeric vector of length nrow(A).
#' @keywords internal
#' @noRd
.pairwise_cor_vec <- function(A, B, cor_method) {
  n <- nrow(A)
  if (n == 0L) return(numeric(0))
  if (cor_method == "spearman" || anyNA(A) || anyNA(B)) {
    return(vapply(seq_len(n), function(k) {
      stats::cor(A[k, , drop = FALSE], B[k, , drop = FALSE],
        method = cor_method, use = "pairwise.complete.obs")[1L, 1L]
    }, numeric(1)))
  }
  Ac <- sweep(A, 1L, rowMeans(A), "-")
  Bc <- sweep(B, 1L, rowMeans(B), "-")
  denom <- sqrt(rowSums(Ac^2) * rowSums(Bc^2))
  out <- rowSums(Ac * Bc) / denom
  out[!is.finite(out)] <- NA_real_
  out
}

#' Internal: Draw a full distance-matched null gene-pair sample
#'
#' Draws \code{bin_counts[[bin]]} random pairs per bin so the null
#' pair-distance distribution matches the observed one.  Per-chromosome
#' sampling is vectorised with a pre-computed eligible-chromosome index per
#' bin.  Genes may repeat across pairs -- the null models the distribution of
#' genomic distances, not of genes.  Returns alternating
#' \code{c(a1, b1, a2, b2, ...)}; pairs that could not be drawn for a bin are
#' filled with \code{NA} entries (the caller skips such permutation draws).
#'
#' @param bin_counts Named integer vector of pairs needed per bin.
#' @param by_chr Named list of per-chromosome data frames (columns
#'   \code{gene_id}, \code{pos}, sorted by \code{pos}).
#' @param eligible_by_bin Named list of per-bin eligible chromosome vectors.
#' @return Character vector (alternating gene pairs), or \code{character(0)}.
#' @keywords internal
#' @noRd
.clique_sample_matched_pairs <- function(bin_counts, by_chr, eligible_by_bin) {
  out_a <- character()
  out_b <- character()
  chrs_all <- names(by_chr)
  bounds_list <- .clique_bin_bounds()
  for (bin in names(bin_counts)) {
    n <- bin_counts[[bin]]
    if (n <= 0L) next
    if (bin == "trans") {
      if (length(chrs_all) < 2L) {
        out_a <- c(out_a, rep(NA_character_, n))
        out_b <- c(out_b, rep(NA_character_, n))
        next
      }
      a <- sample(chrs_all, n, replace = TRUE)
      b <- vapply(a, function(x) {
        rest <- chrs_all[chrs_all != x]
        sample(rest, 1L)
      }, character(1))
      out_a <- c(out_a,
        vapply(a, function(x) sample(by_chr[[x]]$gene_id, 1L), character(1)))
      out_b <- c(out_b,
        vapply(b, function(x) sample(by_chr[[x]]$gene_id, 1L), character(1)))
      next
    }
    elig <- eligible_by_bin[[bin]]
    if (length(elig) == 0L) {
      out_a <- c(out_a, rep(NA_character_, n))
      out_b <- c(out_b, rep(NA_character_, n))
      next
    }
    bounds <- bounds_list[[bin]]
    lo <- unname(bounds[["lo"]])
    hi <- unname(bounds[["hi"]])
    done <- 0L
    attempts <- 0L
    while (done < n && attempts < 40L) {
      attempts <- attempts + 1L
      m <- n - done
      chr_s <- sample(elig, m, replace = TRUE)
      for (ch in unique(chr_s)) {
        d <- by_chr[[ch]]
        pos <- d$pos
        L <- nrow(d)
        k <- which(chr_s == ch)
        ii <- sample.int(L, length(k), replace = TRUE)
        p_i <- pos[ii]
        sj <- findInterval(p_i + lo, pos) + 1L
        ej <- if (is.finite(hi)) findInterval(p_i + hi, pos) else L
        span <- ej - sj + 1L
        span[span < 1L] <- 0L
        jj <- sj + as.integer(stats::runif(length(k)) * span)
        ok <- span > 0L & jj >= sj & jj <= ej & jj != ii
        ok[is.na(jj)] <- FALSE
        if (any(ok)) {
          out_a <- c(out_a, d$gene_id[ii[ok]])
          out_b <- c(out_b, d$gene_id[jj[ok]])
          done <- done + sum(ok)
        }
      }
    }
    if (done < n) {
      out_a <- c(out_a, rep(NA_character_, n - done))
      out_b <- c(out_b, rep(NA_character_, n - done))
    }
  }
  c(rbind(out_a, out_b))
}

#' Internal: Co-expression score for a gene set
#'
#' Mean of all pairwise (off-diagonal) gene-expression correlations within a
#' gene set, computed on the sample-level expression matrix with
#' \code{use = "pairwise.complete.obs"}.  Returns \code{NA_real_} when fewer
#' than two valid genes or no finite correlation is available.
#'
#' @param expr_sub Numeric matrix (genes x samples), genes as rows.
#' @param cor_method Character. \code{"pearson"} or \code{"spearman"}.
#' @return A single numeric score, or \code{NA_real_}.
#' @keywords internal
#' @noRd
.clique_coexpr_score <- function(expr_sub, cor_method) {
  if (is.null(expr_sub) || nrow(expr_sub) < 2L) {
    return(NA_real_)
  }
  cormat <- stats::cor(t(expr_sub), method = cor_method,
    use = "pairwise.complete.obs")
  offdiag <- cormat[lower.tri(cormat)]
  offdiag <- offdiag[!is.na(offdiag)]
  if (length(offdiag) == 0L) {
    return(NA_real_)
  }
  mean(offdiag)
}

#' Internal: Per-clique gene co-expression permutation test
#'
#' For each clique gene set, the observed mean pairwise expression correlation
#' is compared against a null distribution.  Under \code{null_method =
#' "genome"} the null is built by randomly drawing gene sets of equal size
#' from the measurable gene pool.  Under \code{null_method =
#' "distance_matched"} the null is built by drawing gene pairs whose genomic
#' distance distribution matches the observed clique gene pairs (same bin
#' counts per distance bin, including trans-chromosomal pairs), so the test
#' answers whether 3D co-localised genes are more co-expressed than random
#' genes with the same linear proximity.  This answers whether the genes that
#' co-localise in a 3D clique are more transcriptionally coordinated than
#' chance, without assuming any particular expression model.  Like the
#' aggregate peak-enrichment layer, this is a background-matched permutation
#' enrichment -- it does not test individual gene-pair edges and should not be
#' read as per-gene evidence.
#'
#' @param cliques_res List. Output of \code{\link{cluster_enhancer_cliques}}.
#' @param expr_matrix Matrix or data.frame (genes x samples, rownames = gene
#'   IDs) with sample-level expression values, or a character path to such a
#'   file (first column = gene ID).
#' @param gene_role Character vector. Which \code{gene_role} values in
#'   \code{clique_genes} supply the per-clique gene sets.  Typically
#'   \code{"anchor_gene"} (promoter hubs) or \code{"linked_promoter"}
#'   (enhancer cliques).
#' @param n_perm Integer. Number of permutations. Default \code{1000}.
#' @param seed Integer or \code{NULL}. Seed for the permutations, applied with
#'   \code{withr::local_seed} (no RNG leakage).
#' @param cor_method Character. \code{"pearson"} (default) or
#'   \code{"spearman"}.
#' @param min_genes Integer. Minimum number of matched genes required for a
#'   clique to be tested. Default \code{3}.
#' @param quiet Logical. Suppress progress messages.
#' @param null_method Character. \code{"genome"} (default): null gene sets are
#'   drawn uniformly from the measurable gene pool (no genomic-position
#'   matching).  \code{"distance_matched"}: null gene pairs are drawn with the
#'   same same-chromosome distance distribution as the observed clique gene
#'   pairs (and trans-chromosomal pairs matched separately), using
#'   \code{gene_coords}.  The distance-matched null answers whether 3D
#'   co-localised genes are more co-expressed than random genes with the same
#'   linear genomic proximity -- it removes the confound that nearby genes
#'   share local regulatory environments.
#' @param gene_coords Data frame (\code{gene_id}, \code{chr}, \code{start},
#'   \code{end}) or TxDb, required when \code{null_method =
#'   "distance_matched"}.  Ignored otherwise.
#' @return A list with \code{per_clique} (data.frame) and \code{summary}
#'   (named list).
#' @keywords internal
#' @noRd
.clique_gene_coexpression_test <- function(cliques_res, expr_matrix,
                                           gene_role,
                                           n_perm = 1000L, seed = NULL,
                                           cor_method = "pearson",
                                           min_genes = 3L, quiet = FALSE,
                                           null_method = "genome",
                                           gene_coords = NULL) {
  log_message <- .make_log_message(quiet)
  empty <- data.frame(
    clique_id = character(),
    n_genes = integer(), n_matched = integer(), n_pairs = integer(),
    observed_score = numeric(), null_mean = numeric(), null_sd = numeric(),
    z_score = numeric(), p_value = numeric(), n_permutations = integer(),
    stringsAsFactors = FALSE
  )

  if (is.null(cliques_res) || !is.list(cliques_res) ||
    is.null(cliques_res$clique_genes) ||
    !is.data.frame(cliques_res$clique_genes)) {
    stop("`cliques_res` must be the output of cluster_enhancer_cliques() ",
      "containing a `clique_genes` table.",
      call. = FALSE
    )
  }
  cg <- cliques_res$clique_genes
  required_cols <- c("clique_id", "gene", "gene_role")
  missing_cols <- setdiff(required_cols, colnames(cg))
  if (length(missing_cols) > 0L) {
    stop("`clique_genes` is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # --- resolve expression matrix ---
  if (is.character(expr_matrix) && length(expr_matrix) == 1L &&
    nzchar(expr_matrix)) {
    parsed <- .read_expression_matrix(expr_matrix)
    mat <- parsed$matrix
  } else if (is.matrix(expr_matrix) || is.data.frame(expr_matrix)) {
    mat <- as.data.frame(expr_matrix, check.names = FALSE)
    if (is.null(rownames(mat))) {
      stop("`expr_matrix` must have gene IDs as row names.",
        call. = FALSE
      )
    }
    if (!all(vapply(mat, is.numeric, logical(1)))) {
      stop("`expr_matrix` must contain only numeric expression values.",
        call. = FALSE
      )
    }
    # Keep the package-wide non-negative expression contract for both input
    # modes (matrix object and file path go through the same validation).
    if (any(vapply(mat, function(x) any(x < 0, na.rm = TRUE), logical(1)))) {
      stop("Expression matrix contains negative values. This package ",
        "expects non-negative expression (TPM, FPKM, RPKM, raw counts, ",
        "CAGE, or nascent signal).",
        call. = FALSE
      )
    }
    if (any(vapply(mat, function(x) any(is.infinite(x) | is.nan(x), na.rm = TRUE), logical(1)))) {
      stop("Expression matrix contains non-finite values (Inf, -Inf, NaN). ",
        "Replace Inf/NaN with NA or remove those rows.",
        call. = FALSE
      )
    }
  } else {
    stop("`expr_matrix` must be a matrix/data.frame (genes x samples) ",
      "with gene IDs as row names, or a file path.",
      call. = FALSE
    )
  }
  if (ncol(mat) < 2L) {
    stop("Co-expression testing requires at least two sample columns.",
      call. = FALSE
    )
  }

  # --- measurable gene pool ---
  gene_ids <- trimws(rownames(mat))
  if (anyNA(gene_ids) || !all(nzchar(gene_ids))) {
    stop("`expr_matrix` row names contain blank or missing gene IDs.",
      call. = FALSE
    )
  }
  row_ok <- vapply(seq_len(nrow(mat)), function(i) {
    v <- as.numeric(mat[i, ])
    !all(is.na(v)) && stats::var(v, na.rm = TRUE) > 0
  }, logical(1))
  pool_idx <- which(row_ok)
  pool_upper <- toupper(gene_ids[pool_idx])
  if (length(pool_upper) < 2L) {
    return(list(per_clique = empty, summary = list(
      n_cliques = 0L, n_tested = 0L, n_skipped_small = 0L,
      n_skipped_no_match = 0L, n_significant = 0L,
      median_p_value = NA_real_, cor_method = cor_method,
      n_permutations = n_perm, seed = seed
    )))
  }

  # --- distance-matched null preparation (optional) ---
  by_chr <- NULL
  eligible_by_bin <- NULL
  chr_lookup <- NULL
  coord_map <- NULL
  if (null_method == "distance_matched") {
    coords <- .clique_gene_coords(gene_coords)
    coords <- coords[toupper(coords$gene_id) %in% pool_upper, , drop = FALSE]
    if (nrow(coords) < 2L) {
      stop("No gene coordinates matched the measurable gene pool; ",
        "`gene_coords` gene IDs must match expression matrix row names.",
        call. = FALSE
      )
    }
    coords$gene_id <- toupper(coords$gene_id)
    chr_lookup <- stats::setNames(coords$chr, coords$gene_id)
    coord_map <- stats::setNames(coords$mid, coords$gene_id)
    by_chr <- lapply(split(coords, coords$chr), function(d) {
      d <- d[order(d$mid), , drop = FALSE]
      d$pos <- as.numeric(d$mid)
      d[, c("gene_id", "pos")]
    })
    by_chr <- by_chr[vapply(by_chr, nrow, integer(1)) >= 2L]
    bounds_list <- .clique_bin_bounds()
    eligible_by_bin <- lapply(bounds_list, function(b) {
      lo <- unname(b[["lo"]])
      names(by_chr)[vapply(by_chr, function(d) {
        (max(d$pos) - min(d$pos)) >= lo
      }, logical(1))]
    })
  }

  # --- per-clique gene sets ---
  cg_use <- cg[cg$gene_role %in% gene_role, , drop = FALSE]
  if (nrow(cg_use) == 0L) {
    stop("No clique_genes rows match gene_role = ",
      paste(sQuote(gene_role), collapse = ", "),
      ". Did you select the mode that matches the clique type ",
      "('anchor_gene' for promoter hubs, 'linked_promoter' for ",
      "enhancer cliques)?",
      call. = FALSE
    )
  }
  cg_split <- split(trimws(as.character(cg_use$gene)),
    cg_use$clique_id
  )
  clique_ids <- names(cg_split)

  score_fun <- function(genes_upper) {
    m <- match(genes_upper, pool_upper)
    m <- m[!is.na(m)]
    if (length(m) < 2L) {
      return(NA_real_)
    }
    .clique_coexpr_score(as.matrix(mat[pool_idx[m], , drop = FALSE]),
      cor_method
    )
  }

  # Pair-level score for the distance-matched null: mean correlation of a
  # flat vector of alternating gene pairs (a1,b1,a2,b2,...).
  pair_score_fun <- function(pair_genes_upper) {
    if (length(pair_genes_upper) < 2L) return(NA_real_)
    ma <- match(pair_genes_upper[seq.int(1L, length(pair_genes_upper), 2L)],
      pool_upper)
    mb <- match(pair_genes_upper[seq.int(2L, length(pair_genes_upper), 2L)],
      pool_upper)
    if (anyNA(ma) || anyNA(mb)) return(NA_real_)
    cors <- .pairwise_cor_vec(
      as.matrix(mat[pool_idx[ma], , drop = FALSE]),
      as.matrix(mat[pool_idx[mb], , drop = FALSE]),
      cor_method
    )
    cors <- cors[!is.na(cors)]
    if (length(cors) == 0L) NA_real_ else mean(cors)
  }

  bin_names <- c(names(.clique_bin_bounds()), "trans")

  if (!is.null(seed)) {
    withr::local_seed(seed)
  }

  rows <- lapply(clique_ids, function(cid) {
    genes <- unique(na.omit(cg_split[[cid]]))
    genes_upper <- toupper(genes)
    matched <- unique(genes_upper[genes_upper %in% pool_upper])
    n_genes <- length(genes)
    n_matched <- length(matched)
    if (n_matched < min_genes || n_matched < 2L) {
      return(data.frame(
        clique_id = cid, n_genes = n_genes, n_matched = n_matched,
        n_pairs = NA_integer_, observed_score = NA_real_,
        null_mean = NA_real_, null_sd = NA_real_, z_score = NA_real_,
        p_value = NA_real_, n_permutations = n_perm,
        stringsAsFactors = FALSE
      ))
    }
    s_obs <- score_fun(matched)
    if (is.na(s_obs)) {
      return(data.frame(
        clique_id = cid, n_genes = n_genes, n_matched = n_matched,
        n_pairs = NA_integer_, observed_score = NA_real_,
        null_mean = NA_real_, null_sd = NA_real_, z_score = NA_real_,
        p_value = NA_real_, n_permutations = n_perm,
        stringsAsFactors = FALSE
      ))
    }
    if (null_method == "distance_matched" && !is.null(by_chr)) {
      # Observed gene-pair distance bin counts (vectorised)
      idx_ok <- matched[matched %in% coords$gene_id]
      if (length(idx_ok) < 2L) {
        return(data.frame(
          clique_id = cid, n_genes = n_genes, n_matched = n_matched,
          n_pairs = NA_integer_, observed_score = NA_real_,
          null_mean = NA_real_, null_sd = NA_real_, z_score = NA_real_,
          p_value = NA_real_, n_permutations = n_perm,
          stringsAsFactors = FALSE
        ))
      }
      pos_v <- coord_map[idx_ok]
      chr_v <- chr_lookup[idx_ok]
      dmat <- abs(outer(pos_v, pos_v, "-"))
      dmat[outer(chr_v, chr_v, "!=")] <- Inf
      pd <- dmat[lower.tri(dmat)]
      bin_vec <- vapply(pd, .clique_distance_bin, character(1))
      bin_counts <- stats::setNames(integer(length(bin_names)), bin_names)
      tab <- table(bin_vec)
      bin_counts[names(tab)] <- as.integer(tab)
      n_obs_pairs <- sum(bin_counts)
      q_null <- vapply(seq_len(n_perm), function(i) {
        pair_vec <- .clique_sample_matched_pairs(bin_counts, by_chr,
          eligible_by_bin)
        if (length(pair_vec) < 2L) return(NA_real_)
        pair_score_fun(toupper(pair_vec))
      }, numeric(1))
    } else {
      q_null <- vapply(seq_len(n_perm), function(i) {
        s_genes <- sample(pool_upper, size = n_matched, replace = FALSE)
        score_fun(s_genes)
      }, numeric(1))
    }
    q_null <- q_null[!is.na(q_null)]
    if (length(q_null) < 10L) {
      return(data.frame(
        clique_id = cid, n_genes = n_genes, n_matched = n_matched,
        n_pairs = n_matched * (n_matched - 1L) / 2L,
        observed_score = s_obs, null_mean = NA_real_, null_sd = NA_real_,
        z_score = NA_real_, p_value = NA_real_, n_permutations = length(q_null),
        stringsAsFactors = FALSE
      ))
    }
    p <- (sum(q_null >= s_obs) + 1) / (length(q_null) + 1)
    z <- (s_obs - mean(q_null)) / stats::sd(q_null)
    data.frame(
      clique_id = cid, n_genes = n_genes, n_matched = n_matched,
      n_pairs = n_matched * (n_matched - 1L) / 2L,
      observed_score = s_obs, null_mean = mean(q_null),
      null_sd = stats::sd(q_null), z_score = z, p_value = p,
      n_permutations = length(q_null),
      stringsAsFactors = FALSE
    )
  })
  per_clique <- do.call(rbind, rows)
  # Append genes column from the split input
  gene_list <- vapply(per_clique$clique_id, function(cid) {
    g <- cg_split[[cid]]
    if (is.null(g)) return(NA_character_)
    g <- unique(trimws(g[!is.na(g) & nzchar(g)]))
    if (length(g) == 0L) return(NA_character_)
    paste(g, collapse = ";")
  }, character(1))
  per_clique$genes <- gene_list
  per_clique <- per_clique[order(per_clique$p_value, per_clique$clique_id), , drop = FALSE]
  rownames(per_clique) <- NULL

  n_significant <- sum(per_clique$p_value <= 0.05, na.rm = TRUE)
  summary <- list(
    n_cliques = length(clique_ids),
    n_tested = sum(!is.na(per_clique$p_value)),
    n_skipped_small = sum(per_clique$n_matched < min_genes, na.rm = TRUE),
    n_skipped_no_match = sum(
      is.na(per_clique$p_value) & per_clique$n_matched >= min_genes,
      na.rm = TRUE
    ),
    n_significant = n_significant,
    median_p_value = stats::median(per_clique$p_value, na.rm = TRUE),
    cor_method = cor_method,
    n_permutations = n_perm,
    null_method = null_method,
    seed = seed
  )
  log_message("Co-expression test: ", summary$n_tested, "/",
    summary$n_cliques, " cliques tested, ",
    summary$n_significant, " significant (p <= 0.05).")
  list(per_clique = per_clique, summary = summary)
}

#' Test Gene Co-expression within Promoter Hub / Enhancer Clique Gene Sets
#'
#' Background-matched permutation test of whether the genes assigned to each
#' 3D clique (promoter hubs: \code{anchor_gene}; enhancer cliques:
#' \code{linked_promoter}) are more transcriptionally coordinated than random
#' gene sets of equal size drawn from the measurable gene pool.
#'
#' @details
#' \strong{Method.} For each clique the observed statistic is the mean of all
#' pairwise (off-diagonal) gene-expression correlations among the clique's
#' matched genes, computed on sample-level expression values
#' (\code{use = "pairwise.complete.obs"}).  The empirical p-value is the
#' fraction of null scores at least as large as the observed score, with a +1
#' pseudocount (right-sided: we ask whether clique genes are \emph{more}
#' co-expressed than chance).
#'
#' \strong{Null models.} Two null strategies are available:
#' \itemize{
#'   \item \code{null_method = "genome"} (default): each null gene set is
#'     drawn uniformly from the measurable gene pool (genes with finite,
#'     non-zero variance across samples) and the same statistic is recomputed
#'     \code{n_perm} times.  This tests against a genomic-random expectation.
#'   \item \code{null_method = "distance_matched"}: null \emph{gene pairs} are
#'     drawn with the same same-chromosome distance distribution as the
#'     observed clique gene pairs (binned at 10 kb / 100 kb / 1 Mb / 10 Mb,
#'     trans-chromosomal pairs matched separately).  This controls for the
#'     strong genomic prior that linearly neighbouring genes share local
#'     regulatory environments and co-express by position alone, so the test
#'     answers the sharper question: are 3D co-localised genes more
#'     co-expressed than random genes with the same linear proximity?
#'     Requires \code{gene_coords}.
#' }
#'
#' \strong{Scope and interpretation.} This is a background-matched enrichment,
#' not a per-gene or per-gene-pair test.  It supplies the gene-level functional
#' layer for \code{\link{cluster_enhancer_cliques}} output that is symmetric to
#' the peak-enrichment layer.  It is best read as an exploratory validation:
#' cliques whose genes are significantly co-expressed are more likely to reflect
#' coordinated regulation of the assigned genes rather than chance spatial
#' proximity.  No multiple-testing correction is
#' applied in the returned \code{per_clique} p-values; when testing many
#' cliques, apply your own FDR control (e.g. \code{p.adjust(..., "BH")}).
#'
#' \strong{Why not the anchor-level test?} Promoter-anchor positions are
#' determined by gene structure rather than 3D organisation, so position-based
#' enrichment (as used for enhancer cliques) would largely reflect the genomic
#' prior.  Co-expression of the \emph{assigned genes} is the interpretable
#' functional readout for promoter hubs.
#'
#' @param cliques_res List. Output of \code{\link{cluster_enhancer_cliques}}.
#' @param expr_matrix Matrix, data.frame, or character path. Sample-level
#'   expression matrix with gene IDs as row names (first column = gene ID when
#'   a file path).  Non-negative values (TPM/FPKM/counts/CAGE/nascent).
#' @param gene_mode Character. Which clique genes supply the gene sets.
#'   \code{"auto"} (default): uses \code{anchor_gene} rows when present,
#'   otherwise \code{linked_promoter} rows -- this matches the mode used to
#'   build the cliques (promoter hubs vs enhancer cliques).  \code{"anchors"}
#'   forces \code{anchor_gene}; \code{"linked"} forces
#'   \code{linked_promoter}.
#' @param n_perm Integer. Number of permutations. Default \code{1000L}.
#' @param seed Integer or \code{NULL}. Optional seed for the permutations
#'   (applied with \code{withr::local_seed} so the caller's RNG state is
#'   restored).  Default \code{NULL} uses the current RNG state.
#' @param cor_method Character. \code{"pearson"} (default) or
#'   \code{"spearman"}.  Spearman is more robust to outliers and monotone
#'   non-linearity at the cost of power on small gene sets.
#' @param min_genes Integer. Minimum number of expression-matched genes for a
#'   clique to be tested.  Cliques below this are returned with
#'   \code{p_value = NA}. Default \code{3}.
#' @param null_method Character. Null sampling strategy.
#'   \code{"genome"} (default): null gene sets are drawn uniformly from the
#'   measurable gene pool -- no genomic-position matching.  This is the
#'   historical behaviour.
#'   \code{"distance_matched"}: null gene pairs are drawn with the same
#'   same-chromosome distance distribution as the observed clique gene pairs
#'   (trans-chromosomal pairs matched separately), using \code{gene_coords}.
#'   This removes the confound that linearly neighbouring genes share local
#'   regulatory environments, so it answers the sharper question: are 3D
#'   co-localised genes more co-expressed than random genes with the same
#'   linear genomic proximity?  Requires \code{gene_coords}.
#' @param gene_coords Data frame or \code{NULL}. Required when
#'   \code{null_method = "distance_matched"}.  A data frame with columns
#'   \code{gene_id}, \code{chr}, \code{start}, \code{end} (gene IDs matching
#'   the expression matrix row names / clique gene symbols, matched
#'   case-insensitively), or a TxDb object from which coordinates are
#'   extracted via \code{GenomicFeatures::genes()}.  Ignored when
#'   \code{null_method = "genome"}.
#' @param write_output Logical. If \code{TRUE}, write a multi-sheet Excel
#'   workbook to \code{out_dir}. Default \code{FALSE}.
#' @param out_dir Character. Output directory for the Excel file (used when
#'   \code{write_output = TRUE}). Default \code{"./results"}.
#' @param project_name Character. File-name prefix for the exported workbook.
#'   Default \code{"Coexpression"}.
#' @param quiet Logical. Suppress progress messages. Default \code{FALSE}.
#' @return A list:
#' \itemize{
#'   \item \code{per_clique} -- data.frame, one row per clique, ordered by
#'     p-value: \code{clique_id}, \code{n_genes}, \code{n_matched},
#'     \code{n_pairs}, \code{observed_score}, \code{null_mean},
#'     \code{null_sd}, \code{z_score}, \code{p_value},
#'     \code{n_permutations}, \code{genes} (semicolon-delimited gene list).
#'     Cliques with too few matched genes or no finite observed score have
#'     \code{NA} statistics.
#'   \item \code{summary} -- named list: \code{n_cliques}, \code{n_tested},
#'     \code{n_skipped_small}, \code{n_skipped_no_match},
#'     \code{n_significant} (p <= 0.05), \code{median_p_value},
#'     \code{cor_method}, \code{n_permutations}, \code{null_method},
#'     \code{seed}.
#'   \item \code{merged} -- data.frame, \code{per_clique} merged with
#'     \code{cliques_res$cliques} (anchor counts, chromosome, internal loop
#'     count); \code{NULL} when the cliques table is unavailable.  Also
#'     includes an \code{FDR} column (BH-adjusted p-values).
#' }
#' @seealso \code{\link{cluster_enhancer_cliques}},
#'   \code{\link{refine_loop_anchors_by_expression}}
#' @export
#'
#' @examples
#' # Promoter-hub cliques from the enhancer-cliques example
#' pp_df <- data.frame(
#'   loop_ID = c("L1", "L2"),
#'   chr1 = c("chr1", "chr1"),
#'   start1 = c(100L, 210L), end1 = c(200L, 310L),
#'   chr2 = c("chr1", "chr1"),
#'   start2 = c(210L, 320L), end2 = c(310L, 420L),
#'   anchor1_type = c("P", "P"), anchor2_type = c("P", "P"),
#'   anchor1_gene = c("GENE1", "GENE2"),
#'   anchor2_gene = c("GENE2", "GENE3"),
#'   loop_type = "P-P",
#'   a1_id = c("B1", "B2"), a2_id = c("B2", "B3"),
#'   stringsAsFactors = FALSE
#' )
#' ph <- cluster_enhancer_cliques(list(loop_annotation = pp_df),
#'   loop_types = c("P-P", "eP-P"), quiet = TRUE
#' )
#'
#' # Synthetic expression matrix with a co-expressed GENE1/GENE2 pair
#' set.seed(1)
#' fake_expr <- matrix(rnorm(4 * 6), nrow = 4)
#' rownames(fake_expr) <- c("GENE1", "GENE2", "GENE3", "GENE4")
#' colnames(fake_expr) <- paste0("S", 1:6)
#' res <- test_hub_gene_coexpression(ph, fake_expr, n_perm = 100, seed = 1)
#' head(res$per_clique)
test_hub_gene_coexpression <- function(
  cliques_res,
  expr_matrix,
  gene_mode = c("auto", "anchors", "linked"),
  n_perm = 1000L,
  seed = NULL,
  cor_method = c("pearson", "spearman"),
  min_genes = 3L,
  null_method = c("genome", "distance_matched"),
  gene_coords = NULL,
  write_output = FALSE,
  out_dir = "./results",
  project_name = "Coexpression",
  quiet = FALSE
) {
  gene_mode <- match.arg(gene_mode)
  cor_method <- match.arg(cor_method)
  null_method <- match.arg(null_method)
  .assert_scalar_count(n_perm, "n_perm", min = 10L)
  .assert_scalar_count(min_genes, "min_genes", min = 2L)
  if (!is.null(seed)) {
    .assert_scalar_count(seed, "seed", min = 1)
  }
  if (null_method == "distance_matched" && is.null(gene_coords)) {
    stop("`null_method = 'distance_matched'` requires `gene_coords` ",
      "(a data.frame with gene_id, chr, start, end, or a TxDb object).",
      call. = FALSE
    )
  }

  gene_role <- if (gene_mode == "anchors") {
    "anchor_gene"
  } else if (gene_mode == "linked") {
    "linked_promoter"
  } else {
    if (is.null(cliques_res) || !is.list(cliques_res) ||
      is.null(cliques_res$clique_genes) ||
      !is.data.frame(cliques_res$clique_genes)) {
      stop("`cliques_res` must be the output of cluster_enhancer_cliques().",
        call. = FALSE
      )
    }
    roles <- unique(as.character(cliques_res$clique_genes$gene_role))
    if ("anchor_gene" %in% roles) {
      "anchor_gene"
    } else if ("linked_promoter" %in% roles) {
      "linked_promoter"
    } else {
      stop("Could not infer `gene_mode` from `clique_genes`; ",
        "no anchor_gene or linked_promoter rows present. ",
        "Set `gene_mode` explicitly.",
        call. = FALSE
      )
    }
  }

  result <- .clique_gene_coexpression_test(
    cliques_res = cliques_res,
    expr_matrix = expr_matrix,
    gene_role = gene_role,
    n_perm = n_perm,
    seed = seed,
    cor_method = cor_method,
    min_genes = min_genes,
    quiet = quiet,
    null_method = null_method,
    gene_coords = gene_coords
  )

  # Add FDR column
  result$per_clique$FDR <- p.adjust(result$per_clique$p_value, method = "BH")

  # Merge with cliques summary table when available
  merged <- NULL
  if (!is.null(cliques_res$cliques) && is.data.frame(cliques_res$cliques)) {
    merge_cols <- intersect(colnames(cliques_res$cliques),
      c("clique_id", "n_anchors", "chr", "start", "end",
        "n_internal_loops", "anchor_genes", "linked_promoter_genes"))
    if (length(merge_cols) > 1L) {
      merged <- merge(result$per_clique, cliques_res$cliques[, merge_cols, drop = FALSE],
        by = "clique_id", all.x = TRUE)
      merged <- merged[order(merged$FDR, merged$p_value, merged$clique_id), , drop = FALSE]
      rownames(merged) <- NULL
    }
  }
  result$merged <- merged

  # Export to Excel
  if (isTRUE(write_output)) {
    .ensure_out_dir(TRUE, out_dir)
    wb <- openxlsx::createWorkbook()
    if (!is.null(merged)) {
      .add_sheet(wb, "Coexpression_Merged", merged)
    }
    .add_sheet(wb, "Coexpression_PerClique", result$per_clique)
    .add_sheet(wb, "Summary", data.frame(
      metric = names(result$summary),
      value = as.character(unlist(result$summary)),
      stringsAsFactors = FALSE
    ))
    .save_workbook(wb, out_dir, project_name,
      suffix = "_Coexpression.xlsx",
      fail_prefix = "Failed to save coexpression workbook: "
    )
  }

  result
}
