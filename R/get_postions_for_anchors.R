

#' Convert a list of integration anchors into pairwise positions in a dimensional reduction
#'
#' @param anchors An output from FindIntegrationAnchors
#' @param min_score The minimum score for anchor strength to report
#' @param reduction The reduction to use for the position extraction
#'
#' @return A tibble containing pairs of cell IDs along with their xy coordinates
#' @export
#'
#' @examples
get_positions_for_anchors <- function(anchors, min_score=0, reduction="umap") {

  # This is a pain.  The dataset1 and dataset2 values aren't consistent.  If you're
  # comparing 1 and 2 then dataset1 will sometimes be 1 and sometimes 2.  We must
  # therefore group by dataset1 and dataset2 to ensure that we're always dealing
  # with consistently numbered datasets.

  anchors@anchors |>
    dplyr::filter(score>=min_score) |>
    dplyr::group_by(dataset1,dataset2) |>
    dplyr::mutate(cellid1 = colnames(anchors@object.list[[dataset1[1]]])[cell1]) |>
    dplyr::mutate(cellid2 = colnames(anchors@object.list[[dataset2[1]]])[cell2]) |>
    dplyr::ungroup() -> positions



  Seurat::Embeddings(anchors@object.list[[anchors@anchors$dataset1[1]]], reduction = reduction) |>
    tibble::as_tibble(rownames="cellid1") |>
    dplyr::select(1:3) -> cellid1_positions

  colnames(cellid1_positions) <- c("cellid1","cellid1_x", "cellid1_y")

  Seurat::Embeddings(default_anchors@object.list[[default_anchors@anchors$dataset2[1]]], reduction = reduction) |>
    tibble::as_tibble(rownames="cellid2") |>
    dplyr::select(1:3) -> cellid2_positions

  colnames(cellid2_positions) <- c("cellid2","cellid2_x", "cellid2_y")

  positions |>
    dplyr::left_join(cellid1_positions) |>
    dplyr::left_join(cellid2_positions) |>
    tibble::add_column(anchor_weight=1) |>
    tibble::as_tibble() -> positions

  return(positions)

}
