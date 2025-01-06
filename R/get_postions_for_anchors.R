

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
    dplyr::mutate(cellid1_x = Seurat::Embeddings(anchors@object.list[[dataset1[1]]], reduction = reduction)[cellid1,] |> tibble::as_tibble() |> dplyr::pull(1) )|>
    dplyr::mutate(cellid1_y = Seurat::Embeddings(anchors@object.list[[dataset1[1]]], reduction = reduction)[cellid1,] |> tibble::as_tibble() |> dplyr::pull(2) )|>
    dplyr::mutate(cellid2_x = Seurat::Embeddings(anchors@object.list[[dataset2[1]]], reduction = reduction)[cellid2,] |> tibble::as_tibble() |> dplyr::pull(1) )|>
    dplyr::mutate(cellid2_y = Seurat::Embeddings(anchors@object.list[[dataset2[1]]], reduction = reduction)[cellid2,] |> tibble::as_tibble() |> dplyr::pull(2) )|>
    dplyr::ungroup() -> positions


  positions |>
    tibble::add_column(anchor_weight=1) |>
    tibble::as_tibble() -> positions

  return(positions)

}
