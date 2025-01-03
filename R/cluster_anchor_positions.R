#' Summarise integration anchors
#'
#' @param anchor_positions A raw set of anchor positions from get_positions_for_anchors
#' @param min_score The minimum score for the anchors to use
#' @param cluster_threshold The euclidean clustering threshold to use - higher gives fewer clusters
#'
#' @return A tibble with grouped and summarised anchor positions
#' @export
#'
#' @examples
cluster_anchor_positions <- function(anchor_positions, min_score=0, max_score=1, cluster_threshold=5) {

  anchor_positions |>
    dplyr::filter(score>=min_score & score <= max_score) -> anchor_positions


  hclust(dist(anchor_positions %>% select(cellid1_x,cellid1_y, cellid2_x, cellid2_y), method = "euclidean")) -> clustered_anchors

  cutree(clustered_anchors,h = cluster_threshold) -> groups

  anchor_positions %>%
    tibble::add_column(cluster_groups = groups) |>
    dplyr::group_by(cluster_groups) |>
    dplyr::summarise(
      cellid1_x = median(cellid1_x, na.rm = TRUE),
      cellid1_y = median(cellid1_y, na.rm = TRUE),
      cellid2_x = median(cellid2_x, na.rm = TRUE),
      cellid2_y = median(cellid2_y, na.rm = TRUE),
      score=mean(score),
      anchor_weight = sum(anchor_weight)
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(anchor_weight) -> weighted_anchor_positions

  return(weighted_anchor_positions)

}
