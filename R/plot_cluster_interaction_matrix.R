#' Plot a cluster level interaction matrix for anchors
#'
#' @param cluster_anchors
#'
#' @return heatmap of cluster percentage interaction
#' @export
#'
#' @examples
plot_cluster_interaction_matrix <- function(cluster_anchors) {

  cluster_anchors |>
    dplyr::group_by(cluster1,cluster2) |>
    dplyr::count() |>
    dplyr::ungroup() |>
    tidyr::complete(cluster1,cluster2) |>
    dplyr::mutate(n=tidyr::replace_na(n,0)) |>
    dplyr::group_by(cluster1) |>
    dplyr::mutate(percent1 = 100*n/sum(n)) |>
    dplyr::group_by(cluster2) |>
    dplyr::mutate(percent2 = 100*n/sum(n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(percent = (percent1+percent2)/2) |>
    tidyheatmaps::tidy_heatmap(
      rows=cluster1,
      columns=cluster2,
      values=percent,
      cluster_rows=TRUE,
      cluster_cols = TRUE,
      scale="none",
      colors = c("white","red2"),
      border_color = "grey"
    ) -> plot

  return(plot)

}
