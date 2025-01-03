#' Plot integration anchors onto a dimensionality reduction plot
#'
#' @param split_objects The split Seurat object list
#' @param anchor_positions The output of get_positions_for_anchors
#' @param reduction The name of the dimension reduction to use
#' @param min_score The minimum anchor score to use
#'
#' @return A ggplot showing the anchors overlaid on the dimension plot
#' @export
#'
#' @examples
plot_integration_anchors <- function(split_objects, anchor_positions, reduction="umap", min_score=0) {

  anchor_positions |>
    dplyr::filter(score >= min_score) -> anchor_positions

  names(split_objects) -> conditions

  Seurat::Embeddings(split_objects[[1]], reduction="umap.unintegrated") |>
    tibble::as_tibble() |>
    tibble::add_column(
      sample=conditions[1]
    ) -> dimension_reduction1

  Seurat::Embeddings(split_objects[[2]], reduction="umap.unintegrated") |>
    tibble::as_tibble() |>
    tibble::add_column(
      sample=conditions[2]
    ) -> dimension_reduction2


  dplyr::bind_rows(dimension_reduction1,dimension_reduction2) -> plot_data

  expansion_factor <- dplyr::if_else(max(anchor_positions$anchor_weight) == 1, 0.5, 6)

  plot_data |>
    ggplot2::ggplot(ggplot2::aes(x=.data[[(colnames(plot_data)[1])]], y=.data[[(colnames(plot_data)[2])]], colour=sample)) +
    ggplot2::geom_point(size=0.5) +
    ggplot2::annotate(
      geom="segment",
      x=anchor_positions$cellid1_x,
      y=anchor_positions$cellid1_y,
      xend=anchor_positions$cellid2_x,
      yend=anchor_positions$cellid2_y,
      colour="black",
      size=anchor_positions$anchor_weight/max(anchor_positions$anchor_weight)*expansion_factor,
      alpha=0.5
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=5)))-> anchor_plot


  return(anchor_plot)
}
