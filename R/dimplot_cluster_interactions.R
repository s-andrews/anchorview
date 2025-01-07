#' Dimensionality reduction plot with anchored clusters
#'
#' @param cluster_anchors Output of FindIntegrationAnchors on clustered data
#' @param min_score Minimum anchor score to use
#' @param clusters1 Cluster metadata name in dataset 1
#' @param clusters2 Cluster metadata name in dataset 2
#' @param reduction Name of reduction to use
#' @param min_percent Minimum percent anchor linkage to show in plot
#'
#' @return ggplot of dimension plot coloured by cluster showing linkages
#' @export
#'
#' @examples
dimplot_cluster_interactions <- function(cluster_anchors, min_score=0, clusters1="seurat_clusters", clusters2="seurat_clusters", reduction="umap", min_percent=10) {

  clustered_anchors -> cluster_anchors

  # We need to get the values for the percentage overlaps

  extract_cluster_anchors(cluster_anchors, min_score=min_score, clusters1 = clusters1, clusters2=clusters2) -> extracted_clusters

  extracted_clusters |>
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
    dplyr::mutate(percent = (percent1+percent2)/2) -> percentages

  # We need all of the xy positions from the two datasets

  Seurat::Embeddings(cluster_anchors@object.list[[1]], reduction=reduction) -> d1_embeddings
  Seurat::Embeddings(cluster_anchors@object.list[[2]], reduction=reduction) -> d2_embeddings

  colnames(d1_embeddings) -> dimnames

  d1_embeddings |>
    tibble::as_tibble(rownames="cellid") |>
    dplyr::left_join(cluster_anchors@object.list[[1]][[]] |> tibble::as_tibble(rownames="cellid") |> dplyr::select(cellid,{{ clusters1 }})) |>
    dplyr::rename(cluster={{clusters1}}) |>
    dplyr::mutate(cluster=paste0("data1_",cluster)) -> d1_embeddings


  d2_embeddings |>
    tibble::as_tibble(rownames="cellid") |>
    dplyr::left_join(cluster_anchors@object.list[[2]][[]] |> tibble::as_tibble(rownames="cellid") |> dplyr::select(cellid,{{ clusters2 }})) |>
    dplyr::rename(cluster={{clusters2}}) |>
    dplyr::mutate(cluster=paste0("data2_",cluster)) -> d2_embeddings


  # We need the centroid positions for all clusters
  dplyr::bind_rows(d1_embeddings, d2_embeddings) |>
    dplyr::rename(x=dimnames[1], y=dimnames[2]) |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
      x = median(x),
      y = median(y)
    ) -> cluster_centroids

  # Now we need to add the cluster positions to the linkages we're plotting

  percentages |>
    dplyr::select(cluster1,cluster2,percent) |>
    dplyr::mutate(
      cluster1=paste0("data1_",cluster1),
      cluster2=paste0("data2_",cluster2)
    ) |>
    dplyr::left_join(
      cluster_centroids |> dplyr::rename(cluster1=cluster,xstart=x,ystart=y)
    ) |>
    dplyr::left_join(
      cluster_centroids |> dplyr::rename(cluster2=cluster,xend=x,yend=y)
    ) |>
    dplyr::filter(percent > min_percent) -> linkages



  dplyr::bind_rows(d1_embeddings, d2_embeddings) |>
    ggplot2::ggplot(ggplot2::aes(x=.data[[dimnames[1]]], y=.data[[dimnames[2]]], colour=cluster)) +
    ggplot2::geom_point(size=0.5) +
    ggplot2::annotate(
      geom="segment",
      x=linkages$xstart,
      y=linkages$ystart,
      xend=linkages$xend,
      yend=linkages$yend,
      colour="black",
      alpha=0.5,
      size=linkages$percent/30
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=5))) -> plot


  return(plot)

}
