#' Add cluster annotations to anchors
#'
#' @param anchors Output of FindInteractionAnchors
#' @param clusters1 The name of the cluster metadata in dataset1
#' @param clusters2 The name of the cluster metadata in dataset2
#' @param min_score Minimum anchor score to use
#'
#' @return A tibble showing all anchors with cellids and cluster names
#' @export
#'
#' @examples
extract_cluster_anchors <- function(anchors, min_score=0, clusters1="seurat_clusters", clusters2="seurat_clusters") {

  # We need to get a data structure which has dataset1 cellid, dataset1 cluster, dataset2 cellid, dataset2 cluster

  anchors@anchors |>
    dplyr::filter(score>=min_score) |>
    dplyr::filter(dataset1==1) |>
    dplyr::mutate(cellid1 = colnames(anchors@object.list[[1]])[cell1]) |>
    dplyr::mutate(cluster1 = anchors@object.list[[1]][[]][cellid1,clusters1]) |>
    dplyr::mutate(cellid2 = colnames(anchors@object.list[[2]])[cell2]) |>
    dplyr::mutate(cluster2 = anchors@object.list[[2]][[]][cellid2,clusters2]) |>
    tibble::as_tibble() -> clusters12

  anchors@anchors |>
    dplyr::filter(score>=min_score) |>
    dplyr::filter(dataset1==2) |>
    rename(cell2=cell1,cell1=cell2,dataset1=dataset2,dataset2=dataset1) |>
    dplyr::mutate(cellid1 = colnames(anchors@object.list[[1]])[cell1]) |>
    dplyr::mutate(cluster1 = anchors@object.list[[1]][[]][cellid1,clusters1]) |>
    dplyr::mutate(cellid2 = colnames(anchors@object.list[[2]])[cell2]) |>
    dplyr::mutate(cluster2 = anchors@object.list[[2]][[]][cellid2,clusters2]) |>
    tibble::as_tibble() -> clusters21

  dplyr::bind_rows(clusters12,clusters21) -> cluster_interactions

  return(cluster_interactions)

}
