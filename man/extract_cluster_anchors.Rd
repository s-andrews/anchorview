% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_cluster_anchors.R
\name{extract_cluster_anchors}
\alias{extract_cluster_anchors}
\title{Add cluster annotations to anchors}
\usage{
extract_cluster_anchors(
  anchors,
  min_score = 0,
  clusters1 = "seurat_clusters",
  clusters2 = "seurat_clusters"
)
}
\arguments{
\item{anchors}{Output of FindInteractionAnchors}

\item{min_score}{Minimum anchor score to use}

\item{clusters1}{The name of the cluster metadata in dataset1}

\item{clusters2}{The name of the cluster metadata in dataset2}
}
\value{
A tibble showing all anchors with cellids and cluster names
}
\description{
Add cluster annotations to anchors
}
