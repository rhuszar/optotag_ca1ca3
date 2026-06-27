"""
Identification of stereotypical trajectories by Frechet-distance clustering.

Each trajectory is a sequence of 2D positions. Pairwise discrete Frechet
distances are clustered hierarchically, and the indices belonging to the single
dominant (most populated) cluster are returned as the stereotypical set.
"""
import collections
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from similaritymeasures import frechet_dist


def frechet_distance_matrix(trajectories):
    """trajectories: list of (n_points_i, 2) arrays. Returns (n, n) distance matrix."""
    n = len(trajectories)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = frechet_dist(trajectories[i], trajectories[j])
            D[j, i] = D[i, j]
    return D


def cluster_trajectories(distance_matrix, distance_threshold=20, method='complete'):
    """
    distance_matrix: (n, n) symmetric. Returns (labels, Z) where labels are
    contiguous cluster ids ordered by dendrogram leaves (singletons labelled -1),
    and Z is the scipy linkage matrix.
    """
    Z = linkage(distance_matrix[np.triu_indices(len(distance_matrix), 1)],
                method=method, metric='precomputed')
    cluster_ids = fcluster(Z, distance_threshold, criterion='distance')
    cluster_counter = collections.Counter(cluster_ids)
    cluster_labels = -np.ones_like(cluster_ids)

    working_label, previous_cluster = 0, None
    for leaf in leaves_list(Z):
        if cluster_counter[cluster_ids[leaf]] <= 1:
            continue
        if previous_cluster is None:
            previous_cluster = cluster_ids[leaf]
            cluster_labels[leaf] = working_label
            continue
        if cluster_ids[leaf] != previous_cluster:
            working_label += 1
        cluster_labels[leaf] = working_label
        previous_cluster = cluster_ids[leaf]

    return cluster_labels, Z


def most_common_cluster_indices(cluster_labels, min_cluster_size=5):
    """Indices of the most populated non-singleton cluster, or [] if it is smaller than min_cluster_size."""
    mask = cluster_labels >= 0
    if not np.any(mask):
        return np.array([], dtype=int)
    selected = collections.Counter(cluster_labels[mask]).most_common(1)[0][0]
    idx = np.arange(len(cluster_labels))[cluster_labels == selected]
    return idx if len(idx) >= min_cluster_size else np.array([], dtype=int)


def select_stereotypical_trajectories(trajectories, distance_threshold=20,
                                      min_cluster_size=5, subsample=5):
    """
    End-to-end selection. trajectories: list of (n_points_i, 2) arrays.
    Returns the indices (into `trajectories`) forming the dominant stereotypical cluster.
    """
    sub = [np.asarray(t)[::subsample] for t in trajectories]
    D = frechet_distance_matrix(sub)
    labels, _ = cluster_trajectories(D, distance_threshold)
    return most_common_cluster_indices(labels, min_cluster_size)
