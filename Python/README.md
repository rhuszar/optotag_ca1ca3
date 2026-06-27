Core Python routines for the following analyses:

- `dpca.py` — demixed PCA (Kobak et al., 2016, eLife), generalized to an
  arbitrary number of task marginalization axes.
- `trajectory_clustering.py` — identification of stereotypical trajectories by
  Fréchet-distance hierarchical clustering.


## dPCA usage

`X` is a `(n_samples, n_neurons)` matrix of binned, smoothed firing rates.
Condition labels are integer arrays of length `n_samples`, one per axis.

```python
import numpy as np
from dpca import dPCA

# Number of components to keep per marginalization (axis subset).
n_components = {'t': 1, 's': 4, 'r': 1, 'st': 3, 'sr': 2, 'tr': 2, 'str': 2}
model = dPCA(n_components=n_components)
model.fit(X, {'s': space_bins, 't': time_bins, 'r': reward_bins}, use_noise_cov=True)

# Encoder / decoder axes, keyed by axis-subset label ('s', 't', 'r', 'st', ...).
spatial_axis = model.encoders['s']        # (n_neurons, n_components['s'])
drift_expression = (X - X.mean(0)) @ model.encoders['t']
```

## Trajectory clustering usage

`trajectories` is a list of `(n_points_i, 2)` position arrays (one per trial).

```python
from trajectory_clustering import select_stereotypical_trajectories

idx = select_stereotypical_trajectories(
    trajectories,
    distance_threshold=20,   # Fréchet linkage threshold (cm)
    min_cluster_size=5,
    subsample=5,
)
# `idx` are indices into `trajectories` forming the dominant stereotypical cluster.
```
