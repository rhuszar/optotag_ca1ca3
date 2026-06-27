"""
Demixed Principal Component Analysis (dPCA), generalized to an arbitrary
number of task marginalization axes.

Reference implementation of the algorithm of Kobak et al. (2016), eLife,
"Demixed principal component analysis of neural population data".

Input is a matrix of binned, smoothed firing rates together with a set of
integer condition labels (one per axis), supplied as plain numpy arrays.
"""
import numpy as np
from itertools import combinations, product as iter_product


class dPCA:
    def __init__(self, n_components=10, lambda_reg=1e-5):
        self.default_n = n_components if isinstance(n_components, int) else 10
        self.n_components_dict = n_components if isinstance(n_components, dict) else {}
        self.lambda_reg = lambda_reg
        self.mean_ = None
        self.encoders = {}
        self.decoders = {}
        self.axis_names = []

    def _compute_stats(self, X, bins_dict):
        n_neurons = X.shape[1]
        axis_names = list(bins_dict.keys())

        indices_list = []
        n_levels = []
        for name in axis_names:
            u, inv = np.unique(bins_dict[name], return_inverse=True)
            indices_list.append(inv)
            n_levels.append(len(u))

        PSTH_tensor = np.zeros((n_neurons,) + tuple(n_levels))
        noise_cov_accum = np.zeros((n_neurons, n_neurons))
        valid_noise_conds = 0

        for idx in iter_product(*[range(n) for n in n_levels]):
            mask = np.ones(X.shape[0], dtype=bool)
            for ax_i, level in enumerate(idx):
                mask &= (indices_list[ax_i] == level)

            X_cond = X[mask]
            k = len(X_cond)

            if k > 0:
                cond_mean = np.mean(X_cond, axis=0)
                PSTH_tensor[(slice(None),) + idx] = cond_mean

                if k > 1:
                    X_centered = X_cond - cond_mean[None, :]
                    cov_cond = (X_centered.T @ X_centered) / k
                    noise_cov_accum += cov_cond
                    valid_noise_conds += 1

        avg_noise_cov = np.eye(n_neurons)
        if valid_noise_conds > 0:
            avg_noise_cov = noise_cov_accum / valid_noise_conds

        return PSTH_tensor, avg_noise_cov

    def _compute_marginalizations(self, PSTH_centered, axis_names):
        K = len(axis_names)
        num_neurons = PSTH_centered.shape[0]

        all_subsets = []
        for r in range(1, K + 1):
            for combo in combinations(range(K), r):
                all_subsets.append(frozenset(combo))

        raw_marginals = {}
        for subset in all_subsets:
            axes_to_avg = tuple(i + 1 for i in range(K) if i not in subset)
            if axes_to_avg:
                marginal = np.mean(PSTH_centered, axis=axes_to_avg, keepdims=True)
                marginal = np.broadcast_to(marginal, PSTH_centered.shape).copy()
            else:
                marginal = PSTH_centered.copy()
            raw_marginals[subset] = marginal

        interactions = {}
        for subset in sorted(all_subsets, key=len):
            interaction = raw_marginals[subset].copy()
            for other in interactions:
                if other < subset:
                    interaction -= interactions[other]
            interactions[subset] = interaction

        targets = {}
        for subset, tensor in interactions.items():
            label = ''.join(axis_names[i] for i in sorted(subset))
            targets[label] = tensor.reshape(num_neurons, -1)

        return targets

    def fit(self, X, bins_dict_or_s_bins=None, t_bins=None, use_noise_cov=True):
        """
        X : (n_samples, n_neurons) binned, smoothed firing rates.
        Condition labels: either a dict {axis_name: (n_samples,) int bins}, e.g.
        {'s': space_bins, 't': time_bins, 'r': reward_bins}, or the legacy
        positional form fit(X, s_bins, t_bins).
        Sets self.encoders / self.decoders, keyed by axis-subset labels.
        """
        if isinstance(bins_dict_or_s_bins, dict):
            bins_dict = dict(bins_dict_or_s_bins)
        elif bins_dict_or_s_bins is not None and t_bins is not None:
            bins_dict = {'s': bins_dict_or_s_bins, 't': t_bins}
        else:
            raise ValueError("Provide either a bins_dict or both s_bins and t_bins")

        self.axis_names = list(bins_dict.keys())

        PSTH_tensor, noise_cov = self._compute_stats(X, bins_dict)
        num_neurons = PSTH_tensor.shape[0]

        avg_axes = tuple(range(1, PSTH_tensor.ndim))
        self.mean_ = np.mean(PSTH_tensor, axis=avg_axes)
        mean_shape = (num_neurons,) + (1,) * (PSTH_tensor.ndim - 1)
        PSTH_centered = PSTH_tensor - self.mean_.reshape(mean_shape)

        targets = self._compute_marginalizations(PSTH_centered, self.axis_names)

        X_flat = PSTH_centered.reshape(num_neurons, -1)

        cov_signal = X_flat @ X_flat.T
        cov_total = cov_signal + noise_cov if use_noise_cov else cov_signal
        inverse_term = np.linalg.pinv(cov_total + self.lambda_reg * np.eye(num_neurons))

        self.encoders = {}
        self.decoders = {}

        for key, target_matrix in targets.items():
            cross_cov = target_matrix @ X_flat.T
            beta_ols = cross_cov @ inverse_term
            y_hat = beta_ols @ X_flat
            U, S, Vt = np.linalg.svd(y_hat, full_matrices=False)

            k = self.n_components_dict.get(key, self.default_n)
            encoder = U[:, :k]
            decoder = encoder.T @ beta_ols

            self.encoders[key] = encoder
            self.decoders[key] = decoder

        return self

    def reconstruct(self, X, decoder_overrides=None):
        if decoder_overrides is None:
            decoder_overrides = {}

        X_centered = X - self.mean_
        reconstruction = np.zeros_like(X_centered)

        for axis in self.encoders:
            if axis in decoder_overrides:
                enc = decoder_overrides[axis]['encoder']
                dec = decoder_overrides[axis]['decoder']
            else:
                enc = self.encoders[axis]
                dec = self.decoders[axis]

            latent = X_centered @ dec.T
            reconstruction += latent @ enc.T

        return reconstruction + self.mean_

    def embed(self, X, axis='t'):
        X_centered = X - self.mean_
        return X_centered @ self.decoders[axis].T
