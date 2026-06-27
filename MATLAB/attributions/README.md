Related to Figure 6

Assesses GPR position decoders during the learning session from two gcPCA-defined subspaces: one capturing rest2-specific activity patterns (new) and one capturing rest1 patterns (old). 

Attribution scores are then computed per neuron per timepoint to quantify each cell's contribution to the decoder's position prediction. 

Because this is computationally expensive, attributions are batched across cluster jobs and reassembled into a single timeseries, yielding a continuous record of how much each neuron drives position coding in the new vs. old subspace across the session.
