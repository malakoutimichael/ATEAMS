
import numpy as np


def bettis(pairs, essential, breaks):
	r"""
	Computes the betti curves admitted by the provided persistence pairs. The
	ith betti curve is the piecewise function that records the rank of the ith
	homology group of the filtration over time.

	Args:
		pairs (np.array): Two-dimensional integer array (e.g. the one returned by
			`ateams.models.InvasionPercolation.proposal`) representing persistence pairs.
		esssential (np.array): Integer array of all homological percolation times
			(e.g. the one returned by `ateams.models.InvasionPercolation.proposal`).
		breaks (np.array): Two-dimensional integer array where the \(i\)th row
			indicates the first and last indices of \(i\)-cells in the complete
			complex filtration. Can be constructed from the `.breaks` property
			of any `Model`, or given immediately by the `.tranches` property of
			any `Model`.

	Returns:
		A \(k \times N\) array where the \(i\)th row corresponds to the expected
		degree of the \(i\)th homology group at all times \(0 \leq t \leq N\), for
		\(N\) the total number of cells in the complex.
	"""
	# The ith Betti curve of the persistence diagram specified by the
	# persistence pairs is the function that records the rank of the ith
	# homology group of the filtration over time. Since we have to store
	# these values contiguously (somehow), we need to flatten them into
	# two objects: a 2d array of all points on all betti curves, and a
	# 1d array of column indices separating them.

	curves = np.zeros((breaks.shape[0], essential.max()+1), dtype=int)
	indices = np.zeros(breaks.shape, dtype=int)

	for d in range(breaks.shape[0]-1):
		lo = breaks[d,0];
		hi = breaks[d+1,1];
		indices[d,0] = lo;
		indices[d,1] = hi;

	indices[indices.shape[0]-1] = breaks[breaks.shape[0]-1]

	for dimension in range(breaks.shape[0]):
		lo, hi = breaks[dimension];

		subpairs = pairs[:,(lo <= pairs[0]) & (pairs[0] < hi)];
		ess = essential[(lo <= essential) & (essential < hi)];

		# Find the birth times (for this dimension)...
		btimes = np.concatenate([subpairs[0], ess])
		btimes.sort()

		# ... and the death times (for this dimension)...
		dtimes = subpairs[1]
		dtimes.sort()

		# ... then indicate the points defining the piecewise function.
		times = np.concatenate([btimes, dtimes])
		ranks = np.zeros(btimes.shape[0]+dtimes.shape[0], dtype=int)

		ranks[:btimes.shape[0]] = np.arange(1,btimes.shape[0]+1)
		ranks[btimes.shape[0]:] = np.flip(np.arange(btimes.shape[0]-dtimes.shape[0], btimes.shape[0]))

		curves[dimension,times] = ranks
		curves[dimension,times.max():] = ess.shape[0]

		# Fill in missing data points.
		lo, hi = indices[dimension];

		for t in range(lo, hi):
			r = curves[dimension,t]
			
			if r < 1: curves[dimension,t] = last
			if r > 0: last = curves[dimension,t]

	
	return curves

