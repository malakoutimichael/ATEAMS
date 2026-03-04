
import numpy as np
from math import comb

from ..arithmetic import ComputePersistencePairs
from ..common import Matrices, FINT


class InvasionPercolation():
	_name = "InvasionPercolation"
	
	def __init__(self, C, dimension=1, initial=None, full=True, stop=lambda: 1, **kwargs):
		"""
		Initializes an invasion percolation on the given complex.

		Args:
			C (Complex): The `Complex` object on which we'll be running experiments.
			dimension (int=1): The dimension of cells on which we're percolating.
			initial (np.array): A boolean vector indicating which cells are included
				in the initial subcomplex.
			full (bool=True): Do we find *all* the giant cycles?
			stop (function): A function that returns the number of giant cycles
				found before sampling the next configuration.
		"""
		# Object access.
		self.complex = C
		self.dimension = dimension
		self.stop = stop
		self._returns = 4
		self.full = full

		# Check if we're debugging.
		self._DEBUG = kwargs.get("_DEBUG", False)

		# Useful values to have later.
		self.cellCount = len(self.complex.flattened)
		self.cells = len(self.complex.Boundary[self.dimension])
		self.faces = len(self.complex.Boundary[self.dimension-1])
		self.target = np.arange(self.complex.breaks[self.dimension], self.complex.breaks[self.dimension+1])

		# Premake the "occupied cells" array; change the dimension of the complex
		# to correspond to the provided dimension.
		self.rank = comb(len(self.complex.corners), self.dimension)

		# Force-recompute the matrices for a different dimension; creates
		# a set of orientations for fast elementwise products.
		self.matrices = Matrices()
		self.matrices.full = self.complex.matrices.full

		# Delegates computation for persistence and cocycle sampling.
		self._delegateComputation()

		# Seed the random number generator.
		self.RNG = np.random.default_rng()

		# If no initial spin configuration is passed, initialize.
		self.spins = self._initial() if initial is None else initial

	
	def _delegateComputation(self):
		low, high = self.complex.breaks[self.dimension], self.complex.breaks[self.dimension+1]
		times = set(range(self.cellCount))
		
		
		def whittle(births, deaths):
			# Avoid doing this computation twice; TODO come back and fix this later
			_essential = sorted(set(
				e for e in times-(births|deaths)
				if low <= e < high
			))

			if self.full: return np.array(_essential, dtype=int)
			else: return np.array([_essential[self._STOP]], dtype=int)

		
		def persist(filtration):
			pairs = ComputePersistencePairs(self.matrices.full, filtration, self.dimension, self.complex.breaks)

			_births, _deaths = zip(*pairs)
			births = set(_births)
			deaths = set(_deaths)

			total = np.array(sorted(times-(births|deaths)), dtype=int)
			essential = whittle(births, deaths)
			pairs = np.array(pairs).T

			return essential, total, pairs
		
		self.persist = persist


	def _filtrate(self):
		"""
		Constructs a filtration.
		"""
		# Construct the filtration. We could probably get away with having a
		# fixed `self.filtration` instance property and then just modifying it
		# at each step; that might be a later thing, though.
		filtration = np.arange(self.cellCount)
		low = self.complex.breaks[self.dimension]
		high = self.complex.breaks[self.dimension+1]

		indices = np.arange(self.cells)
		shuffled = np.random.permutation(indices)
		filtration[low:high] = self.target[shuffled]

		return filtration, shuffled


	def _initial(self):
		"""
		Computes an initial state for the model's Complex.

		Returns:
			A boolean vector representing a random initial subcomplex.
		"""
		init = self.RNG.uniform(size=self.cells)
		return (init < 1/2).astype(bool)

	

	def proposal(self, time):
		"""
		Proposal scheme for generalized invasion percolation evolution.

		Args:
			time (int): Step in the chain.

		Returns:
			A 4-tuple:

			1. a boolean array where each column corresponds to a \(d\)-cell,
				and each row corresponds to a \(d\)-dimensional homological
				percolation event, so each entry indicates the presence or
				absence of that \(d\)-cell when the \(d\)-dimensional homological
				percolation event occurred;
			2. times at which \(d\)-dimensional homological percolation events
				occurred;
			3. times at which *all* homological percolation events occurred;
			4. a two-dimensional int array encoding persistence pairs.

		"""
		# Set a stopping time.
		self._STOP = self.stop();

		# Construct the filtration and find the essential cycles.
		filtration, shuffled = self._filtrate()
		essential, total, pairs = self.persist(filtration)

		j = 0
		low = self.complex.breaks[self.dimension]

		# If we only computed a desired percolation time, then the matrix of
		# occupied cell indices has one row...
		if not self.full:
			t = essential[0]
			occupied = np.zeros(self.cells, dtype=int)
			occupiedIndices = shuffled[:t-low+1]
			occupied[occupiedIndices] = 1

		# Otherwise, it has `rank` rows.
		else:
			occupied = np.zeros((self.rank, self.cells), dtype=bool)

			for t in essential:
				occupiedIndices = shuffled[:t-low+1]
				occupied[j,occupiedIndices] = 1
				j += 1

		return occupied, essential, total, pairs
	

	def _assign(self, cochain): pass

