
import numpy as np
from math import comb

from ..arithmetic import ComputePersistencePairs
from ..common import Matrices, FINT


class InvasionPercolation():
	_name = "InvasionPercolation"
	
	def __init__(self, C, dimension=1, initial=None, full=False, stop=lambda: 1, **kwargs):
		"""
		Initializes an invasion percolation on the given complex.

		Args:
			C (Complex): The `Complex` object on which we'll be running experiments.
			dimension (int=1): The dimension of cells on which we're percolating.
			initial (np.array): A vector of spin assignments to components.
			full (bool=True): Do we find *all* the homological persistence events,
				or just the target one?
			stop (function): A function that returns the number of essential cycles
				found before sampling the next configuration.
		"""
		# Object access.
		self.complex = C
		self.dimension = dimension
		self.stop = stop
		self._returns = 1
		self.field = 2
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
		if not initial: self.spins = self._initial()
		else: self.spins = (initial%self.field).astype(FINT)

	
	def _delegateComputation(self):
		low, high = self.complex.breaks[self.dimension], self.complex.breaks[self.dimension+1]
		times = set(range(self.cellCount))

		def whittle(pairs):
			_births, _deaths = zip(*pairs)
			births = set(_births)
			deaths = set(_deaths)
			_essential = sorted(set(
				e for e in times-(births|deaths)
				if low <= e < high
			))

			if self.full: return np.array(_essential, dtype=int)
			else: return np.array([_essential[self._STOP]], dtype=int)
		
		
		def persist(filtration):
			essential = ComputePersistencePairs(self.matrices.full, filtration, self.dimension, self.complex.breaks)
			return whittle(essential)

		self.persist = persist


	def _filtrate(self):
		"""
		Constructs a filtration based on the evaluation of the cochain.
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

	

	def _proposal(self, time):
		"""
		Proposal scheme for generalized invaded-cluster evolution on the
		random-cluster model.

		Args:
			time (int): Step in the chain.

		Returns:
			A numpy array representing a vector of spin assignments.
		"""
		# Set a stopping time.
		self._STOP = self.stop();

		# Construct the filtration and find the essential cycles.
		filtration, shuffled = self._filtrate()
		essential = self.persist(filtration)

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
			occupied = np.zeros((self.rank, self.cells), dtype=int)

			for t in sorted(essential):
				occupiedIndices = shuffled[:t-low+1]
				occupied[j,occupiedIndices] = 1
				j += 1

		return occupied
	

	def _assign(self, cocycle):
		"""
		Updates mappings from faces to spins and cubes to occupations.

		Args:
			cocycle (np.array): Cocycle on the subcomplex.
		
		Returns:
			Nothing.
		"""
		self.spins = cocycle

