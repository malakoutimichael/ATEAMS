
import numpy as np
from math import comb

from ..arithmetic import ComputePersistencePairs
from ..common import Matrices, FINT


class Bernoulli():
	_name = "Bernoulli"
	
	def __init__(self, C, p=1/2, dimension=1, **kwargs):
		"""
		Initializes classical Bernoulli percolation on the provided complex,
		detecting percolation in the `dimension`-1th homology group.

		Args:
			C (Complex): The `Complex` object on which we'll be running experiments.
			p (float=1/2): The probability with which cells of the desired dimension
				are added.
			dimension (int=1): The dimension of cells on which we're percolating.
		"""
		# Object access. Set a field.
		self.complex = C
		self.dimension = dimension
		self._returns = 2
		self.field = 2
		self.p = p

		# Set a phantom spins attribute so we don't break the Chain
		self.spins = None

		# Force-recompute the matrices for a different dimension; creates
		# a set of orientations for fast elementwise products.
		self.matrices = Matrices()
		self.matrices.full = self.complex.matrices.full

		boundary, coboundary = self.complex.recomputeBoundaryMatrices(dimension)
		self.matrices.boundary = boundary
		self.matrices.coboundary = coboundary

		# Useful values to have later.
		self.cellCount = len(self.complex.flattened)
		self.cells = len(self.complex.Boundary[self.dimension])
		self.faces = len(self.complex.Boundary[self.dimension-1])
		self.target = np.arange(self.complex.breaks[self.dimension], self.complex.breaks[self.dimension+1])

		# Premake the "occupied cells" array; change the dimension of the complex
		# to correspond to the provided dimension.
		self.rank = comb(len(self.complex.corners), self.dimension)
		self.nullity = len(self.complex.Boundary[self.dimension])

		# Set low and high markers for adding cells.
		self.low = self.complex.breaks[self.dimension]
		self.high = self.complex.breaks[self.dimension+1]

		# Delegates computation for persistence.
		self._delegateComputation()

		# Seed the random number generator.
		self.RNG = np.random.default_rng()

	
	def _delegateComputation(self):
		times = set(range(self.cellCount))
		
		def whittle(births, deaths):
			# Avoid doing this computation twice; TODO come back and fix this later
			_essential = sorted(set(
				e for e in times-(births|deaths)
				if self.low <= e < self.high
			))

			return np.array(_essential, dtype=int)

		
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
		uniform = self.RNG.uniform(size=self.cells)
		include = np.nonzero(uniform < self.p)[0]
		exclude = np.nonzero(~(uniform < self.p))[0]
		m = include.shape[0]

		filtration = np.arange(self.cellCount)

		filtration[self.low:self.low+m] = self.target[include]
		filtration[self.low+m:self.high] = self.target[exclude]

		return filtration, include


	def _initial(self):
		"""
		Computes an initial state for the model's Complex.

		Returns:
			A numpy `np.array` representing a vector of spin assignments.
		"""
		return self.RNG.integers(
			0, high=self.field, dtype=FINT, size=self.faces
		)
	

	def proposal(self, time):
		"""
		Proposal scheme for generalized Bernoulli percolation.

		Args:
			time (int): Step in the chain.

		Returns:
			A 2-tuple:

			1. a boolean array with a \(1\) in each entry corresponding to an
				included \(d\)-cell;
			2. times at which \(d\)-dimensional homological percolation events
				occurred.
		"""
		# Construct the filtration and find the essential cycles.
		filtration, included = self._filtrate()
		essential, __, _ = self.persist(filtration)

		# Find when we stopped adding cells.
		stop = self.low + included.shape[0]
		giants = essential[(self.low < essential) & (essential < stop)]
		occupied = np.zeros(self.nullity, dtype=bool)
		occupied[included] = 1

		return occupied, giants
	

	def _assign(self, cochain): pass

