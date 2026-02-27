
from ateams import Chain, Recorder, Player
from ateams.statistics import critical
from ateams.complexes import Cubical
from ateams.models import Glauber, SwendsenWang, InvadedCluster
import numpy as np


C = Cubical().fromCorners([10]*2)
model = InvadedCluster(C, field=2, dimension=1)
M = Chain(model, steps=1000)

actuals = []

with Recorder().record("data/out.lz", blocksize=23) as rec:
	for (spins, occupied, satisfied) in M.progress():
		rec.store((occupied, satisfied))


saveds = []

with Player().playback("data/out.lz", steps=1000) as play:
	for (occupied, satisfied) in play.progress():
		pass