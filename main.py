import numpy as np
from consts import * #
from fuselage import * #imports fuselage calc function
import math
import matplotlib.pyplot as plt
from intersect import intersection
import inspect

plt.figure(figsize=(15,5))
run = designRun(0.0168)
run.runthatshit()
run.showDashboard()

run2 = designRun(0.02)
run2.runthatshit()
run2.showDashboard()