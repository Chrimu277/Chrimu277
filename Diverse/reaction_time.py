# 4th of April 2021
# Christoph Mühleder
# DONE

import numpy as np
import matplotlib.pyplot as plt
import time
from datetime import datetime as dt

N = 3                     # measure 10 times
reactiontimes = []
np.random.seed(dt.now().microsecond)

input("Press Enter when you are ready.\n"
      "Please Press enter QUICKLY when a new HEY appears!")

for i in range(N):
    sec_wait = 5 * np.random.rand()
    time.sleep(sec_wait)
    print("\nHEY Nr #",i)
    t1 = dt.now()
    input("")
    t2 = dt.now()
    t_diff = t2-t1
    reactiontimes.append(t_diff.total_seconds())

print(reactiontimes)

# plot reaction times
plt.figure()
plt.plot(np.arange(N),reactiontimes)
plt.ylim((0,2))
plt.xlabel("Hey #")
plt.ylabel("reaction time / s")
plt.show()