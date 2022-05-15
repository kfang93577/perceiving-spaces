import matplotlib.pyplot as plt
import numpy as np
from sktime.datasets import load_arrow_head
from sktime.datatypes import convert

X,y = load_arrow_head(return_X_y=True)
X = convert(X, from_type="nested_univ", to_type="numpy3D")
labels, counts = np.unique(y, return_counts=True)
fig, ax = plt.subplots(1, figsize=plt.figaspect(0.25))
for label in labels:
    ax.plot(X[y == label, 0, :][0], label = f"class {label} ")
ax.set(ylabel="Scaled distance from midpoint", xlabel = "Index")

plt.show()