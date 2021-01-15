#%%
from typing import List
from numpy.lib.function_base import iterable
import uproot
import numpy as np
import pandas as pd
import awkward as ak
import matplotlib.pyplot as plt

fn = "~/CMSSW_11_1_6/src/ntupleTree-re.root"
rf = uproot.open(fn)
#%%
for key in rf["/treeMaker;1/tree;1"].keys():
    print(key)
    a = rf["/treeMaker;1/tree;1/" + key].array()
    print(a)
    print(ak.type(a))
    print()

# %%
def ga(vars: List[str], eventNumber: int = 0, library: str = "np") -> np.ndarray:
    valvars = rf["/treeMaker;1/tree;1"]
    outpd = pd.DataFrame({var: valvars[var].array(library=library)[eventNumber] for var in vars })

    return outpd 
#%%
arr=ga(['recHit_y','recHit_x','recHit_z','recHit_layer'])
fig1, ((a1,a2),(a3,a4)) = plt.subplots(ncols=2, nrows=2)
axes=[a1,a2,a3,a4]
print(type(f1_axes[0]))

for layer in range(4):
    sel=arr[arr['recHit_layer']==layer]
    axes[layer].scatter(arr['recHit_x'],arr['recHit_y'])
    axes[layer].set_title(f'Layer {layer}')



# %%

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(xs=ga("recHit_y"), ys=ga("recHit_z"), zs=ga("recHit_x"))
ax.set_xlabel("X Label")
ax.set_ylabel("Y Label")
ax.set_zlabel("Z Label")
plt.savefig("firstLight.png")
# %%
plt.hist(ga("recHit_z"))
# %%
