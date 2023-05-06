import numpy as np
from pathlib import Path

path = str(Path("../Masters_project/").resolve())
conc_2011 = np.load(f"{path}/Data_arrays/north/ice_conc_icr_2011.npy")
conc = np.load(f"{path}/Data_arrays/north/conc_1979-2020.npy")
print(np.nanmean(conc[32, ...]))
print(np.nanmean(conc[31, ...]))
# print(f"conc_2011: {np.nanmean(conc_2011[0, 1, ...])}")
#
# print(conc.shape, conc_2011.shape)
#
# conc[32, ...] = conc_2011
#
# print(np.nanmean(conc[32, ...]))
# print(conc.shape)