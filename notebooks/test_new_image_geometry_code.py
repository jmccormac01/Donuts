"""
Quickly test the new updates with some
NGTS images
"""
import os
import glob as g
from donuts import Donuts

os.chdir('/Users/jmcc/Dropbox/data/ngts/action293930_observeField')
imgs = sorted(g.glob("*.fits"))

ref = Donuts(imgs[0], exposure="EXPOSURE", prescan_width=20,
             overscan_width=20, scan_direction='x', calculation_area_override=[100,1000,100,1000],
             ntiles=10)

for i in imgs[1:]:
    shift = ref.measure_shift(i)
    print(shift.x, shift.y)
