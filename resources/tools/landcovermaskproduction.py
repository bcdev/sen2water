from glob import glob
from pmonitor import PMonitor

landcover_path = '/data/fire/C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.tif'
countslist_path = '/data/fire/lccounts.list'
zeromask_root = 'zero-masks'
lc_mask_root = '/data/fire/lc-masks'

inputs = sorted(glob(f"{zeromask_root}/*tif"))
pm = PMonitor(inputs,
              request='landcovermaskproduction',
              hosts=[('localhost', 12)],
              types=[('landcovergranule.py', 6)],
              logdir='log',
              simulation=False)
for input in inputs:
    # zero-masks/32UME-zeros.tif
    granule = input[input.rfind('/')+1:input.rfind('/')+6]
    pm.execute('landcovergranule.py',
               [ input ],
               [ f"{lc_mask_root}/lc-globalmask-{granule}.tif" ],
               [ landcover_path, countslist_path ])

pm.wait_for_completion()
