from glob import glob
from pmonitor import PMonitor

worldcover_path = '/data/worldcover/worldcoverlandmask.zarr'
globalislands_path = '/data/globalislands/globalislandsinlandmask.zarr'
continentalshoreline_path = '/data/continentalshoreline/continentalshorelinemask.zarr'
countslist_path = '/data/sen2water/globalmaskcounts.list'
zeromask_root = 'zero-masks'
s2w_mask_root = 'global-masks'

inputs = sorted(glob(f"{zeromask_root}/*tif"))
pm = PMonitor(inputs,
              request='globalmaskproduction',
              hosts=[('localhost', 12)],
              types=[('globalmasksimple.py', 6)],
              logdir='log',
              simulation=False)
for input in inputs:
    # zero-masks/32UME-zeros.tif
    granule = input[input.rfind('/')+1:input.rfind('/')+6]
    # globalmasksimple.py /data/worldcover/worldcoverlandmask.zarr \
    #     /data/globalislands/globalislandsinlandmask.zarr \
    #     /data/continentalshoreline/continentalshorelinemask.zarr \
    #     /data/sen2water/globalmaskcounts.list \
    #     zero-masks/32UME-zeros.tif \
    #     global-masks/s2w-globalmask-32UME.tif
    pm.execute('globalmasksimple.py',
               [ input ],
               [ f"{s2w_mask_root}/s2w-globalmask-{granule}.tif" ],
               [ worldcover_path, globalislands_path, continentalshoreline_path, countslist_path ])

pm.wait_for_completion()
