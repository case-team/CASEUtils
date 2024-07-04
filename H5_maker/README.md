# H5 Maker (for semi-lep ttbar)

Make h5's from the PFNano's after applying a preselection. 

To run locally, try something like `python make_h5_local.py --ttbar -i input_file.root -o output_file.h5 -y 201X -f 0 `

The `--ttbar` option tells the processor to use the selection for semi-leptonic ttbar.
The option `-f` sets the truth label of the output. Usually signal is 1, QCD is 0, single top -1, ttbar -2, V+jets -3
Use the `--sys` option if you are running on a MC and want to save the event weights for the systematic variation. 
The options are described a bit in `make_h5_local.py`


## Output
The output is an h5 dataset with several different keys

(Documentation should be updated)




## H5_merge
H5\_merge.py is also useful. It combines different h5 files together like hadd.
Does a weighted average for preselection\_eff
