# Performing the trigger efficiency study

## Producing trigger ntuples

```
python make_trigger_h5_condor.py -q longlunch -y 2016 -o [directory to store ntuples] -n [NanoAOD tools src dir] -c [path to GRID certificate] 
python make_trigger_h5_condor.py -q longlunch -y 2017 -o [directory to store ntuples] -n [NanoAOD tools src dir] -c [path to GRID certificate]
python make_trigger_h5_condor.py -q longlunch -y 2018 -o [directory to store ntuples] -n [NanoAOD tools src dir] -c [path to GRID certificate]
```

## Drawing trigger efficiency vs mjj

```
python plot_trigger_efficiency.py -i [directory to store ntuples] -o trigger_efficiency_run2.pdf -y run2 --eta 2.5 --pt 300 --deta 1.3 --m 0 --combined
```

## Drawing 2D trigger efficiency vs mjj and jet mass

```
python plot_trigger_efficiency_2D.py -i [directory to store ntuples] -o 2D_trigger_efficiency_2016_mjj_m1.pdf -y 2016 --eta 2.5 --pt 300 --deta 1.3 --x_variable mjj --y_variable m1
python plot_trigger_efficiency_2D.py -i [directory to store ntuples] -o 2D_trigger_efficiency_2016_mjj_m2.pdf -y 2016 --eta 2.5 --pt 300 --deta 1.3 --x_variable mjj --y_variable m2
python plot_trigger_efficiency_2D.py -i [directory to store ntuples] -o 2D_trigger_efficiency_2017_mjj_m1.pdf -y 2017 --eta 2.5 --pt 300 --deta 1.3 --x_variable mjj --y_variable m1
python plot_trigger_efficiency_2D.py -i [directory to store ntuples] -o 2D_trigger_efficiency_2017_mjj_m2.pdf -y 2017 --eta 2.5 --pt 300 --deta 1.3 --x_variable mjj --y_variable m2
python plot_trigger_efficiency_2D.py -i [directory to store ntuples] -o 2D_trigger_efficiency_2018_mjj_m1.pdf -y 2018 --eta 2.5 --pt 300 --deta 1.3 --x_variable mjj --y_variable m1
python plot_trigger_efficiency_2D.py -i [directory to store ntuples] -o 2D_trigger_efficiency_2018_mjj_m2.pdf -y 2018 --eta 2.5 --pt 300 --deta 1.3 --x_variable mjj --y_variable m2
```
