cd /d %~dp0
mpiexec -np 11 python spa_mpirun.py --db synapt_example_database.sqlite3 synapt_example_inputs/synapt_example.txt

python spa_db_fit.py --db synapt_example_database.sqlite3 synapt_example_inputs/synapt_example.txt --tpoints 1000 --outfile synapt_outputs.txt --plotfile screen --plotmap screen --weighting combined

python spa_db_plot.py --db synapt_example_database.sqlite3 synapt_example_inputs/synapt_example.txt --tpoints 1000 --xaxis t --yaxis Temp --vin 38

python spa_db_refrag.py --db synapt_example_database.sqlite3 synapt_example_inputs/synapt_example.txt --xaxis z --yaxis frac --dH 261.69 --dS 360.55

python spa_db_listruns.py synapt_example_database.sqlite3
pause