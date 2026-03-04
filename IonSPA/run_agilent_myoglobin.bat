cd /d %~dp0
mpiexec -np 11 python spa_mpirun.py --db agilent_example_database.sqlite3 agilent_example_inputs/agilent_example.txt

python spa_db_fit.py --db agilent_example_database.sqlite3 agilent_example_inputs/agilent_example.txt --tpoints 1000 --outfile agilent_outputs.txt --plotfile screen --plotmap screen --weighting combined

python spa_db_plot.py --db agilent_example_database.sqlite3 agilent_example_inputs/agilent_example.txt --tpoints 1000 --xaxis t --yaxis Temp --vin 29

python spa_db_refrag.py --db agilent_example_database.sqlite3 agilent_example_inputs/agilent_example.txt --xaxis t --yaxis frac --dH 83.89 --dS -43.20

python spa_db_listruns.py agilent_example_database.sqlite3
pause