#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH --qos=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=18
#SBATCH -A w21_magfusion
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=giannis_kx@lanl.gov
#SBATCH --export=ALL

date
module load python/3.6-anaconda-5.0.1

srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/db_tear/dt2-ftle/t-1055/10-6/" "raw_t-1055-e-6.sav" 1e-6 &
srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/db_tear/dt2-ftle/t-1055/10-8/" "raw_t-1055-e-8.sav" 1e-8 &
srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/db_tear/dt2-ftle/t-1055/10-9/" "raw_t-1055-e-9.sav" 1e-9 &
srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/db_tear/dt2-ftle/t-1055/10-10/" "raw_t-1055-e-10.sav" 1e-10 &
#srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/db_tear/dt2-ftle/t-1800/" "raw_t-1800.sav" &
#srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/db_tear/dt2-ftle/t-1900/" "raw_t-1900.sav" &
#srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/11_visc_old_nodiff/ftle/t-500/" "raw_t-500.sav" &
#srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/11_visc_old_nodiff/ftle/t-550/" "raw_t-550.sav" &
#srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/11_visc_old_nodiff/ftle/t-600/" "raw_t-600.sav" &
#srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/11_visc_old_nodiff/ftle/t-620/" "raw_t-620.sav" &
#srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/11_visc_old_nodiff/ftle/t-640/" "raw_t-640.sav" &
#srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/11_visc_old_nodiff/ftle/t-654/" "raw_t-654.sav" &
srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/ftle-chipar/t-920/10-8/" "raw_t-920-e-8.sav" 1e-8 &
srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/ftle-chipar/t-920/10-9/" "raw_t-920-e-9.sav" 1e-9 &
srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/ftle-chipar/t-920/10-10/" "raw_t-920-e-10.sav" 1e-10 &
#srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/ftle-chipar/t-1070/" "raw_t-1070.sav" &
#srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/ftle-chipar/t-1120/" "raw_t-1120.sav" &
#srun --exclusive --nodes 1 --ntasks 1 python readborg3.py "/lustre/scratch4/turquoise/giannis_kx/pixie3d/iter/int_kink/11/ftle-chipar/t-1160/" "raw_t-1160.sav" &
wait