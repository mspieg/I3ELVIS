#! /bin/bash

# a single job can depend on all jobs by the same user with the same name
jid1=$sbatch slurm-submit.sh
jid2=$(sbatch --dependency=afterok:$jid1 slurm-submit.sh)
jid3=$(sbatch --dependency=afterok:$jid2 slurm-submit.sh)
jid4=$(sbatch --dependency=afterok:$jid3 slurm-submit.sh)
jid5=$(sbatch --dependency=afterok:$jid4 slurm-submit.sh)
jid6=$(sbatch --dependency=afterok:$jid5 slurm-submit.sh)
jid7=$(sbatch --dependency=afterok:$jid6 slurm-submit.sh)
jid8=$(sbatch --dependency=afterok:$jid7 slurm-submit.sh)
jid9=$(sbatch --dependency=afterok:$jid8 slurm-submit.sh)
jid10=$(sbatch --dependency=afterok:$jid9 slurm-submit.sh)
jid11=$(sbatch --dependency=afterok:$jid10 slurm-submit.sh)
jid12=$(sbatch --dependency=afterok:$jid11 slurm-submit.sh)
jid13=$(sbatch --dependency=afterok:$jid12 slurm-submit.sh)
jid14=$(sbatch --dependency=afterok:$jid13 slurm-submit.sh)
jid15=$(sbatch --dependency=afterok:$jid14 slurm-submit.sh)
jid16=$(sbatch --dependency=afterok:$jid15 slurm-submit.sh)
jid17=$(sbatch --dependency=afterok:$jid16 slurm-submit.sh)