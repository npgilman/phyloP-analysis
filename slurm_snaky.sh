#!/bin/bash
#SBATCH --job-name=Main_accelerated_regions
#SBATCH --account=kgraim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ngilman@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --time=744:00:00
#SBATCH --output=logs/%j.stdout
#SBATCH --error=logs/%j.stderr

pwd; hostname; date
module load snakemake/7.32.4

snakemake --cluster "sbatch -A {cluster.account} -q {cluster.qos}  -c {cluster.cpus-per-task} -N {cluster.Nodes} -t {cluster.runtime} --mem {cluster.mem} -J {cluster.jobname} \
 --mail-user={cluster.mail} --output {cluster.out} --error {cluster.err}" --cluster-config cluster_config.json --jobs 50 --latency-wait 20 --rerun-incomplete --use-envmodules --keep-going 
