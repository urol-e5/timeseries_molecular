# Location of Analysis
UNITY ```/project/pi_hputnam_uri_edu/hputnam/E5_Orthologs/```


# Obtain the 3 genomes

### Apulchra
[Data Origin](https://osf.io/y8963/)
```
wget https://raw.githubusercontent.com/urol-e5/deep-dive-expression/main/D-Apul/data/Apulchra-genome.pep.faa
```
### Pevermanni
[Data Origin](https://www.genoscope.cns.fr/corals/genomes.html) 
```
wget https://gannet.fish.washington.edu/seashell/snaps/Porites_evermanni_v1.annot.pep.fa 
```
### Ptuahiniensis (our samples) = Pmeandrina (current reference)
[Data Origin](http://cyanophora.rutgers.edu/Pocillopora_meandrina/)
```
wget https://gannet.fish.washington.edu/seashell/snaps/Pocillopora_meandrina_HIv1.genes.pep.faa
```

# Generating Orthologs with Orthofinder
[Orthofinder](https://github.com/davidemms/OrthoFinder)
```
conda create -n orthofinder_env orthofinder -c bioconda -c conda-forge
```
```
nano run_orthofinder.sh
```

```
#!/bin/bash
#SBATCH --job-name=checksum_raw
#SBATCH --nodes=1 --cpus-per-task=8
#SBATCH --mem=250G  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH --time=06:00:00  # Job time limit
#SBATCH -o slurm-run_orthofinder.out  # %j = job ID
#SBATCH -e slurm-run_orthofinder.err  # %j = job ID
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=hputnam@uri.edu #your email to send notifications
#SBATCH -D /project/pi_hputnam_uri_edu/hputnam/E5_Orthologs/

# load modules needed
conda activate orthofinder_env

# run orthofinder
orthofinder -f /project/pi_hputnam_uri_edu/hputnam/E5_Orthologs/Protein_fastas/ -t 8

```
```
sbatch run_orthofinder.sh
```
