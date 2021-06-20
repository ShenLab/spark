#!/bin/bash
#$ -N Crumble
#$ -cwd
#$ -S /bin/bash
#$ -l h_vmem=6G
#$ -j yes 
#$ -o /dev/null

# Reduce BAM/CRAM files using crumble
# The first arg of this script is a CRAM/CRAM file or a BAM/CRAM file list
# Any file that not end with .cram/.bam will be taken as list
# The second arg is output directory.

# When list is provided, we will submit the job to SGE or use parallel to run multiple jobs
# at the same time. When a single BAM/CRAM is provided, we will run crumble directly
set -o pipefail
set -e 


CRAM=$1
OUTDIR=$2
ENGINE=$3

if [[ ! -f $CRAM ]]; then
    echo "Must provide BAM/CRAM file or BAM/CRAM list!" 2>&1
    exit 1
fi

if [[ ! -d $OUTDIR ]]; then
	echo "Create output directory $OUTDIR";
	mkdir -p $OUTDIR
fi

if [[ $CRAM =~ \.bam$ ]]; then
	BASE=$(basename $CRAM)
	crumble -O bam $CRAM $OUTDIR/$BASE
	samtools index $OUTDIR/$BASE
	samtools flagstat $OUTDIR/$BASE > $OUTDIR/$BASE.flagstat
elif [[ $CRAM =~ \.cram$ ]]; then
	BASE=$(basename $CRAM)
	crumble -O cram $CRAM $OUTDIR/$BASE
	samtools index $OUTDIR/$BASE
	samtools flagstat $OUTDIR/$BASE > $OUTDIR/$BASE.flagstat
else
	if [[ ! -z $SGE_TASK_ID ]]; then
		CRAMFILE=$(awk -v NLINE=$SGE_TASK_ID 'NR==NLINE{print $1}' $CRAM)
		#IID=$(awk -v NLINE=$SGE_TASK_ID 'NR==NLINE{print $2}' $CRAM)
		BASE=$(basename $CRAMFILE)
		if [[ $CRAMFILE =~ \.cram$ ]]; then
			crumble -O cram $CRAMFILE $OUTDIR/$BASE
			samtools index $OUTDIR/$BASE
			samtools flagstat $OUTDIR/$BASE > $OUTDIR/$BASE.flagstat
		elif [[ $CRAMFILE =~ \.bam$ ]]; then
			crumble -O bam $CRAMFILE $OUTDIR/$BASE
			samtools index $OUTDIR/$BASE
			samtools flagstat $OUTDIR/$BASE > $OUTDIR/$BASE.flagstat
		fi
	else
		if [[ -z $ENGINE ]]; then
			echo "Must provide ENGINE!" 2>&1
		fi
		if [[ $ENGINE == 'SGE' ]]; then
			NBAMS=$(cat $CRAM | wc -l)
			qsub -t 1:$NBAMS $0 $CRAM $OUTDIR 
		elif [[ $ENGINE == 'BASH' ]]; then
			NCPUS=$(cat /proc/cpuinfo  | grep '^processor' | wc -l)
			NJOBS=$(( NCPUS/2 ))
			awk '{print $1}' $CRAM | parallel --eta --jobs $NJOBS "$0 {} $OUTDIR"
		else 
			echo "Cannot recognize execution engine: $ENGINE" 2>&1
		fi
	fi
fi

