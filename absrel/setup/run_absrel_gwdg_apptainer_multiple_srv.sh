#!/bin/bash
#SBATCH -n 4                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-03:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p scc-cpu           # Partition to submit to
#SBATCH --mem=12G         # Memory pool for all cores (see also --mem-per-cpu)

# wd: /user/mko/u14219/findgene
module purge
# module load rev/23.12 apptainer/1.1.9
module load gcc/13.2.0 apptainer/1.2.5

export singularity_image="$HOME/findgene/hyphy_v2.5.69.sif"

export Project_DIR="$HOME/findgene/nexus_nointernalnodes/"
# align and tree need to be matching, and tree can't have internal node lables

##HYPHY SCRIPT##

run(){
	SAMPLE=$1

	stem=$(basename ${SAMPLE})
	stem=${stem%%.*}
	tmp_dir=$(dirname ${SAMPLE})
	tmp_dir=$tmp_dir/ABSREL/$stem
	mkdir -p $tmp_dir
	cd $tmp_dir

	#set up run
	ALIGN=$SAMPLE
# 	TREE=$HYPHYDIR/$TREE
	OUTPUT=$tmp_dir/$stem.ABSREL_srv.json
	OUTPUTLOG=$tmp_dir/$stem.ABSREL_srv.log

# 	cd $HYPHYDIR # kind of necessary unless found what path I need to define
# 	./hyphy absrel --alignment $ALIGN --output $OUTPUT &> $OUTPUTLOG
# 	./hyphy absrel --alignment $ALIGN --tree $TREE --output $OUTPUT
	apptainer exec ${singularity_image} sh -c "cd $tmp_dir && hyphy absrel --alignment $ALIGN  --output $OUTPUT --srv Yes --code Universal ENV=TOLERATE_NUMERICAL_ERRORS=1" # runs
}
export -f run

find ${Project_DIR} | grep ".nex$" | parallel --will-cite run {}


# /home/mpg08/mko/Tools/hyphy/HYPHYMP BASEPATH=/home/mpg08/mko/Tools/hyphy/res/TemplateBatchFiles/ CPU=1 BranchSiteREL.bf --alignment $ALIGN --output $OUTPUT
