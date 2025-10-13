# download the hyphy def file
```bash
# can mkdir XX; cd XX first
wget https://raw.githubusercontent.com/maggieMCKO/misc/main/absrel/setup/hyphy_man3.def
```

# build an apptainer image
## in scc legacy cluster - software stack `scc-lmod`
```bash
module load rev/23.12 apptainer/1.1.9
# or search first by 
# module spider apptainer

apptainer build hyphy_v2.5.69.sif hyphy_man3.def
```

## in scc new cluster (project id) - software stack `scc-lmod`
```bash
module load gcc/13.2.0 apptainer/1.2.5
# or search first by 
# module spider apptainer

apptainer build hyphy_v2.5.69.sif hyphy_man3.def
```

# testing
```bash
export image="$HOME/hyphy_v2.5.69.sif"
apptainer exec ${image}  hyphy absrel --help

```

# template script
```bash
# Load the required module for Apptainer
module load apptainer

# Define input variables
ALIGN=$SAMPLE
TREE=$TREE
OUTPUT=$tmp_dir/$stem.ABSREL_srv.json
OUTPUTLOG=$tmp_dir/$stem.ABSREL_srv.log

# Set the path to the Apptainer image
export image="$HOME/hyphy_v2.5.69.sif"

# Execute the ABSREL analysis
apptainer exec ${image} sh -c "cd $tmp_dir && hyphy absrel --alignment $ALIGN --tree $TREE --output $OUTPUT --srv Yes"
```
