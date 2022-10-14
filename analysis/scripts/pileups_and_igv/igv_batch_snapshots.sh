set -e
usage(){
    echo -e "usage: $0 sample_dirname" 
    exit 0
}
sample_dirname=$1
if [[ -z "$sample_dirname" ]]; then usage; fi

sample_dirname=$(realpath $sample_dirname)
sample_name=$(basename $sample_dirname)
script_name=$(dirname $sample_dirname)"/igv_batch_scripts/${sample_name}_batch.txt"
snapshot_dirname=$(dirname $sample_dirname)"/igv_snapshots/"

mkdir -p $(dirname ${script_name})
mkdir -p ${snapshot_dirname}

echo -e "snapshotDirectory ${snapshot_dirname}" > ${script_name}
for dirname in $(find ${sample_dirname} -type d -maxdepth 1 | tail -n+2); do
    echo -e "new\ngenome ${dirname}/${sample_name}_induced_ref.fa.gz\n" >> ${script_name}
    echo -e "load ${dirname}/${sample_name}_mapped_to_induced.bam\n viewaspairs\n" >> ${script_name}
    tool_name=$(basename $dirname)
    IFS=$'\n'; for gene_line in $(grep -E "DBLMSP|AMA1|P48_45" "${dirname}/${sample_name}_gene_regions.bed")
    do
        IFS=$'\t'; elems=($gene_line)    
        adj_start=$((${elems[1]} + 1))
        reg="${elems[0]}:${adj_start}-${elems[2]}"
        gene_name=${elems[3]}
        echo -e "goto ${reg}\nsnapshot ${tool_name}_${sample_name}_${gene_name}.png">> ${script_name}
    done
    echo "exit" >> ${script_name}
done

# Below runs igv without GUI (modified from
# https://gist.github.com/stevekm/ac76c0c2fa4ee89db8ce2421cc6fbffc)
prefix=`dirname $(readlink $(which igv.sh) || echo igv.sh)`
(Xvfb :10 &) && DISPLAY=:10 \
    java -showversion --module-path="${prefix}/lib" -Xmx4g \
        @"${prefix}/igv.args" \
        -Dapple.laf.useScreenMenuBar=true \
        -Djava.net.preferIPv4Stack=true \
        --module=org.igv/org.broad.igv.ui.Main -b ${script_name} \
&& killall Xvfb


