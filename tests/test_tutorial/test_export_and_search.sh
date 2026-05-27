#!/usr/bin/env bash

function glob_exists {
    find "$1" -name "$2" -mindepth 1 -maxdepth 1 | grep -q .
    return "$?"
}

search=false
while [ "$#" -gt 0 ]; do
    case "$1" in
	"--search" | "-s")
	    search=true
	    ;;
	*)
	    >&2 echo "Unrecognized argument $1."
	    exit 1
	    ;;
    esac
    shift
done
set -e
if ! which sudo; then
    function sudo {
	"$@"
    }
fi
if [ "$(uname)" = "Darwin" ]; then
    sudo sysctl -w kern.maxfilesperproc=30000
fi
cd rna_clique
export RNA_CLIQUE="$PWD"
cd ..
[ -d "tutorial" ]
cd tutorial
export TUTORIAL_DIR="$PWD"
cd "$RNA_CLIQUE"
eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"
conda activate rna-clique
mkdir "$TUTORIAL_DIR"/ec_genomes
cd "$TUTORIAL_DIR"/ec_genomes
case "$(uname)" in
    Linux)
	wget https://cs.uky.edu/~acta225/rna_clique/ec_genomes.tar.xz
	;;
    Darwin)
	curl -L -O https://cs.uky.edu/~acta225/rna_clique/ec_genomes.tar.xz
	;;
    *)
	>&2 echo "Unrecognized system $(uname)."
	exit 1
	;;
esac
tar xJvf ec_genomes.tar.xz
[ -f "e19_scaffolds-JAFEMN000000000.fasta" ]
[ -f "e4305_Mas339_20200623_Ref_Scaffolds_CLS.fasta" ]
cd "$RNA_CLIQUE"
python export_and_search.py -C "$TUTORIAL_DIR"/rna_clique_out/config.yaml \
                            -Q "$TUTORIAL_DIR"/ec_genomes/*.fasta \
                            -X "$TUTORIAL_DIR"/full_ec_search_out \
                            -e 1e-99
[ -d "$TUTORIAL_DIR/full_ec_search_out/rna_clique_out" ]
glob_exists "$TUTORIAL_DIR/full_ec_search_out/rna_clique_out/export/" "*.fasta"
grep -q "$TUTORIAL_DIR/full_ec_search_out/rna_clique_out/export/"*.fasta \
     -e '^>-.*:SRR.*$'
grep -q \
     "$TUTORIAL_DIR/full_ec_search_out/rna_clique_out/export/all_ideal.fasta" \
     -e '^>.*:SRR.*:ideal_component_.*$'
# "Search results"
for d in "$TUTORIAL_DIR/full_ec_search_out/rna_clique_out/search_e19_scaffolds-JAFEMN000000000" \
             "$TUTORIAL_DIR/full_ec_search_out/rna_clique_out/search_e4305_Mas339_20200623_Ref_Scaffolds_CLS"; do
    [ -d "$d" ]
    [ -f "$d/queries.sam" ]
    [ -f "$d/stats" ]
    [ -f "$d/subjects.fasta" ]
    seqs="$(grep "^>" -c "$d/subjects.fasta")"
    grep -q "\"seqs\": $seqs" "$d/stats"
done
set -o pipefail
components="$(grep "$TUTORIAL_DIR"/full_ec_search_out/*/search_*/subjects.fasta \
		   -e '>' | grep '[^:]*:[^:]*$' -o | sort -u | sort -t: -k2 | \
		   grep 'SRR2321385' | cut -d: -f2)"
set +o pipefail
(
    while read -u 5 -r line; do
	cat "$TUTORIAL_DIR/full_ec_search_out/rna_clique_out/export/$line.fasta"
    done <<< "$components"
) > "$TUTORIAL_DIR"/matching_components.fastas
# TODO: Add IGV-related tests
if [ "$search" = true ]; then
    blastn -query "$TUTORIAL_DIR"/matching_components.fasta -evalue 1e-99 \
	   -remote -db nr -outfmt 6 -out "$TUTORIAL_DIR"/remote_results
    [ -F "$TUTORIAL_DIR"/remote_results ]
fi
# Subset analysis
python export_and_search.py -C "$TUTORIAL_DIR"/infected_subset_out/config.yaml \
                            -Q "$TUTORIAL_DIR"/ec_genomes/*.fasta \
                            -X "$TUTORIAL_DIR"/subset_ec_search_out \
                            -e 1e-99
# TODO: Use jq to check number of components.
