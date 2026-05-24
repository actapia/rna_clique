#!/usr/bin/env bash
assemble=false
if ! which sudo; then
    function sudo {
	"$@"
    }
fi
parallel=true
assemble=false
while [ "$#" -gt 0 ]; do
    case "$1" in
	"--no-parallel" | "-P")
	    parallel=false
	    ;;
	"--assemble" | "-a")
	    assemble=true	    
	    ;;
	*)
	    >&2 echo "Unrecognized argument $1."
	    exit 1
	    ;;
    esac
    shift
done
set -e
cd rna_clique
export RNA_CLIQUE="$PWD"
cd ..
[ -d "tutorial" ]
cd tutorial
export TUTORIAL_DIR="$PWD"
eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"
conda activate rna-clique
cd ..
if [ "$assemble" = true ]; then
    case "$(uname)" in
	Linux)
	    sudo apt update
	    sudo NEEDRESTART_MODE=a apt -y install trinityrnaseq
	    if [ "$parallel" = true ]; then
		sudo NEEDRESTART_MODE=a apt -y install parallel
	    fi
	    if ! compgen -G "sratoolkit*.tar.gz"; then
		wget --no-verbose https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
	    fi
	    ;;
	Darwin)
	    brew tap brewsci/bio
	    brew install trinity parallel
	    if [ "$parallel" = true ]; then
		brew install parallel
	    fi
	    if ! compgen -G "sratoolkit*.tar.gz"; then
		case "$(uname -m)" in
		    x86_64)
			curl -L -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz
			;;
		    arm64)
			curl -L -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac-arm64.tar.gz
			;;
		    *)
			>&2 echo "Unrecognized machine type $(uname -m)."
			exit 1
			;;
		esac
		
	    fi
	    ;;
	*)
	    >&2 echo "Unrecognized system $(uname)."
	    exit 1
	    ;;
    esac
    compgen -G "sratoolkit.current-*.tar.gz"
    tar xzvf sratoolkit.current-*.tar.gz
    conda install -y lxml requests
    if ! [ -d "download_sra" ]; then
	git clone https://github.com/actapia/download_sra
    fi
    export PATH="$PATH:$PWD/download_sra"
    cd "$TUTORIAL_DIR"
    mkdir trinity_assemblies
    cd trinity_assemblies
    if [ "$parallel" = true ]; then
	tail -n+2 "$RNA_CLIQUE/docs/tutorials/reads2tree/tall_fescue_accs.csv" | \
	    cut -d, -f1 | \
	    parallel --jobs "$PARALLEL_JOBS" \
		     "download_sra.sh -j 0 -r {}; Trinity --seqType fq --max_memory $TRINITY_MEMORY --single {}.fastq  --output trinity_{} --CPU $TRINITY_THREADS; rm {}.fastq;"
    else
	while read -u 6 -r line; do
	    download_sra.sh -j 0 -r "$line"
	    Trinity --seqType fq --max-memory "$TRINITY_MEMORY" \
		    --single "$line.fastq" --output "trinity_$line" \
		    --CPU "$TRINITY_THREADS"
	    rm "$line.fastq"
	    #rm -r "trinity_$line"
	done 6< <(tail -n+2 "$RNA_CLIQUE/docs/tutorials/reads2tree/tall_fescue_accs.csv" | cut -d, -f1)
    fi
else
    cd "$TUTORIAL_DIR"    
    mkdir trinity_assemblies
    cd trinity_assemblies
    case "$(uname)" in
	Linux)
	    wget "http://rna-clique-data.s3-website.us-east-2.amazonaws.com/trinity_assemblies.zip"
	    ;;
	macOS)
	    curl -L -O "http://rna-clique-data.s3-website.us-east-2.amazonaws.com/trinity_assemblies.zip"
	    ;;
	*)
	    >&2 echo "Unrecognized system $(uname)."
	    exit 1
	    ;;
    esac
    unzip trinity_assemblies.zip
    for f in trinity_*.Trinity.fasta; do 
	( grep "$f" -e '^>' | grep -v -e 'TRINITY_.*_c.*_g.*_i.*' ) && exit 1
    done
fi
cd "$TUTORIAL_DIR/trinity_assemblies"
for f in trinity_*/; do
    [ -f "$f/salmon_outdir/quant.sf" ]
done
mkdir with_tpm
for f in SRR2321385 SRR2321388 SRR7990321 SRR8003736 \
		    SRR8003761 SRR8003762; do
    python "$RNA_CLIQUE/docs/tutorials/nonspades/add_tpm.py" \
	   "trinity_$f.Trinity.fasta" \
	   "trinity_$f/salmon_outdir/quant.sf" > "with_tpm/$f.fasta"
done
for f in with_tpm/*.fasta; do
    ( grep "%f" -e '^>' | grep -v -e 'tpm.*_TRINITY_.*_c.*_g.*_i.*' ) \
	&& exit 1
done
mkdir integer_ids
for f in SRR2321385 SRR2321388 SRR7990321 SRR8003736 \
		    SRR8003761 SRR8003762; do
    python "$RNA_CLIQUE/docs/tutorials/nonspades/assign_gene_ids.py" \
	   < "with_tpm/$f.fasta" > "integer_ids/$f.fasta"
done
for f in integer_ids/*.fasta; do
    ( grep "%f" -e '^>' | \
	  grep -v -e 'tpm.*_TRINITY_.*_c.*_g.*_gid.*_i.*' ) && exit 1
done
rm -r with_tpm
for f in SRR2321385 SRR2321388 SRR7990321 SRR8003736 \
		    SRR8003761 SRR8003762; do
    mkdir "$f"
    mv "integer_ids/$f.fasta" "$f"/transcripts.fasta
done
rmdir integer_ids
cd "$RNA_CLIQUE"
python rna_clique.py "$TUTORIAL_DIR"/trinity_assemblies/SRR* -n 50000 \
       -O "$TUTORIAL_DIR/trinity_rna_clique_out" \
       -p '^.*tpm([0-9]+(?:\.[0-9]+)).*gid([0-9]+)_i([0-9]+)'

python export_matrix.py --format table \
       --header \
       -O "$TUTORIAL_DIR/trinity_rna_clique_out"
PYTHONPATH="." python docs/tutorials/reads2tree/make_pcoa.py \
                      "$TUTORIAL_DIR"/trinity_rna_clique_out
[ -f "$TUTORIAL_DIR"/trinity_rna_clique_out/pcoa_2d.svg ]
[ -f "$TUTORIAL_DIR"/trinity_rna_clique_out/pcoa_3d.svg ]

