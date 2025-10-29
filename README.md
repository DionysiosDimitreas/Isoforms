# Isoforms

## Building an isoform-aware CellChat database

1. Prepare 3did domain interaction rulebook

```
# Filter out PFAM domain-domain interaction information (copied from CanIsoNet)
less databases/3did/v2022_01/3did_flat.gz | grep "^#=ID" | cut -f4,5 | perl -ane '$F[0]=~ s/.*(PF\d+).*/$1/; $F[1]=~ s/.*(PF\d+).*/$1/; print "$F[0]\t$F[1]\n$F[1]\t$F[0]\n"' | sort -u | gzip > 3did_pfamInteractions_1804.tsv.gz

# Reformat to DIMA format
python3 scripts/3did2dima.py 3did_pfamInteractions_1804.tsv.gz | gzip > 3did_1804.tbl.gz
```

2. Prepare Ensembl fasta files for PFAM domain identification

```
# Create directory structure for Ensembl protein FASTA files
mkdir -p fasta/ensembl_v75/human

# Create ENSEMBL protein FASTA files
zless databases/ensembl/v75/Homo_sapiens.GRCh37.75.pep.all.fa.gz | perl -ne 'if(/^>(ENSP\d+) .*/){$i=$1; $h{$i}=$_; next} else{$h{$i}.=$_} END{foreach $i(sort keys %h){open(O, ">fasta/ensembl_v75/human/$i.fasta"); print O $h{$i}; close(O) }}'
```

3. Run PfamScan for domain identification
This repository includes a script to automate running `hmmscan` on a batch of FASTA files and merging the results.

  **a) Prerequisites:**

You must have [HMMER](http://hmmer.org/) (which includes `hmmscan` and `hmmpress`) installed and available in your system's PATH.

  **b) Run the Script**
First, make the script executable (you only need to do this once):
```bash
chmod +x scripts/run_pfam_scan.sh
```

Now, run the script by providing three arguments:
1.  Your input FASTA directory
2.  The path to your `Pfam-A.hmm` database
3.  A directory name for the results

**Example:**
```bash
./scripts/run_pfam_scan.sh fasta/ensembl_v75/human/ ~/databases/Pfam/Pfam-A.hmm ./pfam_scan_results
```

The script will create the output directory, run `hmmscan` on every `.fasta` file, and merge all results into a single compressed file named `pfam_all_hits.domtbl.gz`.

