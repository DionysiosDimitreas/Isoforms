# Isoforms

## Building an isoform-aware CellChat database

1. Prepare 3did domain interaction file

```
# Filter out PFAM domain-domain interaction information (copied from CanIsoNet)
less databases/3did/v2022_01/3did_flat.gz | grep "^#=ID" | cut -f4,5 | perl -ane '$F[0]=~ s/.*(PF\d+).*/$1/; $F[1]=~ s/.*(PF\d+).*/$1/; print "$F[0]\t$F[1]\n$F[1]\t$F[0]\n"' | sort -u | gzip > 3did_pfamInteractions_1804.tsv.gz

# Reformat to DIMA format
python3 scripts/3did2dima.py 3did_pfamInteractions_1804.tsv.gz | gzip > 3did_1804.tbl.gz
```
