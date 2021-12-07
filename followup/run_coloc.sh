pmid_seq="29892013 24390342 29083406 30929738 26502338 26192919 26192919 26192919 28067908 28067908 28067908"
label_seq="AE RA_GWASmeta_European Allergy ASTHMA sle IBD CD UC ibd cd uc"

set ${label_seq}
for pmid in ${pmid_seq}
do
echo "Rscript --no-restore --no-save script/coloc_gwas_immu.R ${pmid} $1"
echo "Rscript --no-restore --no-save script/coloc_immu.R ${pmid} $1"
shift
done

set ${label_seq}
for pmid in ${pmid_seq}
do
echo "Rscript --no-restore --no-save script/coloc_visual_immu.R ${pmid} $1"
shift
done

