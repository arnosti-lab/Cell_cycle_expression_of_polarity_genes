cd ./Cell_cycle_expression_of_polarity_genes/

mkdir Raw_data

cd Raw_data

# Karaiskos et al 2017, Drosophila embryo  --------------------------------

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95025/suppl/GSE95025_high_quality_cells_digital_expression.txt.gz
gunzip GSE95025_high_quality_cells_digital_expression.txt.gz

# Deng etal. 2019, Drosophila wing disc -----------------------------------

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133204/suppl/GSE133204_RAW.tar
tar -xvf GSE133204_RAW.tar
gunzip GSM3902311_frt82b_normalization_data.csv.gz

