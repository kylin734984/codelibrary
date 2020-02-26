if (!require(BiocManager)){
  install.packages("BiocManager")
}
if (!require(biomaRt)){
  BiocManager::install('biomaRt')
}
build = 38
querysnps = c('rs100094', 'rs10')
attrs <- c('refsnp_id','refsnp_source',"refsnp_source_description" ,'chr_name','chrom_start','chrom_end','minor_allele','minor_allele_freq','minor_allele_count','consequence_allele_string')
hosts = c( "may2009.archive.ensembl.org", "grch37.ensembl.org" , 'ensembl.org')
if (build == 36){
  host <- hosts[1]
}else if (build==37){
  host = hosts[2]
}else if (build==38){
  host = hosts[3]
}
mart <- useMart(biomart="ENSEMBL_MART_SNP", host=host, path = "/biomart/martservice")
hsnp <- useDataset(dataset = 'hsapiens_snp', mart=mart)
getBM(attrs, filters = "", values = querysnps, mart=hsnp)