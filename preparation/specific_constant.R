# check species
check_species <- function(params) {
  species_genes <- c()
  if (tolower(params$species) %in% c("human", "homo sapiens")) {
    species_genes$s.genes <- c(
      "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1",
      "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS",
      "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
      "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN",
      "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"
    )
    species_genes$g2m.genes <- c(
      "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2",
      "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3",
      "PIMREG", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11",
      "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "JPT1", "CDC20", "TTK",
      "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2",
      "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
      "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"
    )

    species_genes$hb.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")

    species_genes$gdt.marker <- c("TRGV9", "TRDV2", "TRAC", "TRGC1", "TRAV1-2")

    return(species_genes)
  } else if (tolower(params$species) %in% c("mouse", "mus musculus")) {
    species_genes$s.genes <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1",
                           "Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells",
                           "Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7",
                           "Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin",
                           "Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8"
    )
    species_genes$g2m.genes <- c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2",
                             "Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3",
                             "Pimreg","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11",
                             "Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Jpt1","Cdc20","Ttk",
                             "Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2",
                             "Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf",
                             "Nek2","G2e3","Gas2l3","Cbx5","Cenpa"
    )
    species_genes$hb.genes <- c("Hba-a1", "Hbb-bs", "Hba-a2", "Hbb-bt", "	Hbb-bh1")
    
    return(species_genes)
    
  }
  else {
    stop(sprintf("ERROR: Unknown or missing species: %s.\n These species are support: %s\n", params$species, base_params$support_species))
  }
  # TODO mouse
}
