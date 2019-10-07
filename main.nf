#!/usr/bin/env nextflow

// Folder with individual  cohort Tumor_Sample_Barcode files (1 file == 1 case-control pair)
Channel
    .fromPath("${params.maf_folder}",  type: 'dir' )
    .ifEmpty { exit 1, "Path to input --maf_folder is incorrect or the folder contains no .maf files." }
    .set { maf_folder_channel}


process merge_maf_files {
  tag "$maf_folder"
  publishDir params.outdir, mode: 'copy'
  container 'lifebitai/mtsgmftls:101'
  echo true

  input:
  file(maf_folder) from maf_folder_channel

  output:
  file("cohort_maf.RData") into merged_maf_file_channel

  script:
  """
  #!/usr/bin/env Rscript

  library(maftools)
  
  cohort_maf   <- maftools::merge_mafs(list.files(paste0("${maf_folder}", "/")))

  # Saving on object in RData format
  #save(cohort_maf, file = "cohort_maf.RData")
  """
}

process run_mutsig_analysis {
  tag "$maf_folder"
  publishDir params.outdir, mode: 'copy'
  container 'lifebitai/mtsgmftls:101'

  input:
  file(maf_folder) from maf_folder_channel

  output:
  file("{MultiQC,multiqc_report.html}") into results

  script:
  """
  # copy the docker bin into pwd
  mkdir bin
  cp -r /opt/conda/envs/mutationsignatures-nf/bin/* bin/

  # copy the docker bin into pwd
  mkdir mafs/
  cp -r "${maf_folder}" mafs/
  mv  mafs/ bin/mafs/ 
  
  cd bin

  # copy the rmarkdown into the workdir
  R -e "rmarkdown::render('maftools_report.Rmd', params = list(maf_folder='../${maf_folder}'), output_file='maftools_report.html')"

  cd ..

  mkdir MultiQC && mv bin/maftools_report.html MultiQC/multiqc_report.html
  """
}