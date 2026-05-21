MPRA_DIR="/pasteur/helix/projects/evo_immuno_pop/MPRA"

source ${MPRA_DIR}/MPRAflow/conda/envs/mpraflow_py36/bin/activate


for LIB in LIB1 LIB2;
  do
  ROOT="assoc_B7505-4_${LIB}_noMapQ_270M_exactMatch"
  python ${MPRA_DIR}/scripts/02_assoc_barcode_lib_2/02b_export_pickleFiles.py ${ROOT}
  gzip ${MPRA_DIR}/MPRAflow/outs/${ROOT}/${ROOT}_coords_to_barcodes.tsv
  done;
