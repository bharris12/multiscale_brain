parallel  \
		-j 4 \
		~/miniconda3/envs/python3_base/bin/python3  ~/python_utils/build_networks.py \
  		--name /data/bharris/single_cell_data/hemberg_pancreas/{1}_raw.loom \
  		--partition-col metacluster \
  		--filter-col metacluster \
  		--filter-value outliers \
  		--genes-fn /home/bharris/biccn_paper/data/pancreas/hemberg_highly_expressed_genes.csv \
  		--dense \
  		--aggregate \
  		--alternate-name {1} \
  		--output-path /home/bharris/biccn_paper/data/pancreas/networks/ ::: \
  		baron lawlor seger muraro


