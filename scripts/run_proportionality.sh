parallel  --j 4 \
	--joblog ~/biccn_paper/data/networks/proportionality/job.log \
	/opt/R/R-4.0.0/bin/Rscript ~/biccn_paper/scripts/compute_agg_proportionality.r \
	--dataset {1} \
	--level {2} :::: \
	<(cat /home/bharris/biccn_paper/data/dataset_dict_biccn_sets_7.csv  | grep andata | tr ',' '\n' | tail -n -7| sed 's/h5ad/loom/g') ::: \
	class_label subclass_label joint_cluster_label
