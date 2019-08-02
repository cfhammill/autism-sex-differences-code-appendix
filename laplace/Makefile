model_selection_assets.rda: generate_model_selection_assets.Rmd
	Rscript -e "rmarkdown::render('generate_model_selection_assets.Rmd')"

model_selection_results.rda: model_selection_assets.rda run_model_selection.Rmd
	Rscript -e "rmarkdown::render('run_model_selection.Rmd')"

spatial_analyses_finished.txt: model_selection_results.rda model_selection_assets.rda run_spatial_analyses.Rmd
	Rscript -e "rmarkdown::render('run_spatial_analyses.Rmd')"

merged_experiments.rds: spatial_analyses_finished.txt
	Rscript -e "rmarkdown::render('process_spatial_analyses.Rmd')"

figures_done.txt: merged_experiments.rds
	Rscript -e "rmarkdown::render('output_results.Rmd')"

demographics_summaries.html: demographic_summaries.Rmd
	Rscript -e "rmarkdown::render('demographic_summaries.Rmd')"

vertexwise_analysis.html: vertexwise_analysis.Rmd
	Rscript -e "rmarkdown::render('$<')"

vertexwise_results.html: vertexwise_results.Rmd
	Rscript -e "rmarkdown::render('$<')"
