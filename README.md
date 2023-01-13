# Automated Report Example

This repository contains code to automatically download, process and analyze a publich data set for the epigenetic signature of Parkinson's disease (GSE111629, Illumina 450k array).

This includes an elastic net regression model to predict Parkinson's disease.

All results are reported in an HTML file automatically generated by RMarkdown.

The regression model doesn't actually work particularly well, but ultimately, I picked one of the first GEO data sets I could find and ran it through a fairly basic pipeline - this is more of a tech demo for the automated reports than a sophisticated machine learning project.

## Notes:
Code was written under R 4.2.1.

Due to the sheer scope of the 450k array, running the code requires over 20GB RAM. The required files will take up roughly 16.5GB of disk space.

## Order of operations:
After cloning the repository, run the files in the following order:

1. are_setup.R (installs dependencies)
2. are_download_and_process_dataR (for required methylation and clinical data)
3. are_analysis.R (this includes all analysis and will source are_analysis.Rmd to generate a report file)

This will produce results/are_analysis_report.html.

You may want to also run are_cleanup.R.
This simply removes all files that aren't needed by are_analysis.R - of course, you can also just manually delete the whole repository when done.

###Alternatively, an example report is already included: results/
