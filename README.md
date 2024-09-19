# multiqc
These scripts combine lab data, Illumina demultiplex stats, fastq multiqc and sample sheet files to provide insight to a sequencing run.


***This is my first ever thing to code or publish on such forum. Please, be kind.***

Only started learning coding and python 2 months ago when my Dad passed away. As a consolation my Mom and I thought that it would be good for both of us to do something productive, so she gave me a couple of lectures on python and after that I googled the rest.


***Aim and goals***

Aim of this project first was only to substitute an already existing R script for evaluating results from an Illumina sequencing (NextSeq). 
As a big picture, ideally, it would be generating a complete report from a sequencing run.

As I had learned fundamentals of R previously and I though I could manage it, so here's the result.

Most of the work was done with VSCode. I also tested the scripts under Ubuntu 22.04 LTS with python 3.12 and they work just fine.
Also found an issue (for some reason os won't work under VSC) which problem does not appear in command line, so comment/uncomment that part according to your preferences.
 
Also I'd like to track my development, therefore I publish all 3 versions of the script and currently working on a fourth which applies an object oriented solution.
In these I had three approaches to handle, clean and merge the input data: a list (Nextseq_run_QC.py), a dictionary (Nextseq_run_QC_dict.py) and a pandas DataFrame(Nextseq_run_QC_df.py).


***Description of the scripts***

All scripts have pandas dependencies, therefore you should install them first if you haven't already.
Plotting script has dash and plotly dependencies. (dash_fundamentals.py)

The script requires:
- 4 .xlsx files from the 4 different sets of samples with pre-sequencing data (data structure matches, but everything else is fictional including sample names or projects),
- a .fastq file, 
- a demultiplex file in .csv,
- a sample sheet in .csv.

The script:
- in the demultiplex stat create a "Yield" column with "PASS" and "FAILED" values for samples over a given threshold of reads (now set to 75000),
- merge the rows with same SampleID in fastq file, rename samples to match with the one in sample sheet, 
- merge reads from Read1 and Read2 and calculates average from GC content sums the others,
- reads and merge .xlsx files and 
- merge the result with the sample sheet, the demultiplex stats and the fastq,
- generates 2 .csv from this
- and a third one with basic statistics.

Note:
An additional json file is generated with the third one for further visualizations (to embed with js on an html file to make it look nicer) which is also under construction.

At the moment you can take a quick look at the report by running this dash_fundamentals with `python dash_fundamentals.py` in command line and visit http://127.0.0.1:8050/ in your web browser.
It's ugly and basic, but hey it works.
