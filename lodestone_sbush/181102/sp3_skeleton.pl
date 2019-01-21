=head
BEFORE USAGE:
ssh sbush@eddie3.ecdf.ed.ac.uk
=cut

use strict;
use warnings;
use POSIX qw(ceil strftime);

# REQUIREMENTS
my $home	= '/home/sbush';
my $progs   = "$home/programs";
my $in_file = "$home/sp3/list_of_run_ids.txt"; # contains one SRA run accession per line

# OUTPUT
my $out_dir = '/exports/eddie/scratch/sbush/sp3_output';
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }

# SOFTWARE USED
my $fatal = 0;
## paths to aligners
my $bwa_path    = "$progs/bwa/bwa";
my $stampy_path = "$progs/stampy-1.0.32/stampy.py";
if (!(-e($bwa_path))) 	 { $fatal++; print "ERROR: cannot find $bwa_path\n";    }
if (!(-e($stampy_path))) { $fatal++; print "ERROR: cannot find $stampy_path\n"; }
## paths to assemblers
my $quast_path  = '/exports/eddie/scratch/sbush/quast-5.0.1/quast.py'; # "$progs/quast-5.0.1/quast.py";
my $circos_dir  = "$progs/circos-0.69-6/bin";
my $circos_path = "$progs/circos-0.69-6/bin/circos";
my $spades_path = "$progs/SPAdes-3.13.0-Linux/bin/spades.py";
if (!(-e($quast_path)))  { $fatal++; print "ERROR: cannot find $quast_path\n";  }
if (!(-d($circos_dir)))  { $fatal++; print "ERROR: cannot find $circos_dir\n";  }
if (!(-e($circos_path))) { $fatal++; print "ERROR: cannot find $circos_path\n"; }
if (!(-e($spades_path))) { $fatal++; print "ERROR: cannot find $spades_path\n"; }
## paths to callers
my $bcftools_path = "$progs/bcftools-1.7/bcftools";
if (!(-e($bcftools_path))) { $fatal++; print "ERROR: cannot find $bcftools_path\n"; }
## paths to classifiers
my $kaiju_db_p   	  = "/exports/eddie/scratch/sbush/db.p/kaiju_db.fmi"; # see http://kaiju.binf.ku.dk/server; specifically, http://kaiju.binf.ku.dk/database/kaiju_index_pg.tgz (the proGenomes database)
my $nodes_dmp_p  	  = "/exports/eddie/scratch/sbush/db.p/nodes.dmp";
my $names_dmp_p 	  = "/exports/eddie/scratch/sbush/db.p/names.dmp";
my $kaiju_db_e   	  = "/exports/eddie/scratch/sbush/db.e/kaiju_db_nr_euk.fmi"; # see http://kaiju.binf.ku.dk/server; specifically, http://kaiju.binf.ku.dk/database/kaiju_index_nr_euk_2018-02-23.tgz (the expanded nr database, containing bacteria, archaea, viruses, fungi and microbial eukaryotes)
my $nodes_dmp_e  	  = "/exports/eddie/scratch/sbush/db.e/nodes.dmp";
my $names_dmp_e 	  = "/exports/eddie/scratch/sbush/db.e/names.dmp";
my $kaiju_path 		  = "$progs/kaiju/bin/kaiju";
my $kaijuReport_path  = "$progs/kaiju/bin/kaijuReport";
my $kaijutokrona_path = "$progs/kaiju/bin/kaiju2krona";
my $mergeOutputs_path = "$progs/kaiju/bin/mergeOutputs";
my $centrifuge_path	  = "$progs/centrifuge/centrifuge";
my $centrifuge_db     = "/exports/eddie/scratch/sbush/phv"; # bacteria, archaea, viruses and human; ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p+h+v.tar.gz
my $centrifuge_index_prefix = "$centrifuge_db/p+h+v"; # a Centrifuge index comprises 3 files of the format prefix.X.cf
if (!(-e($kaiju_db_p)))   	   { $fatal++; print "ERROR: cannot find $kaiju_db_p\n";   	    }
if (!(-e($nodes_dmp_p)))  	   { $fatal++; print "ERROR: cannot find $nodes_dmp_p\n";  	    }
if (!(-e($names_dmp_p))) 	   { $fatal++; print "ERROR: cannot find $names_dmp_p\n";  	    }
if (!(-e($kaiju_path))) 	   { $fatal++; print "ERROR: cannot find $kaiju_path\n"; 	    }
if (!(-e($kaijuReport_path)))  { $fatal++; print "ERROR: cannot find $kaijuReport_path\n";  }
if (!(-e($kaijutokrona_path))) { $fatal++; print "ERROR: cannot find $kaijutokrona_path\n"; }
if (!(-e($mergeOutputs_path))) { $fatal++; print "ERROR: cannot find $mergeOutputs_path\n"; }
if (!(-e($centrifuge_path)))   { $fatal++; print "ERROR: cannot find $centrifuge_path\n";   }
if (!(-d($centrifuge_db)))     { $fatal++; print "ERROR: cannot find $centrifuge_db\n";     }
## paths to ancillary software
my $mash_path				  = "$progs/mash-Linux64-v2.1/mash";
my $ascp_path				  = "/home/sbush/.aspera/connect/bin/ascp";
my $aspera_openssh_path		  = "/home/sbush/.aspera/connect/etc/asperaweb_id_dsa.openssh";
my $bgzip_path   			  = "$progs/bin/bgzip";
my $tabix_path   			  = "$progs/bin/tabix"; # tabix and bgzip are components of the SAMtools package, commonly used to compress and index VCFs; see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/
my $seqtk_path 				  = "$progs/seqtk/seqtk";
my $fastqc_path				  = "$progs/FastQC/fastqc";
my $picard_path   			  = "$progs/picard.jar";
my $efetch_path 			  = "$progs/edirect/efetch"; # for a fuller description of the Entrez Direct tools (used to automate downloads from NCBI), see https://www.ncbi.nlm.nih.gov/books/NBK179288/
my $xtract_path 			  = "$progs/edirect/xtract";
my $esearch_path 			  = "$progs/edirect/esearch";
my $esummary_path 			  = "$progs/edirect/esummary";
my $tetyper_path			  = "/home/sbush/.local/bin/TETyper.py";
my $vcfsort_path 			  = "$progs/bin/vcf-sort"; # a component of VCFtools; see https://vcftools.github.io/index.html
my $samtools_path 			  = "$progs/samtools-1.7/samtools";
my $sambamba_path			  = "$progs/sambamba";
my $cutadapt_path			  = "/home/sbush/.local/bin/cutadapt";
my $trim_galore_path		  = "$progs/TrimGalore-0.5.0/trim_galore";
my $trimmomatic_path 		  = "$progs/Trimmomatic-0.38/trimmomatic-0.38.jar";
my $ktImportText_path		  = "$progs/bin/ktImportText";
my $vcfoverlay_path 		  = "$progs/vcflib/bin/vcfoverlay";
my $vcfallelicprimitives_path = "$progs/vcflib/bin/vcfallelicprimitives";
if (!(-e($mash_path)))    			   { $fatal++; print "ERROR: cannot find $mash_path\n";    				}
if (!(-e($ascp_path)))    			   { $fatal++; print "ERROR: cannot find $ascp_path\n";    				}
if (!(-e($bgzip_path)))   			   { $fatal++; print "ERROR: cannot find $bgzip_path\n";   				}
if (!(-e($tabix_path)))   			   { $fatal++; print "ERROR: cannot find $tabix_path\n";   				}
if (!(-e($seqtk_path)))   			   { $fatal++; print "ERROR: cannot find $seqtk_path\n";   				}
if (!(-e($fastqc_path)))   			   { $fatal++; print "ERROR: cannot find $fastqc_path\n";   			}
if (!(-e($picard_path)))   			   { $fatal++; print "ERROR: cannot find $picard_path\n";   			}
if (!(-e($efetch_path)))   		   	   { $fatal++; print "ERROR: cannot find $efetch_path\n";   			}
if (!(-e($xtract_path)))   		   	   { $fatal++; print "ERROR: cannot find $xtract_path\n";   			}
if (!(-e($esearch_path)))   		   { $fatal++; print "ERROR: cannot find $esearch_path\n";   			}
if (!(-e($esummary_path)))   		   { $fatal++; print "ERROR: cannot find $esummary_path\n";   			}
if (!(-e($tetyper_path)))   		   { $fatal++; print "ERROR: cannot find $tetyper_path\n";   			}
if (!(-e($vcfsort_path)))   		   { $fatal++; print "ERROR: cannot find $vcfsort_path\n";   			}
if (!(-e($samtools_path))) 			   { $fatal++; print "ERROR: cannot find $samtools_path\n"; 			}
if (!(-e($sambamba_path))) 			   { $fatal++; print "ERROR: cannot find $sambamba_path\n"; 			}
if (!(-e($cutadapt_path))) 			   { $fatal++; print "ERROR: cannot find $cutadapt_path\n"; 			}
if (!(-e($trim_galore_path)))		   { $fatal++; print "ERROR: cannot find $trim_galore_path\n";			}
if (!(-e($trimmomatic_path)))		   { $fatal++; print "ERROR: cannot find $trimmomatic_path\n";			}
if (!(-e($ktImportText_path)))		   { $fatal++; print "ERROR: cannot find $ktImportText_path\n";			}
if (!(-e($aspera_openssh_path)))	   { $fatal++; print "ERROR: cannot find $aspera_openssh_path\n";		}
if (!(-e($vcfoverlay_path))) 		   { $fatal++; print "ERROR: cannot find $vcfoverlay_path\n"; 			}
if (!(-e($vcfallelicprimitives_path))) { $fatal++; print "ERROR: cannot find $vcfallelicprimitives_path\n"; }
exit if ($fatal > 0);

# PARAMETERS
my $num_threads = 40;
my $phred_threshold = 20;
my $minimum_read_length = 100;
my $maximum_pc_of_n_bases = ceil(0.05*$minimum_read_length); # i.e. 5% of the bases in a read of a given minimum length are allowed to be N
my $num_of_bases_overlapping_with_adapter_before_trimming = 5; # TrimGalore parameter -s; note the default is the very stringent value of 1
my $run_initial_fastqc = 'yes';
my $run_trimgalore     = 'yes';
my $run_trimmomatic    = 'yes';
my $run_final_fastqc   = 'yes';
my $remove_human_reads_directly_or_indirectly = 'indirectly'; # NOT YET IMPLEMENTED. 'directly' means classifying reads against a Centrifuge database of human sequences and removing any matches; 'indirectly' means we exclude any read that does NOT get classified as bacterial by the speciation step (i.e. by Kaiju using the EMBL proGenomes database)
my $for_reference_genome_pick_default_or_use_mash = 'mash'; # 'default' is the primary accession
my $remove_supplementary_alignments = 'yes'; # see SAM file specification point 1.2 ('chimeric alignment'): https://samtools.github.io/hts-specs/SAMv1.pdf
my $remove_non_primary_alignments   = 'yes'; # see SAM file specification point 1.2 ('multiple mapping'): https://samtools.github.io/hts-specs/SAMv1.pdf
my $output_all_sites   = 'no'; # for variant calling with mpileup; 'yes' to output all sites, 'no' to output only the variant sites
my $regularise_vcf     = 'yes';
my $create_assembly    = 'yes'; # if 'yes', run SPAdes on the QC'd, cleaned and classfied raw reads; if 'no', don't
my $qc_assembly        = 'yes'; # requires that $create_assembly also eq 'yes'
my $run_tetyper		   = 'no';
my $keep_bam		   = 'no';
my $keep_cleaned_fq	   = 'no';
if (($run_tetyper eq 'yes') and ($create_assembly eq 'no')) { $create_assembly = 'yes'; }
if (($create_assembly eq 'no') and ($qc_assembly eq 'yes')) { $qc_assembly 	   = 'no';  }

##########

# OBTAIN LIST OF SRA RUN ACCESSIONS
if (!(-e($in_file)))
	{ print "ERROR: cannot open $in_file\n";
	  exit 1;
	}
my %sra_run_ids = ();
open(IN,$in_file) or die $!;
while(<IN>)
	{ my $run_id = $_; chomp($run_id);
	  $sra_run_ids{$run_id}++;
	}
close(IN) or die $!;

# FOR EACH RUN ID, CREATE A FILE OF COMMANDS THAT (A) DOWNLOAD RAW SEQUENCING DATA FROM THE ENA, (B) QC/PRE-PROCESS READS, (C) PREDICT GENUS/SPECIES USING KAIJU AND A SET OF CORE PROKARYOTES, (D) CLASSIFY READS USING CENTRIFUGE, (E) COMBINE KAIJU AND CENTRIFUGE OUTPUT, ADJUDICATING CLASSIFICATIONS BASED ON LOWEST COMMON ANCESTOR, (F) DEPLETE HUMAN/UNCLASSIFABLE READS (THE LATTER DEFINED MORE OPENLY AS THOSE THAT ARE NON-MICROBIAL, I.E. POTENTIAL CONTAMINANTS, INCL. HUMAN), (G) SELECT, DOWNLOAD AND INDEX AN APPROPRIATE REFERENCE GENOME BASED ON TOP SPECIES PREDICTION, (H) ALIGN CLEANED READS TO REFERENCE GENOME, (I) CALL VARIANTS, (J) REGULARISE VCF, (K) SORT, COMPRESS AND INDEX THE FINAL VCF, (L) ASSEMBLE THE SET OF CLEANED, CLASSIFIED READS, (M) RUN TETYPER, (N) TIDY UP
my $run_ids_seen = 0; my $run_ids_total = scalar keys %sra_run_ids;
foreach my $run_id (sort keys %sra_run_ids)
	{ $run_ids_seen++;
	  next if ($run_id ne 'SRR974838');
	  next if ($run_id !~ /^\w+\d+$/);
	  
	  my $now_gmt  = strftime "%e"."_"."%b"."_"."%Y"."_"."%H"."h_"."%M"."m_"."%S"."s", gmtime; # e.g 20_Oct_2018_08h_23m_58s
	  my $log_file = "$out_dir/log.$run_id.txt";
	  open(LOG,'>',$log_file) or die $!;
	  
	  my $sh_file = "$home/sp3/commands.$run_id.sh";
	  open(SH,'>',$sh_file) or die $!;
	  print SH "#!/bin/bash\n";
	  print SH '#$ -N '; 			  print SH "$run_id\n";
	  print SH '#$ -cwd';			  print SH "\n";
	  print SH '#$ -pe sharedmem 4';  print SH "\n";
	  print SH '#$ -l h_rt=05:00:00'; print SH "\n";
	  print SH '#$ -l h_vmem=20G';    print SH "\n";
	  print SH ". /etc/profile.d/modules.sh\n";
	  print SH "module load java/jdk/1.8.0\n";
	  print SH "module load roslin/R/3.4.3\n"; # required for Quast
	  print SH "module load roslin/python/2.7.13\n";
	  print SH "export C_INCLUDE_PATH=/home/sbush/programs/include\n"; # required for mpileup
	  print SH "export LD_LIBRARY_PATH=/home/sbush/programs/lib\n"; # required for mpileup
	  print SH "PATH=\$PATH:/home/sbush/programs/bin\n"; # required for samtools
	  print SH "PATH=\$PATH:$circos_dir\n"; # required for Quast
	  print SH "export PYTHONPATH=\$PYTHONPATH:/usr/lib64/python2.7/site-packages\n"; # required for Quast (to find the location of matplotlib)
	  print SH "export PERL5LIB=/home/sbush/perl5/lib/perl5\n"; # required for Circos (it's the location of essential Perl modules; see https://stackoverflow.com/questions/1557959/how-can-i-find-out-where-a-perl-module-is-installed)
	  print SH "SECONDS=0\n"; # start the clock
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "echo \"-->\$current_date_time. beginning analysis of $run_id\n\" >> $log_file\n";
	  
 	  # RECORD THE PARAMETERS USED FOR THE PIPELINE
	  print LOG "### Parameters ###\n";
	  print LOG "minimum Phred score for quality trimming: $phred_threshold\n";
	  print LOG "minimum read length for quality trimming: $minimum_read_length\n";
	  print LOG "maximum % of ambiguous (N) bases allowed in a read of the minimum length: $maximum_pc_of_n_bases\n";
	  print LOG "number of bases overlapping with an adapter sequence before any trimming occurs? $num_of_bases_overlapping_with_adapter_before_trimming\n";
	  print LOG "run initial FastQC (recommended)? $run_initial_fastqc\n";
	  print LOG "run TrimGalore (recommended)? $run_trimgalore\n";
	  print LOG "run Trimmomatic (recommended)? $run_trimmomatic\n";
	  print LOG "run post-processing FastQC (recommended)? $run_final_fastqc\n";
	  print LOG "will the VCF record calls for all sites? $output_all_sites (if no, VCF will only record calls at variant sites)\n";
	  print LOG "regularise VCF (recommended)? $regularise_vcf\n";
	  print LOG "will a de novo assembly be made of the QC'd, cleaned and classified reads? $create_assembly\n";
	  print LOG "QC the de novo assembly (requires that it be created)? $qc_assembly\n";
	  print LOG "run TETyper? $run_tetyper\n";
	  print LOG "are we keeping the BAM? $keep_bam\n";
	  print LOG "are we keeping the cleaned fqs? $keep_cleaned_fq\n";

	  # RECORD THE VERSION NUMBERS USED FOR EACH PIECE OF SOFTWARE IN THE PIPELINE (WE KNOW THESE COMMANDS CANNOT FAIL AS WE TESTED FOR THE EXISTENCE OF EACH VARIABLE EARLIER)
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"\n-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. testing software availability before beginning...\" >> $log_file\n";
	  print SH "echo \"\n### Aspera ###\n\" >> $log_file\n";
	  print SH "$ascp_path --version >> $log_file\n";
	  print SH "echo \"\n### FastQC ###\n\" >> $log_file\n";
	  print SH "$fastqc_path --version >> $log_file\n";
	  print SH "echo \"\n### cutadapt ###\n\" >> $log_file\n";
	  print SH "$cutadapt_path --version >> $log_file\n";
	  print SH "echo \"\n### TrimGalore ###\n\" >> $log_file\n";
	  print SH "$trim_galore_path --version >> $log_file\n";
	  print SH "echo \"\n### Trimmomatic ###\n\" >> $log_file\n";
	  print SH "java -jar $trimmomatic_path -version >> $log_file\n";
	  print SH "echo \"\n### Centrifuge ###\n\" >> $log_file\n";
	  print SH "$centrifuge_path -version &>> $log_file\n";
	  print SH "echo \"\n### Kaiju ###\n\" >> $log_file\n";
	  print SH "$kaiju_path -h &>> $log_file\n";
	  print SH "echo \"\n### Mash ###\n\" >> $log_file\n";
	  print SH "$mash_path --version &>> $log_file\n";
	  print SH "echo \"\n### seqtk ###\n\" >> $log_file\n";
	  print SH "$seqtk_path &>> $log_file\n";
	  print SH "echo \"\n### Entrez Direct ###\n\" >> $log_file\n";
	  print SH "$efetch_path --help &>> $log_file\n";
	  print SH "echo \"\n### BWA ###\n\" >> $log_file\n";
	  print SH "$bwa_path &>> $log_file\n";
	  print SH "echo \"\n### SAMtools ###\n\" >> $log_file\n";
	  print SH "$samtools_path &>> $log_file\n";
	  print SH "echo \"\n### BCFtools ###\n\" >> $log_file\n";
	  print SH "$bcftools_path --version &>> $log_file\n";
	  print SH "echo \"\n### Picard Tools ###\n\" >> $log_file\n";
	  print SH "java -jar $picard_path MarkDuplicates --version &>> $log_file\n";
	  print SH "echo \"\n### VCFlib ###\n\" >> $log_file\n";
	  print SH "$vcfoverlay_path --version &>> $log_file\n"; # $vcfoverlay_path is chosen because version numbers are not available for $vcfallelicprimitives_path; see https://github.com/vcflib/vcflib/issues/102
	  print SH "echo \"\n### VCFtools ###\n\" >> $log_file\n";
	  print SH "$vcfsort_path --help &>> $log_file\n";
	  print SH "echo \"\n### bgzip ###\n\" >> $log_file\n";
	  print SH "$bgzip_path --version &>> $log_file\n";
	  print SH "echo \"\n### tabix ###\n\" >> $log_file\n";
	  print SH "$tabix_path --version &>> $log_file\n";
	  print SH "echo \"\n### SPAdes ###\n\" >> $log_file\n";
	  print SH "python $spades_path --version &>> $log_file\n";
	  print SH "echo \"\n### Quast ###\n\" >> $log_file\n";
	  print SH "python $quast_path --version &>> $log_file\n";
	  print SH "echo \"\n### Circos ###\n\" >> $log_file\n";
	  print SH "$circos_path --version &>> $log_file\n";
	  print SH "echo \"\n### TETyper ###\n\" >> $log_file\n";
	  print SH "$tetyper_path -h &>> $log_file\n";
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"\n-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. confirming database availability...\" >> $log_file\n";	  
	  print SH "echo \"\n### Database employed by Centrifuge ###\n\" >> $log_file\n";
	  print SH "date_last_modified=\"`date -r $centrifuge_db`\"\n";
	  print SH "echo \"$centrifuge_db. last modified \$date_last_modified\" >> $log_file\n";	  
	  print SH "echo \"\n### Database employed by Kaiju ###\n\" >> $log_file\n";
	  print SH "date_last_modified=\"`date -r $kaiju_db_p`\"\n";
	  print SH "echo \"$kaiju_db_p. last modified \$date_last_modified\" >> $log_file\n";

	  # START THE PIPELINE, SETTING IT TO TRAP ERRORS
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"\n-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. beginning pipeline...\n\" >> $log_file\n";
	  print SH "set -e\n";
	  print SH "trap 'last_command=\$current_command; current_command=\$BASH_COMMAND' DEBUG\n";
	  print SH "trap 'echo \"the last command run - \${last_command} - finished with exit code \$?.\"' EXIT\n";
	  
	  # (A) DOWNLOAD RAW DATA FROM ENA
	  my $first_3; my $first_6; my $digits;
	  if ($run_id =~ /^(.{3}).*?$/) { $first_3 = $1; }
	  if ($run_id =~ /^(.{6}).*?$/) { $first_6 = $1; }
	  if ($run_id =~ /^.*?(\d+)$/)  { $digits  = $1; }
	  my $number_of_digits = length($digits);
	  my $ena1 = ''; my $ena2 = '';
	  if ($number_of_digits == 6)
		{ $ena1 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$run_id/$run_id"."_1.fastq.gz";
		  $ena2 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$run_id/$run_id"."_2.fastq.gz";
		}
	  elsif ($number_of_digits == 7)
		{ if ($digits =~ /^.+?(\d{1})$/)
			{ my $last_digit = $1;
			  $ena1 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_1.fastq.gz";
			  $ena2 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_2.fastq.gz";
			}
		}
	  elsif ($number_of_digits == 8)
		{ if ($digits =~ /^.+?(\d{2})$/)
			{ my $last_two_digits = $1;
			  $ena1 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_1.fastq.gz";
			  $ena2 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_2.fastq.gz";
			}
		}
	  elsif ($number_of_digits == 9)
		{ if ($digits =~ /^.+?(\d{3})$/)
			{ my $last_three_digits = $1;
			  $ena1 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_1.fastq.gz";
			  $ena2 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_2.fastq.gz";
			}
		}
	  if (($ena1 eq '') or ($ena2 eq '2'))
		{ print LOG "WARNING: unable to parse SRA run ID $run_id; this will be skipped\n";
		}
	  next if (($ena1 eq '') or ($ena2 eq '2'));
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. downloading fastqs for $run_id...\n\" >> $log_file\n";
	  print SH "mkdir -p $out_dir/$run_id\n";
	  print SH "cd $out_dir/$run_id\n";
	  print SH "$ascp_path -QT -l 100m -P33001 -i $aspera_openssh_path $ena1 $out_dir/$run_id &>> $log_file\n"; # note that these commands only obtain PAIRED-END reads (by contrast, single-end read files take the form "$run_id.fastq.gz")
	  print SH "$ascp_path -QT -l 100m -P33001 -i $aspera_openssh_path $ena2 $out_dir/$run_id &>> $log_file\n";
	  my $fq_1 = "$out_dir/$run_id/$run_id"."_1.fastq.gz";
	  my $fq_2 = "$out_dir/$run_id/$run_id"."_2.fastq.gz";
	  my $log_number = 1;
	  
	  # (B) QUALITY INSPECTION AND PRE-PROCESSING
	  # initial quality inspection with FastQC
  	  if ($run_initial_fastqc eq 'yes')
		{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. running initial FastQC on $run_id...\n\" >> $log_file\n";
		  print SH "$fastqc_path $fq_1 $fq_2 &>> $log_file\n";
		  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_initial_fastqc_report\n";
		  print SH "mv $out_dir/$run_id/$run_id"."_1_fastqc.html $out_dir/$run_id/logs/$log_number"."_initial_fastqc_report\n";
		  print SH "mv $out_dir/$run_id/$run_id"."_2_fastqc.html $out_dir/$run_id/logs/$log_number"."_initial_fastqc_report\n";
		  print SH "mv $out_dir/$run_id/$run_id"."_1_fastqc.zip $out_dir/$run_id/logs/$log_number"."_initial_fastqc_report\n";
		  print SH "mv $out_dir/$run_id/$run_id"."_2_fastqc.zip $out_dir/$run_id/logs/$log_number"."_initial_fastqc_report\n";
		  $log_number++;
		}
	  # first round of cleaning (to automatically predict, and then remove, adapter sequence): TrimGalore
	  # note that Trimmomatic and cutadapt have complementary, but distinct, functions. While Trimmomatic can also remove adapters, it requires that the adapter sequence be specified. By contrast, TrimGalore (which employs cutadapt) automatically predicts the adapter - more useful when input is of unknown origin. TrimGalore can also filter reads based on N content: Trimmomatic cannot. Here we use TrimGalore for N and adapter removal only, with quality-trimming later performed by Trimmomatic.
	  if ($run_trimgalore eq 'yes')
		{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. running cutadapt/TrimGalore on $run_id...\n\" >> $log_file\n";
		  print SH "$trim_galore_path --path_to_cutadapt $cutadapt_path --stringency $num_of_bases_overlapping_with_adapter_before_trimming --trim-n --max_n $maximum_pc_of_n_bases --paired --length 0 -o $out_dir/$run_id $fq_1 $fq_2 &>> $log_file\n";
		  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_trimgalore_report\n";
		  print SH "mv $fq_1"."_trimming_report.txt $out_dir/$run_id/logs/$log_number"."_trimgalore_report\n";
		  print SH "mv $fq_2"."_trimming_report.txt $out_dir/$run_id/logs/$log_number"."_trimgalore_report\n";
		  my $trimmed_fq_1 = "$out_dir/$run_id/$run_id"."_1_val_1.fq.gz";
		  my $trimmed_fq_2 = "$out_dir/$run_id/$run_id"."_2_val_2.fq.gz";
		  print SH "mv $trimmed_fq_1 $fq_1\n";
		  print SH "mv $trimmed_fq_2 $fq_2\n";
		  $log_number++;
		}
	  # second round of cleaning (for general read quality): Trimmomatic. This will (i) remove bases from the end of a read if they are below a Phred score of $phred_threshold, (ii) clip the read if the average Phred score within a 4bp sliding window advanced from the 5’ end falls below $phred_threshold, and (iii) impose a minimum read length of $minimum_read_length bp
	  # given 2 fqs as input, Trimmomatic will output 4 files, 2 of which are for singleton reads: read 1 where read 2 failed, and read 2 where read 1 failed. These are combined into one file of de facto single-end reads, which may be aligned separately.
	  my $fq_se = '';
	  if ($run_trimmomatic eq 'yes')
		{ my $run_id_fwd_paired = "$out_dir/$run_id/$run_id.forward_paired.fq.gz";
		  my $run_id_fwd_single = "$out_dir/$run_id/$run_id.forward_unpaired.fq.gz";
		  my $run_id_rev_paired = "$out_dir/$run_id/$run_id.reverse_paired.fq.gz";
		  my $run_id_rev_single = "$out_dir/$run_id/$run_id.reverse_unpaired.fq.gz";
		  $fq_se    	    	= "$out_dir/$run_id/$run_id.se.fq.gz";
		  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_trimmomatic_report\n";
		  my $trim_log			= "$out_dir/$run_id/logs/$log_number"."_trimmomatic_report/$run_id.trimmomatic_log.txt";
		  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. running Trimmomatic on the adapter-trimmed reads from $run_id...\n\" >> $log_file\n";
		  print SH "java -jar $trimmomatic_path PE -phred33 -threads $num_threads $fq_1 $fq_2 $run_id_fwd_paired $run_id_fwd_single $run_id_rev_paired $run_id_rev_single TRAILING:$phred_threshold SLIDINGWINDOW:4:$phred_threshold MINLEN:$minimum_read_length >$trim_log 2>&1\n";
		  print SH "mv $run_id_fwd_paired $fq_1\n";
		  print SH "mv $run_id_rev_paired $fq_2\n";
		  print SH "zcat $run_id_fwd_single $run_id_rev_single | gzip > $fq_se\n";
		  print SH "rm $run_id_fwd_single $run_id_rev_single\n";
		  $log_number++;
		}
	  # final quality inspection with FastQC
	  if ($run_final_fastqc eq 'yes')
		{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. running FastQC on the cleaned and adapter-trimmed reads from $run_id...\n\" >> $log_file\n";
		  if ($fq_se ne '')
			{ print SH "$fastqc_path $fq_1 $fq_2 $fq_se &>> $log_file\n"; }
		  elsif ($fq_se eq '')
			{ print SH "$fastqc_path $fq_1 $fq_2 &>> $log_file\n"; }
		  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_final_fastqc_report\n";
		  print SH "mv $out_dir/$run_id/$run_id"."_1_fastqc.html $out_dir/$run_id/logs/$log_number"."_final_fastqc_report\n";
		  print SH "mv $out_dir/$run_id/$run_id"."_2_fastqc.html $out_dir/$run_id/logs/$log_number"."_final_fastqc_report\n";
		  print SH "mv $out_dir/$run_id/$run_id"."_1_fastqc.zip $out_dir/$run_id/logs/$log_number"."_final_fastqc_report\n";
		  print SH "mv $out_dir/$run_id/$run_id"."_2_fastqc.zip $out_dir/$run_id/logs/$log_number"."_final_fastqc_report\n";
		  print SH "mv $out_dir/$run_id/$run_id.se_fastqc.html $out_dir/$run_id/logs/$log_number"."_final_fastqc_report\n";
		  print SH "mv $out_dir/$run_id/$run_id.se_fastqc.zip $out_dir/$run_id/logs/$log_number"."_final_fastqc_report\n";
		  $log_number++;
		}
	  
	  # (C) CREATE KRONA-FORMATTED KAIJU REPORTS FOR THE PE AND SE READS, ALONGSIDE PREDICTING GENUS AND SPECIES AND REPORTING A LIST OF CLASSIFIED READS (SO AS TO BE ABLE TO EXCLUDE THE UNCLASSIFIED, IF DESIRED)
	  my $kaiju_report_log_number = $log_number;
	  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_kaiju_report\n";
	  for(my $x=0;$x<=1;$x++)
		{ my $pe_or_se = '';
		  if ($x == 0) { $pe_or_se = 'paired-end'; } elsif ($x == 1) { $pe_or_se = 'single-end'; }
		  next if (($pe_or_se eq 'single-end') and ($run_trimmomatic eq 'no'));
		  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. running Kaiju on the cleaned and adapter-trimmed $pe_or_se reads from $run_id...\n\" >> $log_file\n";
		  if ($pe_or_se eq 'paired-end')
			{ print SH "$kaiju_path -z $num_threads -t $nodes_dmp_p -f $kaiju_db_p -i $fq_1 -j $fq_2 -a greedy -e 5 -E 0.05 -o $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.kaiju_report.$pe_or_se.txt\n"; }
		  elsif ($pe_or_se eq 'single-end')
			{ print SH "$kaiju_path -z $num_threads -t $nodes_dmp_p -f $kaiju_db_p -i $fq_se -a greedy -e 5 -E 0.05 -o $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.kaiju_report.$pe_or_se.txt\n"; }
		  print SH "$kaijuReport_path -t $nodes_dmp_p -n $names_dmp_p -i $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.kaiju_report.$pe_or_se.txt -r genus -o $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.kaiju_genus_prediction.$pe_or_se.txt\n";
		  print SH "$kaijuReport_path -t $nodes_dmp_p -n $names_dmp_p -i $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.kaiju_report.$pe_or_se.txt -r species -o $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.kaiju_species_prediction.$pe_or_se.txt\n";
		  print SH "awk '\$1==\"C\" { print \$2 }' $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.kaiju_report.$pe_or_se.txt > $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.classified_read_list.$pe_or_se.txt\n";
		  print SH "$kaijutokrona_path -t $nodes_dmp_p -n $names_dmp_p -i $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.kaiju_report.$pe_or_se.txt -o $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.kaiju_report_for_krona.$pe_or_se.txt &>> $log_file\n";
		  print SH "rm $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.kaiju_report.$pe_or_se.txt\n";
		  print SH "$ktImportText_path -o $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.kaiju_report.$pe_or_se.html $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.kaiju_report_for_krona.$pe_or_se.txt &>> $log_file\n";
		  print SH "rm $out_dir/$run_id/logs/$log_number"."_kaiju_report/$run_id.kaiju_report_for_krona.$pe_or_se.txt\n";
		}
	  $log_number++;
		
=cut
	  # (D) CREATE CENTRIFUGE REPORT FOR THE PE AND SE READS
	  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_centrifuge_report\n";
	  for(my $x=0;$x<=1;$x++)
		{ my $pe_or_se = '';
		  if ($x == 0) { $pe_or_se = 'paired-end'; } elsif ($x == 1) { $pe_or_se = 'single-end'; }
		  next if (($pe_or_se eq 'single-end') and ($run_trimmomatic eq 'no'));
		  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. running Centrifuge on the cleaned and adapter-trimmed $pe_or_se reads from $run_id...\n\" >> $log_file\n";
		  if ($pe_or_se eq 'paired-end')
			{ print SH "$centrifuge_path -x $centrifuge_index_prefix -m1 $fq_1 -m2 $fq_2\n"; }
		  elsif ($pe_or_se eq 'single-end')
			{ print SH "$centrifuge_path -x $centrifuge_index_prefix -r $fq_se\n"; }
		  $log_number++;
		}
		
	  # (E) COMBINE KAIJU AND CENTRIFUGE OUTPUT, ADJUDICATING CLASSIFICATIONS BASED ON LOWEST COMMON ANCESTOR
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. combining Centrifuge and Kaiju reports for $run_id...\n\" >> $log_file\n";
	  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_consensus_kaiju_and_centrifuge_report\n";
	  print SH "$mergeOutputs_path -t $nodes_dmp -c lca -i <(sort -k2,2 centrifuge.out) -j <(sort -k2,2 centrifuge.out) -o combined.out\n";
	  $log_number++;
=cut	
	  # (F) REMOVE UNCLASSIFIABLE READS EITHER DIRECTLY (I.E. THOSE CLASSIFIED AS HUMAN) OR INDIRECTLY (I.E. THOSE NOT CLASSIFIED AS BACTERIAL), USING CENTRIFUGE OR KAIJU, RESPECTIVELY
	  if ($remove_human_reads_directly_or_indirectly eq 'indirectly')
		{ my $fq_1_excl_unclassifiable  = "$out_dir/$run_id/$run_id"."_1.excl_unclassifiable.fastq.gz";
		  my $fq_2_excl_unclassifiable  = "$out_dir/$run_id/$run_id"."_2.excl_unclassifiable.fastq.gz";
		  my $fq_se_excl_unclassifiable = "$out_dir/$run_id/$run_id.se.excl_unclassifiable.fastq.gz";
		  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. removing reads unclassified by Kaiju from the cleaned and adapter-trimmed paired-end reads from $run_id...\n\" >> $log_file\n";
		  print SH "$seqtk_path subseq $fq_1 $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_report/$run_id.classified_read_list.paired-end.txt | gzip > $fq_1_excl_unclassifiable\n";
		  print SH "$seqtk_path subseq $fq_2 $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_report/$run_id.classified_read_list.paired-end.txt | gzip > $fq_2_excl_unclassifiable\n";
		  if ($fq_se ne '')
			{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
			  print SH "elapsed_time=\$SECONDS\n";
			  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. removing reads unclassified by Kaiju from the cleaned and adapter-trimmed single-end reads from $run_id...\n\" >> $log_file\n";
			  print SH "$seqtk_path subseq $fq_se $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_report/$run_id.classified_read_list.single-end.txt > $fq_se_excl_unclassifiable\n";
			}
		  print SH "rm $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_report/$run_id.classified_read_list.paired-end.txt $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_report/$run_id.classified_read_list.single-end.txt\n";
		  print SH "mv $fq_1_excl_unclassifiable $fq_1\n";
		  print SH "mv $fq_2_excl_unclassifiable $fq_2\n";
		  if ($fq_se ne '')
			{ print SH "mv $fq_se_excl_unclassifiable $fq_se\n"; }
		}
		
	  # (G) SELECT, DOWNLOAD AND INDEX AN APPROPRIATE REFERENCE GENOME. BY DEFAULT, THIS WILL BE THE NCBI DESIGNATED PRIMARY ACCESSION (SEE https://www.biostars.org/p/244479/)
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. downloading reference genome...\n\" >> $log_file\n";
	  print SH "top_species_hit_pe=\$(awk -F '\\t' 'NR==3 {print \$3}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_report/$run_id.kaiju_species_prediction.paired-end.txt)\n";
	  print SH "top_species_hit_se=\$(awk -F '\\t' 'NR==3 {print \$3}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_report/$run_id.kaiju_species_prediction.single-end.txt)\n";
	  print SH "if [ \"\$top_species_hit_pe\" != \"\$top_species_hit_se\" ]; then echo \"ERROR: candidate species for $run_id as predicted by top hit for PE reads (\$top_species_hit_pe) != top hit for SE reads (\$top_species_hit_se); unable to proceed with alignment\n\" >> $log_file; fi\n";
	  print SH "if [ \"\$top_species_hit_pe\" != \"\$top_species_hit_se\" ]; then exit 1; fi\n"; # CHECKPOINT: fail if there are discordant Kaiju predictions
	  print SH "GENOME=\$($esearch_path -db genome -query \"\$top_species_hit_pe\"[orgn] | $efetch_path -format docsum | tee \"$out_dir/$run_id/\$top_species_hit_pe.genome.esearch.docsum\")\n";
	  print SH "ACC=`echo \$GENOME | $xtract_path -pattern DocumentSummary  -element Assembly_Accession`\n";
	  print SH "RESULT=\$($esearch_path -db assembly -query \"\$ACC\" | $efetch_path -format docsum | tee \"$out_dir/$run_id/\$top_species_hit_pe.assembly.esearch.docsum\")\n";
	  print SH "FTPP=`echo \$RESULT | $xtract_path -pattern DocumentSummary -element FtpPath_GenBank`\n";
	  print SH "TAXID=`echo \$RESULT | $xtract_path -pattern DocumentSummary -element Taxid`\n";
	  print SH "SPECIESTAXID=`echo \$RESULT | $xtract_path -pattern DocumentSummary -element SpeciesTaxid`\n";
	  print SH "BASENAME=`basename \$FTPP`\n";
	  print SH "FTPPATHG=\$FTPP/\$BASENAME'_genomic.fna.gz'\n";
	  print SH "FTPPATHGFF=\$FTPP/\$BASENAME'_genomic.gff.gz'\n";
	  print SH "rm \"$out_dir/$run_id/\$top_species_hit_pe.genome.esearch.docsum\" \"$out_dir/$run_id/\$top_species_hit_pe.assembly.esearch.docsum\"\n";
	  if ($for_reference_genome_pick_default_or_use_mash eq 'default')
		{ print SH "wget \$FTPPATHG -O $out_dir/$run_id/ref.fa.gz &>> $log_file\n";
		  print SH "gunzip $out_dir/$run_id/ref.fa.gz\n";
		  print SH "wget \$FTPPATHGFF -O $out_dir/$run_id/ref.gff.gz &>> $log_file\n";
		  print SH "gunzip $out_dir/$run_id/ref.gff.gz\n";
		}
	  print SH "echo \"-->species predicted: \$top_species_hit_pe\n\" >> $log_file\n";
	  if ($for_reference_genome_pick_default_or_use_mash eq 'default')
		{ print SH "echo \"-->NCBI taxonomy ID (strain): \$TAXID\n\" >> $log_file\n"; }
	  print SH "echo \"-->NCBI taxonomy ID (species): \$SPECIESTAXID\n\" >> $log_file\n";
	  if ($for_reference_genome_pick_default_or_use_mash eq 'default')
		{ print SH "echo \"-->reference genome downloaded: \$FTPPATHG\n\" >> $log_file\n"; }
	  my $ref = "$out_dir/$run_id/ref.fa";
	  my $gff = "$out_dir/$run_id/ref.gff";
	  
	  # USE THE KAIJU BEST SPECIES HIT TO NARROW DOWN A SEARCH OF "SPECIES SPACE". WE DO THIS BY DOWNLOADING THE LATEST COMPLETE GENOMES FOR EACH SPECIES, CREATING MASH SKETCHES FOR THEM AND THEN PICKING THE SKETCH WITH SHORTEST DISTANCE TO THE READS.
	  # see question 15 ("how can I download RefSeq data for all complete bacterial genomes?") of this FAQ: https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/
	  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_reference_genomes_downloaded\n";
	  my $ref_genome_download_log = "$out_dir/$run_id/logs/$log_number"."_reference_genomes_downloaded/$run_id.reference_genome_download_log.txt";
	  print SH "echo \"TOP KAIJU HIT (SPECIES): \$top_species_hit_pe\" >> $ref_genome_download_log\n";
	  print SH "echo \"NCBI TAXON ID (FOR SPECIES): \$SPECIESTAXID\" >> $ref_genome_download_log\n";
	  if ($for_reference_genome_pick_default_or_use_mash eq 'default')
		{ print SH "echo \"NCBI TAXON ID (FOR PRIMARY ACCESSION OF THIS SPECIES): \$TAXID\" >> $ref_genome_download_log\n";
		  print SH "echo \"URL FOR PRIMARY ACCESSION: \$FTPPATHG\" >> $ref_genome_download_log\n";
		}
	  elsif ($for_reference_genome_pick_default_or_use_mash eq 'mash')
		{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. downloading the set of complete genomes corresponding to the best hit taxonomy ID...\n\" >> $log_file\n";
		  print SH "wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt &>> $log_file\n";
		  print SH "date_last_modified=\"`date -r assembly_summary.txt`\"\n";
		  print SH "echo \"--> downloaded ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt. last modified \$date_last_modified\" >> $log_file\n";
		  print SH "awk -v speciestaxid=\$SPECIESTAXID -F \"\\t\" '\$7==speciestaxid && \$12==\"Complete Genome\" && \$11==\"latest\"{print \$20}' $out_dir/$run_id/assembly_summary.txt > $out_dir/$run_id/ftpdirpaths.txt\n";
		  print SH "awk 'BEGIN{FS=OFS=\"/\";filesuffix=\"genomic.fna.gz\"}{ftpdir=\$0;asm=\$10;file=asm\"_\"filesuffix;print ftpdir,file}' $out_dir/$run_id/ftpdirpaths.txt > $out_dir/$run_id/ftpfilepaths.txt\n";
		  print SH "rm $out_dir/$run_id/ftpdirpaths.txt\n";
		  print SH "number_of_complete_genomes=\"`wc -l < $out_dir/$run_id/ftpfilepaths.txt`\"\n";
		  print SH "echo \"--> assembly_summary.txt parsed to identify \$number_of_complete_genomes complete genomes for species taxon ID \$SPECIESTAXID\" >> $log_file\n";
		  print SH "echo \"NUMBER OF COMPLETE GENOMES ASSIGNED THIS SPECIES TAXON ID, OBTAINED AFTER PARSING ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt (DOWNLOADED \$date_last_modified): \$number_of_complete_genomes\" >> $ref_genome_download_log\n";
		  print SH "echo \"GENOME URLS:\" >> $ref_genome_download_log\n";
		  print SH "cat $out_dir/$run_id/ftpfilepaths.txt >> $ref_genome_download_log\n";
		  print SH "mkdir $out_dir/$run_id/temp_dir_genome_storage\n";
		  print SH "wget -i $out_dir/$run_id/ftpfilepaths.txt -P $out_dir/$run_id/temp_dir_genome_storage &>> $log_file\n";
		  print SH "rm $out_dir/$run_id/ftpfilepaths.txt\n";
		  print SH "find $out_dir/$run_id/temp_dir_genome_storage -maxdepth 1 -type f > $out_dir/$run_id/fa_filepaths_for_mash_sketch.txt\n";
		  if ($fq_se ne '')
			{ print SH "cat $fq_1 $fq_2 $fq_se > $out_dir/$run_id/$run_id.concatenated_reads.fq.gz\n"; }
		  elsif ($fq_se eq '')
			{ print SH "cat $fq_1 $fq_2 > $out_dir/$run_id/$run_id.concatenated_reads.fq.gz\n"; }
		  print SH "$mash_path sketch -p $num_threads -l $out_dir/$run_id/fa_filepaths_for_mash_sketch.txt -o $out_dir/$run_id/$run_id.fa &>> $log_file\n";
		  print SH "$mash_path sketch -p $num_threads -o $out_dir/$run_id/$run_id.fq $out_dir/$run_id/$run_id.concatenated_reads.fq.gz &>> $log_file\n";
		  print SH "rm $out_dir/$run_id/$run_id.concatenated_reads.fq.gz\n";
		  print SH "$mash_path dist -p $num_threads $out_dir/$run_id/$run_id.fq.msh $out_dir/$run_id/$run_id.fa.msh > $out_dir/$run_id/$run_id.mash-dist_output.tsv\n";
		  print SH "closest_mash_hit=\$(basename \$(sort -g -k3 $out_dir/$run_id/$run_id.mash-dist_output.tsv | head -1 | cut -f2))\n";
		  print SH "closest_mash_hit_root=\${closest_mash_hit%_genomic.fna.gz}\n";
		  print SH "path_to_closest_mash_hit=\$(grep \$closest_mash_hit_root $out_dir/$run_id/assembly_summary.txt | head -1 | awk '{print \$(NF)}')\n";
		  print SH "FTPPATHG=\$path_to_closest_mash_hit/\$closest_mash_hit_root'_genomic.fna.gz'\n";
		  print SH "FTPPATHGFF=\$path_to_closest_mash_hit/\$closest_mash_hit_root'_genomic.gff.gz'\n";
		  print SH "wget \$FTPPATHG -O $out_dir/$run_id/ref.fa.gz &>> $log_file\n";
		  print SH "gunzip $out_dir/$run_id/ref.fa.gz\n";
		  print SH "wget \$FTPPATHGFF -O $out_dir/$run_id/ref.gff.gz &>> $log_file\n";
		  print SH "gunzip $out_dir/$run_id/ref.gff.gz\n";
		  print SH "echo \"URL FOR CLOSEST MASH HIT: \$FTPPATHG\" >> $ref_genome_download_log\n";
		  print SH "mv $out_dir/$run_id/$run_id.mash-dist_output.tsv $out_dir/$run_id/logs/$log_number"."_reference_genomes_downloaded/$run_id.mash-dist_output.tsv\n";
		  print SH "rm $out_dir/$run_id/assembly_summary.txt\n";
		  print SH "rm $out_dir/$run_id/fa_filepaths_for_mash_sketch.txt\n";
		  print SH "rm $out_dir/$run_id/$run_id.fq.msh $out_dir/$run_id/$run_id.fa.msh\n";
		  print SH "rm -r $out_dir/$run_id/temp_dir_genome_storage\n";		  
		}
	  $log_number++;
	  
	  # INDEX THE REFERENCE GENOME
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. indexing reference genome...\n\" >> $log_file\n";
	  print SH "$bwa_path index $ref 2>> $log_file\n";
	  
	  # (H) ALIGN CLEANED, CLASSIFIED, READS TO REFERENCE GENOME USING BWA. THE BAM IS POST-PROCESSED WITH PICARD TO COORDINATE-SORT, FIX MATE INFORMATION (IF RELEVANT) AND TO MARK DUPLICATES.
	  for(my $x=0;$x<=1;$x++) # this loop aligns, in turn, the set of PE and SE reads (those which, after the above trimming/QC, have lost their mate)
		{ my $pe_or_se = '';
		  if ($x == 0) { $pe_or_se = 'paired-end'; } elsif ($x == 1) { $pe_or_se = 'single-end'; }
		  next if (($pe_or_se eq 'single-end') and ($run_trimmomatic eq 'no'));
		  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. aligning the quality- and adapter-trimmed $pe_or_se reads from $run_id...\n\" >> $log_file\n";
		  if ($pe_or_se eq 'paired-end')
			{ #print SH "$bwa_path mem -R '\@RG\\tID:group_$run_id\\tSM:sample_$run_id\\tPL:Illumina\\tLIB:lib_$run_id\\tPU:unit_$run_id' -t $num_threads -M $ref $fq_1 $fq_2 2>> $log_file | $samtools_path view -Shb - > $out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam\n";
			  print SH "$bwa_path mem -R '\@RG\\tID:group_$run_id\\tSM:sample_$run_id\\tPL:Illumina\\tLIB:lib_$run_id\\tPU:unit_$run_id' -t $num_threads -M $ref $fq_1 $fq_2 2>> $log_file | $sambamba_path view -S -f bam -t $num_threads -o $out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam /dev/stdin\n";
			  if (($create_assembly eq 'no') and ($keep_cleaned_fq eq 'no')) { print SH "rm $fq_1 $fq_2\n"; }
			}
		  elsif ($pe_or_se eq 'single-end')
			{ #print SH "$bwa_path mem -R '\@RG\\tID:group_$run_id\\tSM:sample_$run_id\\tPL:Illumina\\tLIB:lib_$run_id\\tPU:unit_$run_id' -t $num_threads -M $ref $fq_se 2>> $log_file | $samtools_path view -Shb - > $out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam\n";
			  print SH "$bwa_path mem -R '\@RG\\tID:group_$run_id\\tSM:sample_$run_id\\tPL:Illumina\\tLIB:lib_$run_id\\tPU:unit_$run_id' -t $num_threads -M $ref $fq_se 2>> $log_file | $sambamba_path view -S -f bam -t $num_threads -o $out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam /dev/stdin\n";
			  if (($create_assembly eq 'no') and ($keep_cleaned_fq eq 'no')) { print SH "rm $fq_se\n"; }
			}
		  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. cleaning BAM of the aligned, classified, quality- and adapter-trimmed $pe_or_se reads from $run_id...\n\" >> $log_file\n";
		  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam OUTPUT=$out_dir/$run_id/$run_id.unsorted.cleaned.$pe_or_se.bam TMP_DIR=$out_dir/$run_id &>> $log_file\n";
		  print SH "rm $out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam\n";
		  print SH "mv $out_dir/$run_id/$run_id.unsorted.cleaned.$pe_or_se.bam $out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam\n";
		  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. fixing mate information in the BAM of the aligned, classified, quality- and adapter-trimmed $pe_or_se reads from $run_id...\n\" >> $log_file\n";
		  print SH "java -jar $picard_path FixMateInformation INPUT=$out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam OUTPUT=$out_dir/$run_id/$run_id.unsorted.fixmated.$pe_or_se.bam TMP_DIR=$out_dir/$run_id &>> $log_file\n";
		  print SH "rm $out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam\n";
		  print SH "mv $out_dir/$run_id/$run_id.unsorted.fixmated.$pe_or_se.bam $out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam\n";		  
		  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. sorting BAM of the aligned, classified, quality- and adapter-trimmed $pe_or_se reads from $run_id...\n\" >> $log_file\n";
		  print SH "java -jar $picard_path SortSam INPUT=$out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam OUTPUT=$out_dir/$run_id/$run_id.sorted.$pe_or_se.bam SORT_ORDER=coordinate TMP_DIR=$out_dir/$run_id &>> $log_file\n";
		  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. duplicate-marking BAM of the aligned, classified, quality- and adapter-trimmed $pe_or_se reads from $run_id...\n\" >> $log_file\n";
		  print SH "java -jar $picard_path MarkDuplicates INPUT=$out_dir/$run_id/$run_id.sorted.$pe_or_se.bam OUTPUT=$out_dir/$run_id/$run_id.$pe_or_se.bam METRICS_FILE=$out_dir/$run_id/$run_id.$pe_or_se.metrics ASSUME_SORTED=true TMP_DIR=$out_dir/$run_id &>> $log_file\n";
		  if (($remove_supplementary_alignments eq 'yes') and ($remove_non_primary_alignments eq 'yes'))
			{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
			  print SH "elapsed_time=\$SECONDS\n";
			  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. removing supplementary and non-primary alignments from the BAM of the aligned, classified, quality- and adapter-trimmed $pe_or_se reads from $run_id...\n\" >> $log_file\n";
			  print SH "$samtools_path view -b -F 2048 -F 256 $out_dir/$run_id/$run_id.$pe_or_se.bam > $out_dir/$run_id/$run_id.$pe_or_se.rm_supplementary.bam 2>> $log_file\n";
			  print SH "mv $out_dir/$run_id/$run_id.$pe_or_se.rm_supplementary.bam $out_dir/$run_id/$run_id.$pe_or_se.bam\n";
			}
		  elsif (($remove_supplementary_alignments eq 'yes') and ($remove_non_primary_alignments eq 'no'))
			{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
			  print SH "elapsed_time=\$SECONDS\n";
			  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. removing supplementary alignments from the BAM of the aligned, classified, quality- and adapter-trimmed $pe_or_se reads from $run_id...\n\" >> $log_file\n";
			  print SH "$samtools_path view -b -F 2048 $out_dir/$run_id/$run_id.$pe_or_se.bam > $out_dir/$run_id/$run_id.$pe_or_se.rm_supplementary.bam 2>> $log_file\n";
			  print SH "mv $out_dir/$run_id/$run_id.$pe_or_se.rm_supplementary.bam $out_dir/$run_id/$run_id.$pe_or_se.bam\n";
			}
		  elsif (($remove_supplementary_alignments eq 'no') and ($remove_non_primary_alignments eq 'yes'))
			{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
			  print SH "elapsed_time=\$SECONDS\n";
			  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. removing non-primary alignments from the BAM of the aligned, classified, quality- and adapter-trimmed $pe_or_se reads from $run_id...\n\" >> $log_file\n";
			  print SH "$samtools_path view -b -F 256 $out_dir/$run_id/$run_id.$pe_or_se.bam > $out_dir/$run_id/$run_id.$pe_or_se.rm_primary.bam 2>> $log_file\n";
			  print SH "mv $out_dir/$run_id/$run_id.$pe_or_se.rm_primary.bam $out_dir/$run_id/$run_id.$pe_or_se.bam\n";
			}
		  print SH "java -jar $picard_path BuildBamIndex INPUT=$out_dir/$run_id/$run_id.$pe_or_se.bam &>> $log_file\n";
		  print SH "rm $out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam $out_dir/$run_id/$run_id.sorted.$pe_or_se.bam $out_dir/$run_id/$run_id.$pe_or_se.metrics\n";
		}
	  # merge the PE and SE BAMs
	  if ($run_trimmomatic eq 'yes')
		{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. merging PE and SE BAMs for the aligned, classified, quality- and adapter-trimmed reads from $run_id...\n\" >> $log_file\n";
		  print SH "java -jar $picard_path MergeSamFiles SORT_ORDER='coordinate' INPUT=$out_dir/$run_id/$run_id.paired-end.bam INPUT=$out_dir/$run_id/$run_id.single-end.bam OUTPUT=$out_dir/$run_id/$run_id.bam TMP_DIR=$out_dir/$run_id &>> $log_file\n";
		  print SH "java -jar $picard_path BuildBamIndex INPUT=$out_dir/$run_id/$run_id.bam &>> $log_file\n";
		  print SH "rm $out_dir/$run_id/$run_id.paired-end.bam $out_dir/$run_id/$run_id.paired-end.bai $out_dir/$run_id/$run_id.single-end.bam $out_dir/$run_id/$run_id.single-end.bai\n";
		}
	  elsif ($run_trimmomatic eq 'no')
		{ print SH "mv $out_dir/$run_id/$run_id.paired-end.bam $out_dir/$run_id/$run_id.bam\n";
		  print SH "java -jar $picard_path BuildBamIndex INPUT=$out_dir/$run_id/$run_id.bam &>> $log_file\n";
		}
	  # obtain summary statistics for the final, merged, BAM
	  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_bam_statistics\n";
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. obtaining summary statistics for the merged BAM of $run_id...\n\" >> $log_file\n";
	  print SH "$samtools_path stats $out_dir/$run_id/$run_id.bam > $out_dir/$run_id/logs/$log_number"."_bam_statistics/$run_id.samtools_stats_output.txt\n";
	  print SH "$samtools_path flagstat $out_dir/$run_id/$run_id.bam > $out_dir/$run_id/logs/$log_number"."_bam_statistics/$run_id.samtools_flagstat_output.txt\n";
	  $log_number++;
	  
	  # (I) CALL VARIANTS
	  # variant calling using mpileup
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. calling variants from $run_id using mpileup...\n\" >> $log_file\n";
	  if ($output_all_sites eq 'no') # if no, we are only interested in calling the variant sites
		{ print SH "$bcftools_path mpileup -Ou -f $ref $out_dir/$run_id/$run_id.bam 2>> $log_file | $bcftools_path call --threads $num_threads --ploidy 1 -mv -Ov -o $out_dir/$run_id/$run_id.vcf\n"; } # the "v" in "-mv" specifies that the output file only contain variant sites; omit this to output all sites
	  elsif ($output_all_sites eq 'yes')
		{ print SH "$bcftools_path mpileup -Ou -f $ref $out_dir/$run_id/$run_id.bam 2>> $log_file | $bcftools_path call --threads $num_threads --ploidy 1 -m -Ov -o $out_dir/$run_id/$run_id.vcf\n"; }
      # variant calling using cortex (which first requires that the reference genome be indexed with Stampy)
#	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
#	  print SH "elapsed_time=\$SECONDS\n";
#	  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. calling variants from $run_id using cortex...\n\" >> $log_file\n";
#	  print SH "$stampy_path -G ref $ref 2>> $log_file\n";
#	  print SH "$stampy_path -g ref -H ref 2>> $log_file\n";
#	  print SH "echo \"$ref\" > $out_dir/$run_id/cortex.in.index_ref.fofn\n";
#	  print SH "echo -e \"sample\t$out_dir/$run_id/cortex.in.index_ref.fofn\t.\t.\" > $out_dir/$run_id/cortex.in.index\n";
	  
	  # CLOCKWORK: USE MINOS TO ADJUDICATE BETWEEN VCFs PRODUCED BY BWA/MPILEUP AND BWA/CORTEX, SO CREATING ONE FINAL VCF
	  
	  # (J) REGULARISE VCF SO THAT DIFFERENT REPRESENTATIONS OF THE SAME INDEL OR COMPLEX VARIANT ARE NOT COUNTED AS DIFFERENT VARIANTS (THIS IS USEFUL IF LATER COMPARING VCFs PRODUCED BY DIFFERENT TOOLS)
	  if ($regularise_vcf eq 'yes')
		{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. regularising VCF for $run_id...\n\" >> $log_file\n";
		  print SH "$vcfallelicprimitives_path $out_dir/$run_id/$run_id.vcf > $out_dir/$run_id/$run_id.regularised.vcf 2>> $log_file\n";
		  print SH "rm --interactive=never $out_dir/$run_id/$run_id.vcf\n";
		  print SH "mv $out_dir/$run_id/$run_id.regularised.vcf $out_dir/$run_id/$run_id.vcf\n";
		}
	  
	  # (K) SORT, COMPRESS AND INDEX THE FINAL VCF
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. sorting, compressing and indexing VCF for $run_id...\n\" >> $log_file\n";
	  print SH "$vcfsort_path $out_dir/$run_id/$run_id.vcf > $out_dir/$run_id/$run_id.sorted.vcf 2>> $log_file\n";
	  print SH "mv $out_dir/$run_id/$run_id.sorted.vcf $out_dir/$run_id/$run_id.vcf\n";
	  print SH "$bgzip_path $out_dir/$run_id/$run_id.vcf\n";
	  print SH "$tabix_path -p vcf $out_dir/$run_id/$run_id.vcf.gz\n";
	  
	  # (L) ASSEMBLE THE SET OF CLEANED, CLASSIFIED READS, AND QC THE RESULTING SET OF SCAFFOLDS
	  ## note that the SPAdes parameters listed here are those also recommended for use by TEtyper, with the exception of --disable-rr (not used here as it prevents the creation of scaffolds.fasta). These parameters are also used by a previous study (https://www.nature.com/articles/s41598-018-20384-3).
	  if ($create_assembly eq 'yes')
		{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. creating de novo assembly for $run_id...\n\" >> $log_file\n";
		  if ($fq_se ne '')
			{ print SH "python $spades_path --cov-cutoff auto --careful -1 $fq_1 -2 $fq_2 -s $fq_se -t $num_threads -o $out_dir/$run_id/spades &>> $log_file\n"; }
		  elsif ($fq_se eq '')
			{ print SH "python $spades_path --cov-cutoff auto --careful -1 $fq_1 -2 $fq_2 -t $num_threads -o $out_dir/$run_id/spades &>> $log_file\n"; }
		  print SH "mv $out_dir/$run_id/spades/scaffolds.fasta $out_dir/$run_id/$run_id.fa\n";
		  print SH "rm -r $out_dir/$run_id/spades\n";
		  # create a mash sketch of the scaffolds
		  if ($qc_assembly eq 'yes')
			{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
			  print SH "elapsed_time=\$SECONDS\n";
			  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. creating summary statistics for the de novo assembly of $run_id...\n\" >> $log_file\n";
			  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_assembly_statistics\n";
			  if ($fq_se ne '')
				{ print SH "python $quast_path -r $ref -g $gff -t $num_threads --circos --gene-finding --rna-finding --conserved-genes-finding --pe1 $fq_1 --pe2 $fq_2 --single $fq_se -o $out_dir/$run_id/logs/$log_number"."_assembly_statistics $out_dir/$run_id/$run_id.fa &>> $log_file\n"; }
			  elsif ($fq_se eq '')
				{ print SH "python $quast_path -r $ref -g $gff -t $num_threads --circos --gene-finding --rna-finding --conserved-genes-finding --pe1 $fq_1 --pe2 $fq_2 -o $out_dir/$run_id/logs/$log_number"."_assembly_statistics $out_dir/$run_id/$run_id.fa &>> $log_file\n"; }
			}
		  if ($keep_cleaned_fq eq 'no')
			{ if 	($fq_se ne '') { print SH "rm $fq_1 $fq_2 $fq_se\n"; }
			  elsif ($fq_se eq '') { print SH "rm $fq_1 $fq_2\n"; 		 }
			}
		}
	  
	  # (M) RUN TETYPER
	  ## to be added ##
	  
	  # RUN RAXML
	  ## to be added ##
	  
	  # (N)	TIDY UP
	  # delete the now-redundant reference genome and associated indexing files
	  print SH "rm $ref $ref.amb $ref.ann $ref.bwt $ref.fai $ref.pac $ref.sa $gff\n"; # $ref.stidx $ref.sthash
	  # delete the now-redundant BAM
	  if ($keep_bam eq 'no')
		{ print SH "rm $out_dir/$run_id/$run_id.bam $out_dir/$run_id/$run_id.bai\n"; }	  
	  # make md5s of the core pipeline output: the BAM, VCF and FA
	  if ($keep_bam eq 'yes')
		{ print SH "md5sum $out_dir/$run_id/$run_id.bam > $out_dir/$run_id/$run_id.bam.md5sum\n";
		  print SH "md5sum $out_dir/$run_id/$run_id.bai > $out_dir/$run_id/$run_id.bai.md5sum\n";
		}
	  if ($keep_cleaned_fq eq 'yes')
		{ print SH "md5sum $fq_1 > $fq_1.md5sum\n";
		  print SH "md5sum $fq_2 > $fq_2.md5sum\n";
		  if ($fq_se ne '')
			{ print SH "md5sum $fq_se > $fq_se.md5sum\n"; }
		}
	  if ($create_assembly eq 'yes')
		{ print SH "md5sum $out_dir/$run_id/$run_id.fa > $out_dir/$run_id/$run_id.fa.md5sum\n";
		}
	  print SH "md5sum $out_dir/$run_id/$run_id.vcf.gz > $out_dir/$run_id/$run_id.vcf.gz.md5sum\n";
	  print SH "md5sum $out_dir/$run_id/$run_id.vcf.gz.tbi > $out_dir/$run_id/$run_id.vcf.gz.tbi.md5sum\n";
	
	  # HOW LONG DID IT TAKE TO COMPLETE THIS PIPELINE?
	  print SH "duration=\$SECONDS\n";
	  print SH "echo \"--> total time taken to process $run_id: \$((\$duration / 60)) minutes and \$((\$duration % 60)) seconds\n\" >> $log_file\n";
	  print SH "exit 0\n";
	  close(SH) or die $!;
	  close(LOG) or die $!;
	}

exit 1;