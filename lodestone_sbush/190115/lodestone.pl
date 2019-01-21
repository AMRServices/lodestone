=head
TO CREATE CENTRIFUGE DATABASE FOR HUMAN READ REMOVAL (NOTE THAT 'HUMAN' INCLUDES BOTH THE HUMA GENOME AND HIV):
centrifuge-download -P 6 -o taxonomy taxonomy
centrifuge-download -P 6 -o library -d "vertebrate_mammalian" -a "Chromosome" -t 9606 -c 'reference genome' refseq >> seqid2taxid.map # human
centrifuge-download -P 6 -o library -d "viral" -t 11676 refseq >> seqid2taxid.map # human immunodeficiency virus 1
centrifuge-download -P 6 -o library -d "viral" -t 11709 refseq >> seqid2taxid.map # human immunodeficiency virus 2
cat library/*/*.fna > input-sequences.fna
centrifuge-build -p 6 --conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp input-sequences.fna db_human_and_hiv
rm -r library
rm -r taxonomy
rm input-sequences.fna
rm seqid2taxid.map
=cut

use strict;
use warnings;
use POSIX qw(ceil strftime);

# REQUIREMENTS
my $where = 'ndm'; # One of two options: ndm or eddie. If using eddie, sh formatted for Grid Engine; if NDM, not.
my $home  = '';
if 	  ($where eq 'eddie') { $home = '/home/sbush'; 				  } # ssh sbush@eddie3.ecdf.ed.ac.uk
elsif ($where eq 'ndm')	  { $home = '/mnt/microbio/HOMES/steveb'; }
my $in_dir   = "$home/sp3";
my $progs    = "$home/programs";
my $in_file  = "$in_dir/list_of_run_ids.txt"; # a 4-column tab-separated table, with the second through fourth columns optional. Column 1: SRA run accessions. Column 2: full path to a reference genome (.fa or .fna; used by BWA for alignment). Column 3: full path to the associated GFF (used by SPAdes for assembly). Column 4: full path to the associated .faa (used by Prokka for assembly annotation).
my $test_dir = '';
if    ($where eq 'eddie') { $test_dir = "/exports/eddie/scratch/sbush/test_fqs"; }
elsif ($where eq 'ndm')   { $test_dir = "$in_dir/test_fqs"; 					 } # from make_lodestone_test_fqs.pl; this directory contains strain-specific dwgsim-simulated reads from a range of species. If $run_test_set == 1 (see 'parameter presets', below), then the contents of this directory will be used in lieu of the fastqs specified in $in_file

# OUTPUT
my $scripts_dir = "$in_dir/scripts";
if (!(-d($scripts_dir))) { mkdir $scripts_dir or die $!; }
if ($where eq 'eddie')
	{ my $qsub_commands = "$scripts_dir/qsub_commands.sh";
	  open(QSUB,'>',$qsub_commands) or die $!;
	  print QSUB "#!/bin/bash\n";
	}
my $out_dir = '';
if 	  ($where eq 'eddie') { $out_dir = '/exports/eddie/scratch/sbush/sp3_output'; }
elsif ($where eq 'ndm')   { $out_dir = "$home/sp3_output";   					  }
if (!(-d($out_dir))) { mkdir $out_dir or die $!; } # to store the output of analyses: VCFs, assemblies, and logs
my $essentials_dir = "$in_dir/essentials";
if (!(-d($essentials_dir))) 		   { mkdir $essentials_dir 			  or die $!; } # to store mash sketches of downloaded genomes (in a subdirectory, so that they do not need to be downloaded again) along with other essential content, such as the taxon ID lookup tables
if (!(-d("$essentials_dir/sketches"))) { mkdir "$essentials_dir/sketches" or die $!; }
my $taxonomy_file 				   			 = "$essentials_dir/taxonomic_relationships.tsv"; # despite its presence in the 'essentials' dir, this file is actually for reference only; its creation below is optional
my $species_id_to_genus_id_lookup  			 = "$essentials_dir/species_id_to_genus_id_lookup.tsv";
my $species_and_genus_id_to_family_id_lookup = "$essentials_dir/species_and_genus_id_to_family_id_lookup.tsv";
my $species_and_genus_id_to_order_id_lookup  = "$essentials_dir/species_and_genus_id_to_order_id_lookup.tsv";

# SOFTWARE USED
my $fatal = 0;
## paths to aligners
my $bwa_dir  = ''; if ($where eq 'eddie') { $bwa_dir  = "$progs/bwa";     } elsif ($where eq 'ndm') { $bwa_dir  = "$progs/bwa-0.7.17";     }
my $bwa_path = ''; if ($where eq 'eddie') { $bwa_path = "$progs/bwa/bwa"; } elsif ($where eq 'ndm')	{ $bwa_path = "$progs/bwa-0.7.17/bwa"; }
my $stampy_path = "$progs/stampy-1.0.32/stampy.py";
if (!(-d($bwa_dir))) 	 { $fatal++; print "ERROR: cannot find $bwa_dir\n";     }
if (!(-e($bwa_path))) 	 { $fatal++; print "ERROR: cannot find $bwa_path\n";    }
if (!(-e($stampy_path))) { $fatal++; print "ERROR: cannot find $stampy_path\n"; }
## paths to assemblers and associated software, such as annotators
my $quast_path    = ''; if ($where eq 'eddie') { $quast_path = '/exports/eddie/scratch/sbush/quast-5.0.1/quast.py'; } elsif ($where eq 'ndm') { $quast_path = "$progs/quast-5.0.1/quast.py"; }
my $circos_dir    = "$progs/circos-0.69-6/bin";
my $circos_path   = "$progs/circos-0.69-6/bin/circos";
my $spades_path   = "$progs/SPAdes-3.13.0-Linux/bin/spades.py";
my $prokka_path	  = "$progs/prokka/bin/prokka";
my $barrnap_dir	  = "$progs/barrnap/bin";
my $barrnap_path  = "$progs/barrnap/bin/barrnap";
if (!(-e($quast_path)))    { $fatal++; print "ERROR: cannot find $quast_path\n";    }
if (!(-d($circos_dir)))    { $fatal++; print "ERROR: cannot find $circos_dir\n";    }
if (!(-e($circos_path)))   { $fatal++; print "ERROR: cannot find $circos_path\n";   }
if (!(-e($spades_path)))   { $fatal++; print "ERROR: cannot find $spades_path\n";   }
if (!(-e($prokka_path)))   { $fatal++; print "ERROR: cannot find $prokka_path\n";   }
if (!(-d($barrnap_dir)))   { $fatal++; print "ERROR: cannot find $barrnap_dir\n";   }
if (!(-e($barrnap_path)))  { $fatal++; print "ERROR: cannot find $barrnap_path\n";  }
## paths to callers
my $bcftools_dir  = ''; if ($where eq 'eddie') { $bcftools_dir  = "$progs/bcftools-1.7"; 		  } elsif ($where eq 'ndm') { $bcftools_dir  = "$progs/bcftools"; 		   }
my $bcftools_path = ''; if ($where eq 'eddie') { $bcftools_path = "$progs/bcftools-1.7/bcftools"; } elsif ($where eq 'ndm') { $bcftools_path = "$progs/bcftools/bcftools"; }
if (!(-d($bcftools_dir)))  { $fatal++; print "ERROR: cannot find $bcftools_dir\n";  }
if (!(-e($bcftools_path))) { $fatal++; print "ERROR: cannot find $bcftools_path\n"; }
## paths to classifiers and their associated databases
my $atlas_path 		  = ''; if ($where eq 'eddie') { $atlas_path = '/exports/eddie3_homes_local/sbush/.local/bin/atlas'; } elsif ($where eq 'ndm') { $atlas_path = '/home/ndm.local/steveb/.local/bin/atlas'; }
my $mykrobe_path 	  = ''; if ($where eq 'eddie') { $mykrobe_path = '/exports/eddie3_homes_local/sbush/.local/bin/mykrobe'; } elsif ($where eq 'ndm') { $mykrobe_path = '/home/ndm.local/steveb/.local/bin/mykrobe'; }
my $mccortex_path 	  = ''; if ($where eq 'eddie') { $mccortex_path = '/exports/eddie/scratch/sbush/mccortex/bin/mccortex31'; } elsif ($where eq 'ndm')	{ $mccortex_path = "$progs/mccortex/bin/mccortex31"; }
my $kaiju_db_p 		  = ''; if ($where eq 'eddie') { $kaiju_db_p = "/exports/eddie/scratch/sbush/db.p/kaiju_db.fmi"; } elsif ($where eq 'ndm') { $kaiju_db_p = "$progs/kaiju/db.p/kaiju_db.fmi"; } # see http://kaiju.binf.ku.dk/server; specifically, http://kaiju.binf.ku.dk/database/kaiju_index_pg.tgz (the proGenomes database, dated 16th May 2017)
my $kaiju_db_e   	  = ''; if ($where eq 'eddie') { $kaiju_db_e = "/exports/eddie/scratch/sbush/db.e/kaiju_db_nr_euk.fmi"; } elsif ($where eq 'ndm') { $kaiju_db_e = "$progs/kaiju/db.e/kaiju_db_nr_euk.fmi"; } # see http://kaiju.binf.ku.dk/server; specifically, http://kaiju.binf.ku.dk/database/kaiju_index_nr_euk.tgz (the expanded nr database, containing bacteria, archaea, viruses, fungi and microbial eukaryotes, dated 16th May 2017)
my $nodes_dmp  		  = ''; if ($where eq 'eddie') { $nodes_dmp = "/exports/eddie/scratch/sbush/nodes.dmp"; } elsif ($where eq 'ndm') { $nodes_dmp = "$progs/kaiju/nodes.dmp"; } # if not using the nodes.dmp packaged with the above pre-computed databases, then obtain using: wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | tar xvzf taxdump.tar.gz nodes.dmp
my $names_dmp         = ''; if ($where eq 'eddie') { $names_dmp = "/exports/eddie/scratch/sbush/names.dmp"; } elsif ($where eq 'ndm') { $names_dmp = "$progs/kaiju/names.dmp"; } # if not using the names.dmp packaged with the above pre-computed databases, then obtain using: wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | tar xvzf taxdump.tar.gz names.dmp
my $kaiju_path 		  = "$progs/kaiju/bin/kaiju";
my $kaijuReport_path  = "$progs/kaiju/bin/kaijuReport";
my $kaijutokrona_path = "$progs/kaiju/bin/kaiju2krona";
my $mergeOutputs_path = "$progs/kaiju/bin/mergeOutputs";
my $centrifuge_path	  = "$progs/centrifuge/centrifuge";
my $centrifuge_db     = ''; if ($where eq 'eddie') { $centrifuge_db = "/exports/eddie/scratch/sbush/db.human_and_hiv"; } elsif ($where eq 'ndm') { $centrifuge_db = "$progs/centrifuge/db.human_and_hiv"; } # manually created
my $centrifuge_index_prefix = "$centrifuge_db/db_human_and_hiv"; # a Centrifuge index comprises files of the format prefix.X.cf
if (!(-e($atlas_path)))        { $fatal++; print "ERROR: cannot find $atlas_path\n";   		}
if (!(-e($mykrobe_path)))      { $fatal++; print "ERROR: cannot find $mykrobe_path\n";   	}
if (!(-e($mccortex_path))) 	   { $fatal++; print "ERROR: cannot find $mccortex_path\n"; 	}
if (!(-e($kaiju_db_p)))   	   { $fatal++; print "ERROR: cannot find $kaiju_db_p\n";   	    }
if (!(-e($kaiju_db_e)))   	   { $fatal++; print "ERROR: cannot find $kaiju_db_e\n";   	    }
if (!(-e($nodes_dmp)))  	   { $fatal++; print "ERROR: cannot find $nodes_dmp\n";  	    }
if (!(-e($names_dmp))) 	   	   { $fatal++; print "ERROR: cannot find $names_dmp\n";  	    }
if (!(-e($kaiju_path))) 	   { $fatal++; print "ERROR: cannot find $kaiju_path\n"; 	    }
if (!(-e($kaijuReport_path)))  { $fatal++; print "ERROR: cannot find $kaijuReport_path\n";  }
if (!(-e($kaijutokrona_path))) { $fatal++; print "ERROR: cannot find $kaijutokrona_path\n"; }
if (!(-e($mergeOutputs_path))) { $fatal++; print "ERROR: cannot find $mergeOutputs_path\n"; }
if (!(-e($centrifuge_path)))   { $fatal++; print "ERROR: cannot find $centrifuge_path\n";   }
if (!(-d($centrifuge_db)))     { $fatal++; print "ERROR: cannot find $centrifuge_db\n";     }
## paths to ancillary software
my $ascp_path				  = ''; if ($where eq 'eddie') { $ascp_path = "/home/sbush/.aspera/connect/bin/ascp"; } elsif ($where eq 'ndm') { $ascp_path = "/home/ndm.local/steveb/.aspera/connect/bin/ascp"; }
my $aspera_openssh_path		  = ''; if ($where eq 'eddie') { $aspera_openssh_path = "/home/sbush/.aspera/connect/etc/asperaweb_id_dsa.openssh"; } elsif ($where eq 'ndm') { $aspera_openssh_path = "/home/ndm.local/steveb/.aspera/connect/etc/asperaweb_id_dsa.openssh"; }
my $mash_path				  = "$progs/mash-Linux64-v2.1/mash";
my $bgzip_path   			  = "$progs/bin/bgzip";
my $tabix_path   			  = "$progs/bin/tabix"; # tabix and bgzip are components of the SAMtools package, commonly used to compress and index VCFs; see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/
my $seqtk_path 				  = "$progs/seqtk/seqtk";
my $fastqc_path				  = "$progs/FastQC/fastqc";
my $picard_path   			  = "$progs/picard.jar";
my $efetch_path 			  = "$progs/edirect/efetch"; # for a fuller description of the Entrez Direct tools (used to automate downloads from NCBI), see https://www.ncbi.nlm.nih.gov/books/NBK179288/
my $xtract_path 			  = "$progs/edirect/xtract";
my $esearch_path 			  = "$progs/edirect/esearch";
my $esummary_path 			  = "$progs/edirect/esummary";
my $blast_dir 				  = ''; if ($where eq 'eddie') { $blast_dir = "/exports/eddie/scratch/sbush/ncbi-blast-2.2.25+/bin"; } elsif ($where eq 'ndm') { $blast_dir = "$progs/ncbi-blast-2.2.25+/bin"; }
my $blastn_path				  = ''; if ($where eq 'eddie') { $blastn_path = "/exports/eddie/scratch/sbush/ncbi-blast-2.2.25+/bin/blastn"; } elsif ($where eq 'ndm') { $blastn_path = "$progs/ncbi-blast-2.2.25+/bin/blastn"; } # IMPORTANT: required for TETyper. Note that this expects blastn on the $PATH but that more recent versions of the BLAST+ toolkit have deprecated this command.
my $tetyper_path			  = ''; if ($where eq 'eddie') { $tetyper_path = "/home/sbush/.local/bin/TETyper.py"; } elsif ($where eq 'ndm') { $tetyper_path = "/home/ndm.local/steveb/.local/bin/TETyper.py"; }
my $tetyper_ref 			  = "$progs/TETyper/Tn4401b-1.fasta"; # required for TETyper, but not automatically obtained with pip installation. Manually downloaded from https://github.com/aesheppard/TETyper/raw/master/Tn4401b-1.fasta
my $tetyper_snp_profiles      = "$progs/TETyper/snp_profiles.txt"; # required for TETyper, but not automatically obtained with pip installation. Manually downloaded from https://github.com/aesheppard/TETyper/raw/master/snp_profiles.txt
my $tetyper_struct_profiles   = "$progs/TETyper/struct_profiles.txt"; # required for TETyper, but not automatically obtained with pip installation. Manually downloaded from https://github.com/aesheppard/TETyper/raw/master/struct_profiles.txt
my $vcfsort_path 			  = "$progs/bin/vcf-sort"; # a component of VCFtools; see https://vcftools.github.io/index.html
my $samtools_dir			  = "$progs/samtools-1.7";
my $samtools_path 			  = "$progs/samtools-1.7/samtools";
my $sambamba_path			  = "$progs/sambamba";
my $cutadapt_path			  = ''; if ($where eq 'eddie'){ $cutadapt_path = "/home/sbush/.local/bin/cutadapt"; } elsif ($where eq 'ndm') { $cutadapt_path = "/home/ndm.local/steveb/.local/bin/cutadapt"; }
my $trim_galore_path		  = "$progs/TrimGalore-0.5.0/trim_galore";
my $trimmomatic_path 		  = "$progs/Trimmomatic-0.38/trimmomatic-0.38.jar";
my $ktImportText_path		  = ''; if ($where eq 'eddie') { $ktImportText_path = "$progs/bin/ktImportText"; } elsif ($where eq 'ndm') { $ktImportText_path = "$progs/KronaTools-2.7/bin/ktImportText"; }
my $vcffilter_path 		  	  = "$progs/vcflib/bin/vcffilter";
my $vcfoverlay_path 		  = "$progs/vcflib/bin/vcfoverlay"; # not actually used in the pipeline itself and only listed here to test the version number of VCFlib (for the --version command is not applicable to all VCFlib modules)
my $vcfallelicprimitives_path = "$progs/vcflib/bin/vcfallelicprimitives";
if (!(-e($mash_path)))    			   { $fatal++; print "ERROR: cannot find $mash_path\n";    				}
if (!(-e($ascp_path)))    			   { $fatal++; print "ERROR: cannot find $ascp_path\n";    				}
if (!(-d($blast_dir)))    			   { $fatal++; print "ERROR: cannot find $blast_dir\n";    				}
if (!(-e($bgzip_path)))   			   { $fatal++; print "ERROR: cannot find $bgzip_path\n";   				}
if (!(-e($tabix_path)))   			   { $fatal++; print "ERROR: cannot find $tabix_path\n";   				}
if (!(-e($seqtk_path)))   			   { $fatal++; print "ERROR: cannot find $seqtk_path\n";   				}
if (!(-e($fastqc_path)))   			   { $fatal++; print "ERROR: cannot find $fastqc_path\n";   			}
if (!(-e($picard_path)))   			   { $fatal++; print "ERROR: cannot find $picard_path\n";   			}
if (!(-e($efetch_path)))   		   	   { $fatal++; print "ERROR: cannot find $efetch_path\n";   			}
if (!(-e($xtract_path)))   		   	   { $fatal++; print "ERROR: cannot find $xtract_path\n";   			}
if (!(-e($esearch_path)))   		   { $fatal++; print "ERROR: cannot find $esearch_path\n";   			}
if (!(-e($esummary_path)))   		   { $fatal++; print "ERROR: cannot find $esummary_path\n";   			}
if (!(-e($blastn_path)))	  		   { $fatal++; print "ERROR: cannot find $blastn_path\n";   			}
if (!(-e($tetyper_path)))   		   { $fatal++; print "ERROR: cannot find $tetyper_path\n";   			}
if (!(-e($tetyper_ref)))   		   	   { $fatal++; print "ERROR: cannot find $tetyper_ref\n";   			}
if (!(-e($tetyper_snp_profiles)))      { $fatal++; print "ERROR: cannot find $tetyper_snp_profiles\n";   	}
if (!(-e($tetyper_struct_profiles)))   { $fatal++; print "ERROR: cannot find $tetyper_struct_profiles\n";   }
if (!(-e($vcfsort_path)))   		   { $fatal++; print "ERROR: cannot find $vcfsort_path\n";   			}
if (!(-d($samtools_dir))) 			   { $fatal++; print "ERROR: cannot find $samtools_dir\n"; 				}
if (!(-e($samtools_path))) 			   { $fatal++; print "ERROR: cannot find $samtools_path\n"; 			}
if (!(-e($sambamba_path))) 			   { $fatal++; print "ERROR: cannot find $sambamba_path\n"; 			}
if (!(-e($cutadapt_path))) 			   { $fatal++; print "ERROR: cannot find $cutadapt_path\n"; 			}
if (!(-e($trim_galore_path)))		   { $fatal++; print "ERROR: cannot find $trim_galore_path\n";			}
if (!(-e($trimmomatic_path)))		   { $fatal++; print "ERROR: cannot find $trimmomatic_path\n";			}
if (!(-e($ktImportText_path)))		   { $fatal++; print "ERROR: cannot find $ktImportText_path\n";			}
if (!(-e($aspera_openssh_path)))	   { $fatal++; print "ERROR: cannot find $aspera_openssh_path\n";		}
if (!(-e($vcffilter_path))) 		   { $fatal++; print "ERROR: cannot find $vcffilter_path\n"; 			}
if (!(-e($vcfoverlay_path))) 		   { $fatal++; print "ERROR: cannot find $vcfoverlay_path\n"; 			}
if (!(-e($vcfallelicprimitives_path))) { $fatal++; print "ERROR: cannot find $vcfallelicprimitives_path\n"; }
exit if ($fatal > 0);

# PARAMETER PRESETS
# These presets will override any options manually set below. You can only pick ONE.
my $run_test_set  = 0;
my $quick_varcall = 0;
my $clean_run_of_everything = 1;

# PARAMETERS
my $num_threads = 6;
my $phred_threshold = 20;
my $minimum_read_length = 100;
my $maximum_pc_of_n_bases = ceil(0.05*$minimum_read_length); # i.e. 5% of the bases in a read of a given minimum length are allowed to be N
my $num_of_bases_overlapping_with_adapter_before_trimming = 5; # TrimGalore parameter -s; note the default is the very stringent value of 1
my $min_pc_of_reads_for_being_confident_in_species_classification = 50; # if lower than this, we don't take the top hit species classification; we take the genus
my $min_pc_of_reads_for_being_confident_in_genus_classification   = 50; # if lower than this, we don't take the top hit genus classification; we take the family
my $min_pc_of_reads_for_being_confident_in_family_classification  = 50; # if lower than this, we don't take the top hit family classification; we take the order
my $min_pc_of_reads_for_being_confident_in_order_classification   = 50; # if lower that this, we abort - unless $fail_beyond_order eq 'no' - and report that we do not have confidence in a pure, single-species, input
my $tetyper_flank_length = 20; # essential TETyper parameter --flank_len
my $run_initial_fastqc = 'yes';
my $run_trimgalore     = 'yes';
my $run_trimmomatic    = 'yes';
my $run_final_fastqc   = 'yes';
my $use_the_provided_ref_genome = 'yes'; # if doing so, we will not then need to run Kaiju & Mash
my $run_kaiju_on_se_reads = 'no'; # if 'yes', this allows for an internal test of the consistency of classification, as the same top hit species/genus/family/order should be called by both the PE and SE reads (the latter inevitably a smaller, noiser, set compared to the former). This requires that $run_trimmomatic eq 'yes', as de facto SE reads are only created by this step.
my $run_kaiju_on_db_p = 'yes';
my $run_kaiju_on_db_e = 'yes';
my $keep_read_classifications = 'no'; # if 'yes', we retain the read classification files of Kaiju and Centrifuge (they are otherwise considered interim files and deleted after condensed summary files are made); note these files are large as they are lists of each read ID and associated prediction of each taxon ID
my $fail_beyond_order  = 'no'; # if 'yes', we will fail if too low a % of reads ($min_pc_of_reads_for_being_confident_in_order_classification %) can be assigned to a certain taxonomic rank (in this case, the comparatively broad rank of order)
my $run_mykrobe_if_tb  = 'no'; # if 'yes', run "Mykrobe predict" to predict strain from TB reads
my $remove_supplementary_alignments = 'yes'; # see SAM file specification point 1.2 ('chimeric alignment'): https://samtools.github.io/hts-specs/SAMv1.pdf
my $remove_non_primary_alignments   = 'yes'; # see SAM file specification point 1.2 ('multiple mapping'): https://samtools.github.io/hts-specs/SAMv1.pdf
my $output_all_sites   = 'no'; # for variant calling with mpileup; 'yes' to output all sites, 'no' to output only the variant sites
my $regularise_vcf     = 'yes';
my $filter_vcf		   = 'no';
my $keep_original_vcf  = 'no';
my $max_number_of_alleles_at_site = 1; # for filtering the VCF on INFO field 'AN' (i.e. removing variants that fail this filter): total number of alleles in called genotypes. Set to arbitrarily high if you want to capture everything. A value of 1 captures bimodal variants only.
my $min_MQ  = 10; # for filtering the VCF on VCF INFO field 'MQ' (i.e. removing variants that fail this filter): average mapping quality (Phred)
my $min_DP4 = 10; # for filtering the VCF on VCF INFO field 'DP4' (i.e. removing variants that fail this filter): total number of high-quality ref-forward, ref-reverse, alt-forward and alt-reverse bases
my $create_assembly    = 'no'; # if 'yes', run SPAdes on the QC'd, cleaned and classified raw reads; if 'no', don't
my $qc_assembly        = 'no'; # requires that $create_assembly also eq 'yes'
my $annotate_assembly  = 'no'; # if 'yes', run Prokka on the SPAdes assembly
my $run_tetyper		   = 'no';
my $keep_bam		   = 'no';
my $keep_cleaned_fq	   = 'yes';
if (($run_tetyper 		eq 'yes') and ($create_assembly eq 'no'))  { $create_assembly = 'yes'; }
if (($annotate_assembly eq 'yes') and ($create_assembly eq 'no'))  { $create_assembly = 'no';  }
if (($create_assembly 	eq 'no')  and ($qc_assembly     eq 'yes')) { $qc_assembly 	  = 'no';  }
## apply parameter pre-sets, if so requested (these override any options manually chosen above)
if ($quick_varcall == 1) # runs the lowest-overhead Kaiju (database P) only on the PE reads; does not create assemblies or FastQC reports and does not store intermediary files (e.g. cleaned fqs, BAM, pre-filtered VCF)
	{ $run_initial_fastqc	 	   = 'no';
	  $run_final_fastqc		 	   = 'no';
	  $use_the_provided_ref_genome = 'yes';
	  $run_kaiju_on_se_reads 	   = 'no';
	  $run_kaiju_on_db_p 	 	   = 'yes';
	  $run_kaiju_on_db_e 	 	   = 'no';
	  $keep_read_classifications   = 'yes';
	  $fail_beyond_order     	   = 'yes';
	  $regularise_vcf			   = 'yes';
	  $filter_vcf				   = 'no';
	  $keep_original_vcf 	 	   = 'no';
	  $keep_bam		   	     	   = 'yes';
	  $keep_cleaned_fq   	 	   = 'yes';
	  $run_tetyper		     	   = 'no';
	  $create_assembly   	 	   = 'no';
	  $qc_assembly				   = 'no';
	  $annotate_assembly 	 	   = 'no';
	}
if ($clean_run_of_everything == 1)
	{ $run_initial_fastqc	 	   	   = 'yes';
	  $run_final_fastqc		 	   	   = 'yes';
	  $use_the_provided_ref_genome 	   = 'no';
	  $run_kaiju_on_se_reads 	   	   = 'yes';
	  $run_kaiju_on_db_p 	 	   	   = 'yes';
	  $run_kaiju_on_db_e 	 	   	   = 'yes';
	  $keep_read_classifications   	   = 'no';
	  $fail_beyond_order     	   	   = 'yes';
	  $remove_supplementary_alignments = 'yes';
	  $remove_non_primary_alignments   = 'yes';
	  $regularise_vcf			   	   = 'yes';
	  $filter_vcf				   	   = 'no';
	  $keep_original_vcf 	 	   	   = 'no';
	  $keep_bam		   	     	   	   = 'no';
	  $keep_cleaned_fq   	 	   	   = 'no';
	  $run_tetyper		     	   	   = 'no';
	  $create_assembly   	 	   	   = 'yes';
	  $qc_assembly				   	   = 'yes';
	  $annotate_assembly 	 	   	   = 'yes';
	}

##########

# BEFORE RUNNING THE PIPELINE, PARSE NAMES.DMP AND NODES.DMP TO OBTAIN THE TAXONOMIC RELATIONSHIPS ASSOCIATED WITH EACH TAXON ID. THIS WILL ASSOCIATE EACH ORDER, FAMILY AND GENUS WITH A LIST OF SPECIES IDS. WE NEED THIS INFORMATION IN ORDER TO LATER PARSE ASSEMBLY_SUMMARY.TXT TO DOWNLOAD ALL GENOMES OF A PARTICULAR TAXONOMIC RANK. GENOMES IN THIS FILE ARE LISTED ONLY BY *SPECIES* TAXON ID. WE NEED TO LOOK UP THE BROADER TAXONOMIC RANKS THAT THEY ARE ALSO MEMBERS OF, IN ORDER TO THEN LOOK UP THE LIST OF *ALL* GENOMES ALSO AT THAT RANK.
if ( (!(-e($species_id_to_genus_id_lookup))) || (!(-e($species_and_genus_id_to_family_id_lookup))) || (!(-e($species_and_genus_id_to_order_id_lookup))) )
	{ print "creating taxonomic relationship file; please wait...\n";
	  open(TAXON,'>',$taxonomy_file) or die $!; open(GENUS,'>',$species_id_to_genus_id_lookup) or die $!; open(FAMILY,'>',$species_and_genus_id_to_family_id_lookup) or die $!; open(ORDER,'>',$species_and_genus_id_to_order_id_lookup) or die $!;
	  print TAXON "Taxon ID (1)\tTaxon name (1)\tTaxon rank (1)\tTaxon ID (2)\tTaxon name (2)\tTaxon rank (2)\tTaxon ID (3)\tTaxon name (3)\tTaxon rank (3)\tTaxon ID (4)\tTaxon name (4)\tTaxon rank (4)\tTaxon ID (5)\tTaxon name (5)\tTaxon rank (5)\tTaxon ID (6)\tTaxon name (6)\tTaxon rank (6)\tTaxon ID (7)\tTaxon name (7)\tTaxon rank (7)\tTaxon ID (8)\tTaxon name (8)\tTaxon rank (8)\tTaxon ID (9)\tTaxon name (9)\tTaxon rank (9)\tTaxon ID (10)\tTaxon name (10)\tTaxon rank (10)\n";
	  my %taxon_names = ();
	  open(IN,$names_dmp) or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line); $line =~ s/^\s+//;
		  $line =~ s/\s/\t/;
		  my @line = split(/\t/,$line);
		  my $taxon_id = $line[0]; my $taxon_name = $line[2]; my $entry_type = $line[6];
		  if ($entry_type eq 'scientific name')
			{ $taxon_names{$taxon_id} = $taxon_name;
			}
		}
	  close(IN) or die $!;
	  my %ranks_per_taxon_id = (); my %parents_per_taxon_id = ();
	  open(IN,$nodes_dmp) or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line); $line =~ s/^\s+//;
		  $line =~ s/\s/\t/;
		  my @line = split(/\t/,$line);
		  my $taxon_id = $line[0]; my $parent_taxon_id = $line[2]; my $rank = $line[4];
		  $ranks_per_taxon_id{$taxon_id} = $rank;
		  $parents_per_taxon_id{$taxon_id} = $parent_taxon_id;
		}
	  close(IN) or die $!;
	  my @taxon_ids = ();
	  while((my $taxon_id,my $irrel)=each(%taxon_names))
		{ push(@taxon_ids,$taxon_id); }
	  my @sorted_taxon_ids = sort {$a <=> $b} @taxon_ids;
	  my %ranks_for_each_taxon_id = ();
	  my %species_id_to_genus_id_lookup = (); my %species_and_genus_id_to_family_id_lookup = (); my %species_and_genus_id_to_order_id_lookup = ();
	  foreach my $taxon_id (@sorted_taxon_ids)
		{ my $taxon_name = $taxon_names{$taxon_id};
		  my $taxon_rank = $ranks_per_taxon_id{$taxon_id};
		  $ranks_for_each_taxon_id{$taxon_id}{$taxon_rank}{$taxon_name}++;
		  my $out_line = '';
		  my $last_input_id = '';
		  for(my $x=1;$x<=10;$x++)
			{ my $parental_taxon_id = 'NA'; my $parental_taxon_name = 'NA'; my $parental_taxon_rank = 'NA';
			  my $input_id = '';
			  if ($x == 1)
				{ $input_id = $taxon_id; }
			  else
				{ $input_id = $last_input_id unless ($last_input_id =~ //); }
			  if (exists($parents_per_taxon_id{$input_id}))
				{ $parental_taxon_id   = $parents_per_taxon_id{$input_id};
				  $parental_taxon_name = $taxon_names{$parental_taxon_id};
				  $parental_taxon_rank = $ranks_per_taxon_id{$parental_taxon_id};
				  $last_input_id = $parental_taxon_id;
				}
			  else
				{ $last_input_id = ''; }
			  $out_line .= "$parental_taxon_id\t$parental_taxon_name\t$parental_taxon_rank\t";
			}
		  $out_line =~ s/\t$//;
		  print TAXON "$taxon_id\t$taxon_name\t$taxon_rank\t$out_line\n";
		  my @out_line = split(/\t/,$out_line);
		  if ($taxon_rank eq 'species')
			{ for(my $x=0;$x<@out_line;$x++)
				{ my $parental_taxon_id = $out_line[$x]; my $parental_taxon_name = $out_line[$x+1]; my $parental_taxon_rank = $out_line[$x+2];
				  if ($parental_taxon_rank eq 'genus')	{ $species_id_to_genus_id_lookup{$taxon_id}{$parental_taxon_id}++;  		  }
				  if ($parental_taxon_rank eq 'family') { $species_and_genus_id_to_family_id_lookup{$taxon_id}{$parental_taxon_id}++; }
				  if ($parental_taxon_rank eq 'order')  { $species_and_genus_id_to_order_id_lookup{$taxon_id}{$parental_taxon_id}++;  }
				  $x = $x+2;
				}
			}
		  if ($taxon_rank eq 'genus')
			{ for(my $x=0;$x<@out_line;$x++)
				{ my $parental_taxon_id = $out_line[$x]; my $parental_taxon_name = $out_line[$x+1]; my $parental_taxon_rank = $out_line[$x+2];
				  if ($parental_taxon_rank eq 'family') { $species_and_genus_id_to_family_id_lookup{$taxon_id}{$parental_taxon_id}++; }
				  if ($parental_taxon_rank eq 'order')  { $species_and_genus_id_to_order_id_lookup{$taxon_id}{$parental_taxon_id}++;  }
				  $x = $x+2;
				}
			}
		}
	  close(TAXON) or die $!;
	  foreach my $taxon_id (sort {$a <=> $b} keys %species_id_to_genus_id_lookup)
		{ foreach my $genus_taxon_id (sort {$a <=> $b} keys %{$species_id_to_genus_id_lookup{$taxon_id}})
			{ print GENUS "$taxon_id\t$genus_taxon_id\n";
			}
		}
	  close(GENUS) or die $!;
	  foreach my $taxon_id (sort {$a <=> $b} keys %species_and_genus_id_to_family_id_lookup)
		{ foreach my $family_taxon_id (sort {$a <=> $b} keys %{$species_and_genus_id_to_family_id_lookup{$taxon_id}})
			{ print FAMILY "$taxon_id\t$family_taxon_id\n";
			}
		}
	  close(FAMILY) or die $!;
	  foreach my $taxon_id (sort {$a <=> $b} keys %species_and_genus_id_to_order_id_lookup)
		{ foreach my $order_taxon_id (sort {$a <=> $b} keys %{$species_and_genus_id_to_order_id_lookup{$taxon_id}})
			{ print ORDER "$taxon_id\t$order_taxon_id\n";
			}
		}
	  close(ORDER) or die $!;
	}

# OBTAIN LIST OF SRA RUN ACCESSIONS (FOR RETRIEVAL VIA THE ENA), OR THE LOCATIONS OF SIMULATED FQS FOR TESTING PURPOSES. THE FULL PATH TO A REFERENCE GENOME AND GFF MAY ALSO BE SPECIFIED IN THE SECOND AND THIRD COLUMNS, RESPECTIVELY. IF SO, ALIGNMENTS CAN BE MADE AGAINST THE GIVEN REFERENCE GENOME DIRECTLY, INSTEAD OF HAVING TO PREDICT IT USING KAIJU/MASH. NOTE THAT THE GFF IS NOT NEEDED FOR ALIGNMENT AND VARIANT CALLING BUT ONLY FOR QUAST, I.E. IF YOU ARE QC'ING A SPADES ASSEMBLY.
my %sra_run_ids = ();
if ($run_test_set == 0)
	{ if (!(-e($in_file)))
		{ print "ERROR: cannot open $in_file\n";
		  exit 1;
		}
	  open(IN,$in_file) or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line); $line =~ s/\r//g;
		  my @line = split(/\t/,$line);
		  my $run_id = $line[0]; my $ref_genome = ''; my $ref_gff = ''; my $ref_faa = '';
		  if (defined($line[1])) { $ref_genome = $line[1]; }
		  if (defined($line[2])) { $ref_gff    = $line[2]; }
		  if (defined($line[3])) { $ref_faa    = $line[3]; }
		  $sra_run_ids{$run_id}{fa}  = $ref_genome;
		  $sra_run_ids{$run_id}{gff} = $ref_gff;
		  $sra_run_ids{$run_id}{faa} = $ref_faa;
		}
	  close(IN) or die $!;
	}
elsif ($run_test_set == 1)
	{ if (-d($test_dir))
		{ opendir(DIR,$test_dir) or die $!;
		  my @files = readdir(DIR);
		  closedir(DIR) or die $!;
		  my %simulated_reads = ();
		  foreach my $file (@files)
			{ if ($file =~ /^(.*?)\.(.*?)\.(\d{1})\.read(\d{1})\.fq\.gz$/)
				{ my $species = $1; my $strain = $2; my $replicate_num = $3; my $read_number = $4;
				  $simulated_reads{$species}{$strain}{$replicate_num}{$read_number} = $file;
				}
			}
		  while((my $species,my $irrel)=each(%simulated_reads))
			{ while((my $strain,my $irrel)=each(%{$simulated_reads{$species}}))
				{ while((my $replicate_num,my $irrel)=each(%{$simulated_reads{$species}{$strain}}))
					{ if ( (exists($simulated_reads{$species}{$strain}{$replicate_num}{'1'})) and (exists($simulated_reads{$species}{$strain}{$replicate_num}{'2'})) )
						{ my $read1_path = $simulated_reads{$species}{$strain}{$replicate_num}{'1'};
						  my $read2_path = $simulated_reads{$species}{$strain}{$replicate_num}{'2'};
						  $sra_run_ids{"$species.$strain.$replicate_num"}{read1_path} = "$test_dir/$read1_path";
						  $sra_run_ids{"$species.$strain.$replicate_num"}{read2_path} = "$test_dir/$read2_path";
						}
					}
				}
			}
		}
	}

# FOR EACH RUN ID, CREATE A FILE OF COMMANDS THAT IMPLEMENTS A COMPLETE QC-SPECIATE-ALIGN-VARCALL-ASSEMBLY PIPELINE
my $run_ids_seen = 0; my $run_ids_total = scalar keys %sra_run_ids;
foreach my $run_id (sort keys %sra_run_ids)
	{ $run_ids_seen++;
	  #next if (($run_test_set == 0) and ($run_id ne 'SRR2667440')); # FOR TESTING PURPOSES: reads are from Mycobacterium bovis. Will fail as Kaiju cannot classify this even at the level of order - unless we overrule Kaiju with Mykrobe
	  #next if (($run_test_set == 0) and ($run_id ne 'DRR076046')); # FOR TESTING PURPOSES: the reads for this ID are single-end so this should trigger an unable-to-download failure
	  #next if (($run_test_set == 0) and ($run_id ne 'SRR2965690')); # FOR TESTING PURPOSES: reads are from Citrobacter freundii CAV1321, Kaiju will download based on genus (if using only database P, and no SE reads), and the closest reference assembly is expected to be GCF_001022155.1_ASM102215v1
	  #next if (($run_test_set == 0) and ($run_id ne 'SRR2965739')); # FOR TESTING PURPOSES: reads are from Citrobacter freundii CAV1741, Kaiju will download based on genus (if using only database P, and no SE reads), and the closest reference assembly is expected to be GCF_001022275.1_ASM102227v1
	  next if (($run_test_set == 0) and ($run_id ne 'SRR974838')); # FOR TESTING PURPOSES: reads are from Mycobacterium tuberculosis F11
	  #next if (($run_test_set == 0) and ($run_id ne 'SRR2965590')); # FOR TESTING PURPOSES: reads are from Citrobacter amalonaticus, Kaiju will download based on family (if using only database P, and no SE reads)
	  next if (($run_test_set == 0) and ($run_id !~ /^\w+\d+$/)); # CHECKPOINT: skip unrecognisable SRA run IDs
	  
	  my $log_file = "$out_dir/$run_id.log"; # OPTIONAL: we can time-stamp log files by appending an appropriate variable, e.g. my $now_gmt  = strftime "%e"."_"."%b"."_"."%Y"."_"."%H"."h_"."%M"."m_"."%S"."s", gmtime; # e.g 20_Oct_2018_08h_23m_58s
	  next if (-d("$out_dir/$run_id")); # CHECKPOINT: we have processed data from this run ID already
	  
	  open(LOG,'>',$log_file) or die $!;
	  
	  my $sh_file = "$scripts_dir/commands.$run_id.sh";
	  open(SH,'>',$sh_file) or die $!;
	  print SH "#!/bin/bash\n";
	  if ($where eq 'eddie')
		{ print SH '#$ -N '; 			  print SH "$run_id\n";
		  print SH '#$ -cwd';			  print SH "\n";
		  print SH '#$ -pe sharedmem 6';  print SH "\n";
		  print SH '#$ -l h_rt=05:00:00'; print SH "\n";
		  print SH '#$ -l h_vmem=20G';    print SH "\n";
		  print SH ". /etc/profile.d/modules.sh\n";
		  print SH "module load java/jdk/1.8.0\n";
		  print SH "module load roslin/R/3.4.3\n"; # required for Quast
		  print SH "module load roslin/python/2.7.13\n";
		  print SH "export C_INCLUDE_PATH=/home/sbush/programs/include\n"; # required for mpileup
		  print SH "export LD_LIBRARY_PATH=/home/sbush/programs/lib\n"; # required for mpileup
		  print SH "PATH=\$PATH:/home/sbush/programs/bin\n"; # required for samtools
		}
	  print SH "PATH=\$PATH:$circos_dir\n"; # required for Quast
	  print SH "PATH=\$PATH:$barrnap_dir\n"; # required for Prokka
	  print SH "PATH=\$PATH:$samtools_dir\n"; # required for TETyper
	  print SH "PATH=\$PATH:$bcftools_dir\n"; # required for TETyper
	  print SH "PATH=\$PATH:$blast_dir\n"; # required for TETyper
	  print SH "PATH=\$PATH:$bwa_dir\n"; # required for TETyper
	  if ($where eq 'eddie')
		{ print SH "export PYTHONPATH=\$PYTHONPATH:/usr/lib64/python2.7/site-packages\n"; # required for Quast (to find the location of matplotlib)
		  print SH "export PERL5LIB=/home/sbush/perl5/lib/perl5\n"; # required for Circos (it's the location of essential Perl modules; see https://stackoverflow.com/questions/1557959/how-can-i-find-out-where-a-perl-module-is-installed)
		}
	  print SH "SECONDS=0\n"; # start the clock
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "echo \"-->\$current_date_time. beginning analysis of $run_id\n\" >> $log_file\n";
	  
 	  # (1) RECORD THE PARAMETERS USED FOR THE PIPELINE
	  print LOG "### Parameters ###\n";
	  print LOG "minimum Phred score for quality trimming: $phred_threshold\n";
	  print LOG "minimum read length for quality trimming: $minimum_read_length\n";
	  print LOG "maximum % of ambiguous (N) bases allowed in a read of the minimum length: $maximum_pc_of_n_bases\n";
	  print LOG "number of bases overlapping with an adapter sequence before any trimming occurs? $num_of_bases_overlapping_with_adapter_before_trimming\n";
	  print LOG "run initial FastQC (recommended)? $run_initial_fastqc\n";
	  print LOG "run TrimGalore (recommended)? $run_trimgalore\n";
	  print LOG "run Trimmomatic (recommended)? $run_trimmomatic\n";
	  print LOG "run post-processing FastQC (recommended)? $run_final_fastqc\n";
	  print LOG "will we use a given reference genome, if specified? $use_the_provided_ref_genome (if no, we will run Kaiju & Mash instead. If a genome is specified but unavailable, we will fail)\n";
	  print LOG "will we run Kaiju using database P (if using multiple databases, Kaiju will merge output)? $run_kaiju_on_db_p\n";
	  print LOG "will we run Kaiju using database E (if using multiple databases, Kaiju will merge output)? $run_kaiju_on_db_e\n";
	  print LOG "will we run Kaiju independently on PE and SE reads, as an internal test of the consistency of classification? (de facto SE reads being created by Trimmomatic) $run_kaiju_on_se_reads\n";
	  print LOG "will we fail if Kaiju cannot assign >= $min_pc_of_reads_for_being_confident_in_order_classification% of reads to the rank of order? $fail_beyond_order\n";
	  print LOG "will the VCF record calls for all sites? $output_all_sites (if no, VCF will only record calls at variant sites)\n";
	  print LOG "regularise VCF (recommended)? $regularise_vcf\n";
	  print LOG "filter VCF to remove low-confidence variants (recommended)? $filter_vcf\n";
	  print LOG "will a de novo assembly be made of the QC'd, cleaned and classified reads? $create_assembly\n";
	  print LOG "QC the de novo assembly (requires that it be created)? $qc_assembly\n";
	  print LOG "will we annotate the de novo assembly (requires that it be created)? $annotate_assembly\n";
	  print LOG "run TETyper? $run_tetyper\n";
	  print LOG "are we keeping the BAM? $keep_bam\n";
	  print LOG "are we keeping the cleaned fqs? $keep_cleaned_fq\n";

	  # (2) RECORD THE VERSION NUMBERS OF EACH PIECE OF SOFTWARE IN THE PIPELINE, IN THE APPROXIMATE ORDER IN WHICH THEY ARE USED (WE KNOW THESE COMMANDS CANNOT FAIL AS WE TESTED FOR THE EXISTENCE OF EACH VARIABLE EARLIER)
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
	  print SH "echo \"\n### McCortex ###\n\" >> $log_file\n";
	  print SH "$mccortex_path -h &>> $log_file\n";
	  print SH "echo \"\n### Atlas ###\n\" >> $log_file\n";
	  print SH "$atlas_path --version &>> $log_file\n";
	  print SH "echo \"\n### Mykrobe ###\n\" >> $log_file\n";
	  print SH "$mykrobe_path --version &>> $log_file\n";
	  print SH "echo \"\n### Mash ###\n\" >> $log_file\n";
	  print SH "$mash_path --version &>> $log_file\n";
	  print SH "echo \"\n### seqtk ###\n\" >> $log_file\n";
	  print SH "$seqtk_path &>> $log_file\n";
	  print SH "echo \"\n### Entrez Direct ###\n\" >> $log_file\n";
	  print SH "$efetch_path --help &>> $log_file\n";
	  print SH "echo \"\n### BWA ###\n\" >> $log_file\n";
	  print SH "$bwa_path &>> $log_file\n";
	  print SH "echo \"\n### sambamba ###\n\" >> $log_file\n";
	  print SH "$sambamba_path --version &>> $log_file\n";
	  print SH "echo \"\n### SAMtools ###\n\" >> $log_file\n";
	  print SH "$samtools_path &>> $log_file\n";
	  print SH "echo \"\n### BCFtools ###\n\" >> $log_file\n";
	  print SH "$bcftools_path --version &>> $log_file\n";
	  print SH "echo \"\n### Picard Tools ###\n\" >> $log_file\n";
	  print SH "java -jar $picard_path MarkDuplicates --version &>> $log_file\n";
	  print SH "echo \"\n### VCFlib ###\n\" >> $log_file\n";
	  print SH "$vcfoverlay_path --version &>> $log_file\n"; # $vcfoverlay_path is chosen because version numbers are not available for $vcffilter_path or $vcfallelicprimitives_path; see https://github.com/vcflib/vcflib/issues/102
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
	  print SH "echo \"\n### Prokka ###\n\" >> $log_file\n";
	  print SH "$prokka_path --version &>> $log_file\n";
	  
	  print SH "echo \"\n### TETyper ###\n\" >> $log_file\n";
	  print SH "$tetyper_path -h &>> $log_file\n";
	  print SH "echo \"\n### BLAST+ ###\n\" >> $log_file\n";
	  print SH "$blastn_path -version &>> $log_file\n";
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"\n-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. confirming database availability...\" >> $log_file\n";	  
	  print SH "echo \"\n### Database employed by Centrifuge ###\n\" >> $log_file\n";
	  print SH "date_last_modified=\"`date -r $centrifuge_db`\"\n";
	  print SH "echo \"$centrifuge_db. last modified \$date_last_modified\" >> $log_file\n";
	  print SH "echo \"\n### Databases employed by Kaiju ###\n\" >> $log_file\n";
	  print SH "date_last_modified=\"`date -r $kaiju_db_p`\"\n";
	  print SH "echo \"$kaiju_db_p. last modified \$date_last_modified\" >> $log_file\n";
	  print SH "date_last_modified=\"`date -r $kaiju_db_e`\"\n";
	  print SH "echo \"$kaiju_db_e. last modified \$date_last_modified\" >> $log_file\n";
	  
	  # (3) FORMALLY INITIATE PIPELINE
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"\n-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. beginning pipeline...\n\" >> $log_file\n";
	  print SH "set -e\n";
	  print SH "trap 'last_command=\$current_command; current_command=\$BASH_COMMAND' DEBUG\n";
	  print SH "trap 'echo \"the last command run - \${last_command} - finished with exit code \$?.\"' EXIT\n";
	  
	  # (4) CREATE AN OUTPUT DIRECTORY AND WITHIN IT, DOWNLOAD THE RAW DATA
	  ## (4i) download reads from the ENA, unless we're instead running a systems test using locally stored benchmarking datasets
	  my $ena1 = ''; my $ena2 = ''; my $url1 = ''; my $url2 = '';
	  if ($run_test_set == 0)
		{ my $first_3; my $first_6; my $digits;
		  if ($run_id =~ /^(.{3}).*?$/) { $first_3 = $1; }
		  if ($run_id =~ /^(.{6}).*?$/) { $first_6 = $1; }
		  if ($run_id =~ /^.*?(\d+)$/)  { $digits  = $1; }
		  my $number_of_digits = length($digits);
		  if ($number_of_digits == 6)
			{ $ena1 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$run_id/$run_id"."_1.fastq.gz";
			  $ena2 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$run_id/$run_id"."_2.fastq.gz";
			  $url1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$run_id/$run_id"."_1.fastq.gz";
			  $url2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$run_id/$run_id"."_2.fastq.gz";
			}
		  elsif ($number_of_digits == 7)
			{ if ($digits =~ /^.+?(\d{1})$/)
				{ my $last_digit = $1;
				  $ena1 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_1.fastq.gz";
				  $ena2 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_2.fastq.gz";
				  $url1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_1.fastq.gz";
				  $url2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_2.fastq.gz";
				}
			}
		  elsif ($number_of_digits == 8)
			{ if ($digits =~ /^.+?(\d{2})$/)
				{ my $last_two_digits = $1;
				  $ena1 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_1.fastq.gz";
				  $ena2 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_2.fastq.gz";
				  $url1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_1.fastq.gz";
				  $url2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_2.fastq.gz";
				}
			}
		  elsif ($number_of_digits == 9)
			{ if ($digits =~ /^.+?(\d{3})$/)
				{ my $last_three_digits = $1;
				  $ena1 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_1.fastq.gz";
				  $ena2 = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_2.fastq.gz";
				  $url1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_1.fastq.gz";
				  $url2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_2.fastq.gz";
				}
			}
		  if (($ena1 eq '') or ($ena2 eq '2'))
			{ print LOG "WARNING: unable to parse SRA run ID $run_id; this will be skipped\n";
			}
		}
	  next if ( ($run_test_set == 0) && (($ena1 eq '') or ($ena2 eq '2')) ); # CHECKPOINT: skip unrecognisable SRA run IDs
	  print SH "mkdir -p $out_dir/$run_id\n";
	  print SH "cd $out_dir/$run_id\n";
	  my $fq_1 = ''; my $fq_2 = '';
	  if ($run_test_set == 0)
		{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. downloading fastqs for $run_id...\n\" >> $log_file\n";
		  # using ascp
		  #print SH "if curl --head --fail --silent \"$url1\" >/dev/null; then $ascp_path -QT -l 100m -P33001 -i $aspera_openssh_path $ena1 $out_dir/$run_id &>> $log_file; else echo \"ERROR: unable to find file at URL $url1; is this because $run_id does not contain PE fqs?\" >> $log_file && exit 1; fi\n"; # CHECKPOINT: fail if this file does not exist on the ENA, and so ascp will be unable to download fqs for this run ID (see https://stackoverflow.com/questions/12199059/how-to-check-if-an-url-exists-with-the-shell-and-probably-curl). Note also that the ascp command only obtains PAIRED-END reads (by contrast, single-end read files take the form "$run_id.fastq.gz")
		  #print SH "if curl --head --fail --silent \"$url2\" >/dev/null; then $ascp_path -QT -l 100m -P33001 -i $aspera_openssh_path $ena2 $out_dir/$run_id &>> $log_file; else echo \"ERROR: unable to find file at URL $url2; is this because $run_id does not contain PE fqs?\" >> $log_file && exit 1; fi\n";
		  # using wget
		  print SH "if curl --head --fail --silent \"$url1\" >/dev/null; then wget $url1 &>> $log_file; else echo \"ERROR: unable to find file at URL $url1; is this because $run_id does not contain PE fqs?\" >> $log_file && exit 1; fi\n"; # CHECKPOINT: fail if this file does not exist on the ENA, and so ascp will be unable to download fqs for this run ID (see https://stackoverflow.com/questions/12199059/how-to-check-if-an-url-exists-with-the-shell-and-probably-curl). Note also that the ascp command only obtains PAIRED-END reads (by contrast, single-end read files take the form "$run_id.fastq.gz")
		  print SH "if curl --head --fail --silent \"$url2\" >/dev/null; then wget $url2 &>> $log_file; else echo \"ERROR: unable to find file at URL $url2; is this because $run_id does not contain PE fqs?\" >> $log_file && exit 1; fi\n";
		  $fq_1 = "$out_dir/$run_id/$run_id"."_1.fastq.gz";
		  $fq_2 = "$out_dir/$run_id/$run_id"."_2.fastq.gz";
		  print SH "if [ ! -e \$fq_1 ]; then echo \"ERROR: expected PE reads but unable to find $fq_1\" >> $log_file && exit 1; fi\n"; # CHECKPOINT: fail if we have not downloaded files corresponding to paired-end reads
		  print SH "if [ ! -e \$fq_2 ]; then echo \"ERROR: expected PE reads but unable to find $fq_2\" >> $log_file && exit 1; fi\n";
		}
	  elsif ($run_test_set == 1)
		{ my $fq_1_path = $sra_run_ids{$run_id}{read1_path};
		  my $fq_2_path = $sra_run_ids{$run_id}{read2_path};
		  $fq_1 = "$out_dir/$run_id/$run_id"."_1.fq.gz";
		  $fq_2 = "$out_dir/$run_id/$run_id"."_2.fq.gz";
		  print SH "cp $fq_1_path $fq_1\n";
		  print SH "cp $fq_2_path $fq_2\n";
		}
	  # (4ii) if the user wishes to provide their own reference genome for alignment (plus GFF for assembly and FAA for assembly annotation), then we will obtain them now. These should be gzipped files from NCBI Genome, e.g. under the 'genome', 'protein' and 'gff' tabs of https://www.ncbi.nlm.nih.gov/genome/?term=mycobacterium+tuberculosis%5Borgn%5D. If no URLs are provided, we will resort to the default option: of using Kaiju & Mash to predict (and download) a reference.
	  my $ref = ''; my $gff = ''; my $faa = '';
	  if ($use_the_provided_ref_genome eq 'yes')
		{ my $provided_genome = ''; my $provided_gff = ''; my $provided_faa = '';
		  if (exists($sra_run_ids{$run_id}{fa}))  { $provided_genome = $sra_run_ids{$run_id}{fa};  }
		  if (exists($sra_run_ids{$run_id}{gff})) { $provided_gff    = $sra_run_ids{$run_id}{gff}; }
		  if (exists($sra_run_ids{$run_id}{faa})) { $provided_faa    = $sra_run_ids{$run_id}{faa}; }
		  if (($provided_genome ne '') and ($provided_gff ne '') and ($provided_faa ne ''))
			{ print SH "if curl --head --fail --silent \"$provided_genome\" >/dev/null; then wget $provided_genome -O $out_dir/$run_id/ref.fa.gz &>> $log_file; else echo \"ERROR: a reference genome was specified but we were unable to find the file at URL $provided_genome\" >> $log_file && exit 1; fi\n";
			  print SH "if curl --head --fail --silent \"$provided_gff\" >/dev/null; then wget $provided_gff -O $out_dir/$run_id/ref.gff.gz &>> $log_file; else echo \"ERROR: a reference GFF was specified but we were unable to find the file at URL $provided_gff\" >> $log_file && exit 1; fi\n";
			  print SH "if curl --head --fail --silent \"$provided_faa\" >/dev/null; then wget $provided_faa -O $out_dir/$run_id/ref.faa.gz &>> $log_file; else echo \"ERROR: a reference FAA was specified but we were unable to find the file at URL $provided_faa\" >> $log_file && exit 1; fi\n";
			  print SH "gunzip $out_dir/$run_id/ref.fa.gz\n";
			  print SH "gunzip $out_dir/$run_id/ref.gff.gz\n";
			  print SH "gunzip $out_dir/$run_id/ref.faa.gz\n";
			  $ref = "$out_dir/$run_id/ref.fa";
			  $gff = "$out_dir/$run_id/ref.gff";
			  $faa = "$out_dir/$run_id/ref.faa";
			}
		  else
			{ $use_the_provided_ref_genome = 'no'; # while the user wished to use the provided genomes, none were specified. Consequently, we set this parameter to 'no' to ensure that we'll run Kaiju & Mash after all.
			}
		}
	  my $log_number = 1;
	  
	  # (5) QUALITY INSPECTION AND PRE-PROCESSING
	  # (5i) initial quality inspection with FastQC
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
	  # (5ii) first round of cleaning (to automatically predict, and then remove, adapter sequence): TrimGalore
	  # NOTE: Trimmomatic and cutadapt have complementary, but distinct, functions. While Trimmomatic can also remove adapters, it requires that the adapter sequence be specified. By contrast, TrimGalore (which employs cutadapt) automatically predicts the adapter - more useful when input is of unknown origin. TrimGalore can also filter reads based on N content: Trimmomatic cannot. Here we use TrimGalore for N and adapter removal only, with quality-trimming later performed by Trimmomatic.
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
	  # (5iii) second round of cleaning (for general read quality): Trimmomatic. This will (i) remove bases from the end of a read if they are below a Phred score of $phred_threshold, (ii) clip the read if the average Phred score within a 4bp sliding window advanced from the 5’ end falls below $phred_threshold, and (iii) impose a minimum read length of $minimum_read_length bp
	  # NOTE: given 2 fqs as input, Trimmomatic will output 4 files, 2 of which are for singleton reads: read 1 where read 2 failed, and read 2 where read 1 failed. These are combined into one file of de facto single-end reads, which may be aligned separately.
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
		  print SH "cat $run_id_fwd_single $run_id_rev_single > $fq_se\n";
		  print SH "rm $run_id_fwd_single $run_id_rev_single\n";
		  $log_number++;
		}
	  # (5iv) final quality inspection with FastQC
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
	  
	  # (6) CLASSIFICATION, PART 1: HUMAN READ DETECTION. ACCOMPLISHED BY RUNNING CENTRIFUGE AGAINST A HUMAN-SPECIFIC DATABASE, TAKING THE COMBINED SET OF PE AND SE READS AS INPUT (UNLIKE KAIJU, CENTRIFUGE CAN TAKE BOTH SETS OF READS AS INPUT SIMULTANEOUSLY)
	  # IMPORTANT: the Kaiju databases, below, are prepared using protein sequences of representative genomes. This means that if we were to interpret the task of human read removal as being one of removing Kaiju's unclassified reads (that is, removing human reads essentially as collateral damage) then the proportion removed will be large - not because they are (all) contaminant reads (see the removal step below), but because they come from untranslated regions/un-annotated regions of the genome. Discarding them from further downstream steps does not make sense. In doing so, we will be calling SNPs only in coding regions and the assemblies will be broken. Consequently, to predict human (contaminant) content, we use a nucleotide-based classifier, Centrifuge, directly.
	  # (6i) predict human and HIV read content
	  my $centrifuge_non_human_read_list = '';
	  if ($run_test_set == 0) # human read removal is irrelevant using the simulated test data as we have a priori knowledge that no human content is present
		{ my $centrifuge_report_log_number = $log_number;
		  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_centrifuge_report_on_number_of_human_reads\n";
		  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. running Centrifuge on the cleaned and adapter-trimmed reads from $run_id...\n\" >> $log_file\n";
		  my $centrifuge_summary 			 = "$out_dir/$run_id/logs/$log_number"."_centrifuge_report_on_number_of_human_reads/$run_id.number_of_human_reads.txt";
		  my $centrifuge_read_classification = "$out_dir/$run_id/logs/$log_number"."_centrifuge_report_on_number_of_human_reads/$run_id.centrifuge_read_classification";
		  $centrifuge_non_human_read_list    = "$out_dir/$run_id/logs/$log_number"."_centrifuge_report_on_number_of_human_reads/$run_id.non_human_read_list";
		  if ($fq_se ne '')
			{ print SH "$centrifuge_path -p $num_threads -x $centrifuge_index_prefix -1 $fq_1 -2 $fq_2 -U $fq_se --report-file $centrifuge_summary -S $centrifuge_read_classification &>> $log_file\n"; }
		  elsif ($fq_se eq '')
			{ print SH "$centrifuge_path -p $num_threads -x $centrifuge_index_prefix -1 $fq_1 -2 $fq_2 --report-file $centrifuge_summary -S $centrifuge_read_classification &>> $log_file\n"; }
		  print SH "awk '\$2==\"unclassified\" { print \$1 }' $centrifuge_read_classification > $centrifuge_non_human_read_list\n";
		  print SH "rm $centrifuge_read_classification\n" unless ($keep_read_classifications eq 'yes');
		  # (6ii) filter the Centrifuge output report to erase all mention of HIV content
		  my $centrifuge_summary_without_hiv = "$out_dir/$run_id/logs/$log_number"."_centrifuge_report_on_number_of_human_reads/$run_id.number_of_human_reads_excl_hiv.txt";
		  print SH "awk 'NR==1 { print \$0 }' $centrifuge_summary > $centrifuge_summary_without_hiv\n"; # we will always save the first line of the original Centrifuge output: it's just the header line
		  print SH "grep \"Homo sapiens\" $centrifuge_summary >> $centrifuge_summary_without_hiv\n";
		  print SH "mv $centrifuge_summary_without_hiv $centrifuge_summary\n";
		  $log_number++;
		}
	  
	  # (7) REMOVE READS CLASSIFIED AS HUMAN BY CENTRIFUGE. WE DO THIS BY SUBSETTING THE LIST OF 'UNCLASSIFIED' (THAT IS, NON-HUMAN) READS AFTER MAPPING THE READS TO A HUMAN-SPECIFIC CENTRIFUGE DATABASE.
	  if ($run_test_set == 0)
		{ my $fq_1_excl_human  = "$out_dir/$run_id/$run_id"."_1.excl_human.fastq.gz";
		  my $fq_2_excl_human  = "$out_dir/$run_id/$run_id"."_2.excl_human.fastq.gz";
		  my $fq_se_excl_human = "$out_dir/$run_id/$run_id.se.excl_human.fastq.gz";
		  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. removing reads classified as human from the cleaned and adapter-trimmed paired-end reads from $run_id...\n\" >> $log_file\n";
		  print SH "$seqtk_path subseq $fq_1 $centrifuge_non_human_read_list | gzip > $fq_1_excl_human\n";
		  print SH "$seqtk_path subseq $fq_2 $centrifuge_non_human_read_list | gzip > $fq_2_excl_human\n";
		  if ($fq_se ne '')
			{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
			  print SH "elapsed_time=\$SECONDS\n";
			  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. removing reads classified as human from the cleaned and adapter-trimmed single-end reads from $run_id...\n\" >> $log_file\n";
			  print SH "$seqtk_path subseq $fq_se $centrifuge_non_human_read_list > $fq_se_excl_human\n";
			}
		  print SH "rm $centrifuge_non_human_read_list\n";
		  print SH "mv $fq_1_excl_human $fq_1\n";
		  print SH "mv $fq_2_excl_human $fq_2\n";
		  if ($fq_se ne '')
			{ print SH "mv $fq_se_excl_human $fq_se\n"; }
		}
	  
	  # (8, 9, 10, 11, 12, 13) CLASSIFICATION, PART 2: SPECIATION. CREATE KRONA-FORMATTED KAIJU REPORTS FOR THE PE AND SE READS, ALONGSIDE PREDICTING ORDER, FAMILY, GENUS AND SPECIES
	  if (($run_test_set == 1) or ($use_the_provided_ref_genome eq 'no'))
		{ my $kaiju_report_log_number = $log_number;
		  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_kaiju_reports\n";
		  for(my $y=0;$y<=1;$y++)
			{ my $kaiju_db = ''; my $db_type = '';
			  if ($y == 0)
				{ $kaiju_db = $kaiju_db_p;
				  $db_type  = 'P';
				}
			  elsif ($y == 1)
				{ $kaiju_db = $kaiju_db_e;
				  $db_type  = 'E';
				}
			  next if (($db_type eq 'P') and ($run_kaiju_on_db_p eq 'no'));
			  next if (($db_type eq 'E') and ($run_kaiju_on_db_e eq 'no'));
			  for(my $x=0;$x<=1;$x++)
				{ my $pe_or_se = '';
				  if ($x == 0) { $pe_or_se = 'paired-end'; } elsif ($x == 1) { $pe_or_se = 'single-end'; }
				  next if ( (($pe_or_se eq 'single-end') and ($run_trimmomatic eq 'no')) or (($pe_or_se eq 'single-end') and ($run_kaiju_on_se_reads eq 'no')) );
				  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
				  print SH "elapsed_time=\$SECONDS\n";
				  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. running Kaiju (database $db_type) on the cleaned and adapter-trimmed $pe_or_se reads from $run_id...\n\" >> $log_file\n";
				  if ($pe_or_se eq 'paired-end')
					{ print SH "$kaiju_path -z $num_threads -t $nodes_dmp -f $kaiju_db -i $fq_1 -j $fq_2 -a greedy -e 5 -E 0.05 -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.$db_type.$pe_or_se.txt\n"; }
				  elsif ($pe_or_se eq 'single-end')
					{ print SH "$kaiju_path -z $num_threads -t $nodes_dmp -f $kaiju_db -i $fq_se -a greedy -e 5 -E 0.05 -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.$db_type.$pe_or_se.txt\n"; }
				  print SH "$kaijuReport_path -t $nodes_dmp -n $names_dmp -i $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.$db_type.$pe_or_se.txt -r species -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_species_prediction.$db_type.$pe_or_se.txt\n";
				  print SH "$kaijuReport_path -t $nodes_dmp -n $names_dmp -i $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.$db_type.$pe_or_se.txt -r genus -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_genus_prediction.$db_type.$pe_or_se.txt\n";
				  print SH "$kaijuReport_path -t $nodes_dmp -n $names_dmp -i $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.$db_type.$pe_or_se.txt -r family -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_family_prediction.$db_type.$pe_or_se.txt\n";
				  print SH "$kaijuReport_path -t $nodes_dmp -n $names_dmp -i $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.$db_type.$pe_or_se.txt -r order -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_order_prediction.$db_type.$pe_or_se.txt\n";
				  print SH "$kaijutokrona_path -t $nodes_dmp -n $names_dmp -i $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.$db_type.$pe_or_se.txt -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report_for_krona.$db_type.$pe_or_se.txt &>> $log_file\n";
				  print SH "$ktImportText_path -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.$db_type.$pe_or_se.html $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report_for_krona.$db_type.$pe_or_se.txt &>> $log_file\n";
				  print SH "rm $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report_for_krona.$db_type.$pe_or_se.txt\n";
				}
			}

		  # (9) COMBINE KAIJU OUTPUT FROM MULTIPLE DATABASES, ADJUDICATING CLASSIFICATIONS BASED ON THE LOWEST RANK IF THEY ARE WITHIN THE SAME LINEAGE, ELSE THE LOWEST COMMON ANCESTOR
		  if (($run_kaiju_on_db_p eq 'yes') and ($run_kaiju_on_db_e eq 'yes'))
			{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
			  print SH "elapsed_time=\$SECONDS\n";
			  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. combining multiple Kaiju reports from $run_id...\n\" >> $log_file\n";
			  my $report_P_pe = "$out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.P.paired-end.txt";
			  my $report_E_pe = "$out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.E.paired-end.txt";
			  my $combined_report_pe = "$out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.paired-end.txt";
			  my $combined_report_se = '';
			  print SH "$mergeOutputs_path -t $nodes_dmp -c lowest -i <(sort -k2,2 $report_P_pe) -j <(sort -k2,2 $report_E_pe) -o $combined_report_pe\n";
			  print SH "$kaijuReport_path -t $nodes_dmp -n $names_dmp -i $combined_report_pe -r species -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_species_prediction.paired-end.txt\n";
			  print SH "$kaijuReport_path -t $nodes_dmp -n $names_dmp -i $combined_report_pe -r genus -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_genus_prediction.paired-end.txt\n";
			  print SH "$kaijuReport_path -t $nodes_dmp -n $names_dmp -i $combined_report_pe -r family -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_family_prediction.paired-end.txt\n";
			  print SH "$kaijuReport_path -t $nodes_dmp -n $names_dmp -i $combined_report_pe -r order -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_order_prediction.paired-end.txt\n";
			  print SH "rm $report_P_pe $report_E_pe\n" unless ($keep_read_classifications eq 'yes');
			  if (($fq_se ne '') and ($run_kaiju_on_se_reads eq 'yes'))
				{ my $report_P_se = "$out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.P.single-end.txt";
				  my $report_E_se = "$out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.E.single-end.txt";
				  $combined_report_se = "$out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.single-end.txt";
				  print SH "$mergeOutputs_path -t $nodes_dmp -c lowest -i <(sort -k2,2 $report_P_se) -j <(sort -k2,2 $report_E_se) -o $combined_report_se\n";
				  print SH "$kaijuReport_path -t $nodes_dmp -n $names_dmp -i $combined_report_se -r species -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_species_prediction.single-end.txt\n";
				  print SH "$kaijuReport_path -t $nodes_dmp -n $names_dmp -i $combined_report_se -r genus -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_genus_prediction.single-end.txt\n";
				  print SH "$kaijuReport_path -t $nodes_dmp -n $names_dmp -i $combined_report_se -r family -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_family_prediction.single-end.txt\n";
				  print SH "$kaijuReport_path -t $nodes_dmp -n $names_dmp -i $combined_report_se -r order -o $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_order_prediction.single-end.txt\n";
				  print SH "rm $report_P_se $report_E_se\n" unless ($keep_read_classifications eq 'yes');
				}
			  print SH "rm $combined_report_pe\n" unless ($keep_read_classifications eq 'yes');
			  if (($fq_se ne '') and ($run_kaiju_on_se_reads eq 'yes'))
				{ print SH "rm $combined_report_se\n" unless ($keep_read_classifications eq 'yes'); }
			}
		  elsif (($run_kaiju_on_db_p eq 'yes') and ($run_kaiju_on_db_e eq 'no'))
			{ for(my $x=0;$x<=1;$x++)
				{ my $pe_or_se = '';
				  if ($x == 0) { $pe_or_se = 'paired-end'; } elsif ($x == 1) { $pe_or_se = 'single-end'; }
				  next if ( (($pe_or_se eq 'single-end') and ($run_trimmomatic eq 'no')) or (($pe_or_se eq 'single-end') and ($run_kaiju_on_se_reads eq 'no')) );
				  print SH "mv $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_species_prediction.P.$pe_or_se.txt $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_species_prediction.$pe_or_se.txt\n";
				  print SH "mv $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_genus_prediction.P.$pe_or_se.txt $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_genus_prediction.$pe_or_se.txt\n";
				  print SH "mv $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_family_prediction.P.$pe_or_se.txt $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_family_prediction.$pe_or_se.txt\n";
				  print SH "mv $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_order_prediction.P.$pe_or_se.txt $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_order_prediction.$pe_or_se.txt\n";
				  print SH "mv $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.P.$pe_or_se.html $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.$pe_or_se.html\n";
				  print SH "mv $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.P.$pe_or_se.txt $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.$pe_or_se.txt\n" unless ($keep_read_classifications eq 'yes');
				}
			}
		  elsif (($run_kaiju_on_db_p eq 'no') and ($run_kaiju_on_db_e eq 'yes'))
			{ for(my $x=0;$x<=1;$x++)
				{ my $pe_or_se = '';
				  if ($x == 0) { $pe_or_se = 'paired-end'; } elsif ($x == 1) { $pe_or_se = 'single-end'; }
				  next if ( (($pe_or_se eq 'single-end') and ($run_trimmomatic eq 'no')) or (($pe_or_se eq 'single-end') and ($run_kaiju_on_se_reads eq 'no')) );
				  print SH "mv $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_species_prediction.E.$pe_or_se.txt $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_species_prediction.$pe_or_se.txt\n";
				  print SH "mv $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_genus_prediction.E.$pe_or_se.txt $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_genus_prediction.$pe_or_se.txt\n";
				  print SH "mv $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_family_prediction.E.$pe_or_se.txt $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_family_prediction.$pe_or_se.txt\n";
				  print SH "mv $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_order_prediction.E.$pe_or_se.txt $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_order_prediction.$pe_or_se.txt\n";
				  print SH "mv $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.E.$pe_or_se.html $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.$pe_or_se.html\n";
				  print SH "mv $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.E.$pe_or_se.txt $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.$pe_or_se.txt\n" unless ($keep_read_classifications eq 'yes');
				}
			}
		  if (($run_kaiju_on_db_p eq 'no') or ($run_kaiju_on_db_e eq 'no'))
			{ for(my $x=0;$x<=1;$x++)
				{ my $pe_or_se = '';
				  if ($x == 0) { $pe_or_se = 'paired-end'; } elsif ($x == 1) { $pe_or_se = 'single-end'; }
				  next if ( (($pe_or_se eq 'single-end') and ($run_trimmomatic eq 'no')) or (($pe_or_se eq 'single-end') and ($run_kaiju_on_se_reads eq 'no')) );
				  print SH "rm $out_dir/$run_id/logs/$log_number"."_kaiju_reports/$run_id.kaiju_report.$pe_or_se.txt\n" unless ($keep_read_classifications eq 'yes');
				}
			}
		  $log_number++;

		  # (10) SELECT, DOWNLOAD AND INDEX AN APPROPRIATE REFERENCE GENOME. THIS TAKES THE BEST HIT PREDICTION FOR A GIVEN TAXONOMIC RANK AND THEN DOWNLOADS ALL GENOMES (FROM NCBI ASSEMBLY: https://www.ncbi.nlm.nih.gov/pubmed/26578580) WITHIN A CERTAIN TAXONOMIC 'SPACE'. AFTER THIS, A SPECIFIC REFERENCE GENOME IS SELECTED BASED ON THE CLOSEST MASH DISTANCE BETWEEN THESE GENOMES AND THE READS.
		  print SH "top_species_hit_pe=\$(awk -F '\\t' 'NR==3 {print \$3}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_species_prediction.paired-end.txt)\n";
		  print SH "pc_top_species_hit_pe=\$(awk -F '\\t' 'NR==3 {print \$1}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_species_prediction.paired-end.txt)\n";
		  if (($fq_se ne '') and ($run_kaiju_on_se_reads eq 'yes'))
			{ print SH "top_species_hit_se=\$(awk -F '\\t' 'NR==3 {print \$3}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_species_prediction.single-end.txt)\n";
			  print SH "pc_top_species_hit_se=\$(awk -F '\\t' 'NR==3 {print \$1}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_species_prediction.single-end.txt)\n";
			  print SH "if [ \"\$top_species_hit_pe\" != \"\$top_species_hit_se\" ]; then echo \"ERROR: candidate species for $run_id as predicted by top hit for PE reads (\$top_species_hit_pe) != top hit for SE reads (\$top_species_hit_se); unable to select a reference genome and proceed with alignment\n\" >> $log_file && exit 1; fi\n";
			} # CHECKPOINT: fail if there are discordant Kaiju species predictions (only applicable if there are PE and de-facto SE reads for comparison)
		  print SH "top_genus_hit_pe=\$(awk -F '\\t' 'NR==3 {print \$3}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_genus_prediction.paired-end.txt)\n";
		  print SH "pc_top_genus_hit_pe=\$(awk -F '\\t' 'NR==3 {print \$1}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_genus_prediction.paired-end.txt)\n";
		  if (($fq_se ne '') and ($run_kaiju_on_se_reads eq 'yes'))
			{ print SH "top_genus_hit_se=\$(awk -F '\\t' 'NR==3 {print \$3}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_genus_prediction.single-end.txt)\n";
			  print SH "pc_top_genus_hit_se=\$(awk -F '\\t' 'NR==3 {print \$1}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_genus_prediction.single-end.txt)\n";
			  print SH "if [ \"\$top_genus_hit_pe\" != \"\$top_genus_hit_se\" ]; then echo \"ERROR: candidate genus for $run_id as predicted by top hit for PE reads (\$top_genus_hit_pe) != top hit for SE reads (\$top_genus_hit_se); unable to select a reference genome and proceed with alignment\n\" >> $log_file && exit 1; fi\n";
			} # CHECKPOINT: fail if there are discordant Kaiju genus predictions (only applicable if there are PE and de-facto SE reads for comparison)
		  print SH "top_family_hit_pe=\$(awk -F '\\t' 'NR==3 {print \$3}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_family_prediction.paired-end.txt)\n";
		  print SH "pc_top_family_hit_pe=\$(awk -F '\\t' 'NR==3 {print \$1}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_family_prediction.paired-end.txt)\n";
		  if (($fq_se ne '') and ($run_kaiju_on_se_reads eq 'yes'))
			{ print SH "top_family_hit_se=\$(awk -F '\\t' 'NR==3 {print \$3}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_family_prediction.single-end.txt)\n";
			  print SH "pc_top_family_hit_se=\$(awk -F '\\t' 'NR==3 {print \$1}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_family_prediction.single-end.txt)\n";
			  print SH "if [ \"\$top_family_hit_pe\" != \"\$top_family_hit_se\" ]; then echo \"ERROR: candidate family for $run_id as predicted by top hit for PE reads (\$top_family_hit_pe) != top hit for SE reads (\$top_family_hit_se); unable to select a reference genome and proceed with alignment\n\" >> $log_file && exit 1; fi\n";
			} # CHECKPOINT: fail if there are discordant Kaiju family predictions (only applicable if there are PE and de-facto SE reads for comparison)
		  print SH "top_order_hit_pe=\$(awk -F '\\t' 'NR==3 {print \$3}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_order_prediction.paired-end.txt)\n";
		  print SH "pc_top_order_hit_pe=\$(awk -F '\\t' 'NR==3 {print \$1}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_order_prediction.paired-end.txt)\n";
		  if (($fq_se ne '') and ($run_kaiju_on_se_reads eq 'yes'))
			{ print SH "top_order_hit_se=\$(awk -F '\\t' 'NR==3 {print \$3}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_order_prediction.single-end.txt)\n";
			  print SH "pc_top_order_hit_se=\$(awk -F '\\t' 'NR==3 {print \$1}' $out_dir/$run_id/logs/$kaiju_report_log_number"."_kaiju_reports/$run_id.kaiju_order_prediction.single-end.txt)\n";
			  print SH "if [ \"\$top_order_hit_pe\" != \"\$top_order_hit_se\" ]; then echo \"ERROR: candidate order for $run_id as predicted by top hit for PE reads (\$top_order_hit_pe) != top hit for SE reads (\$top_order_hit_se); unable to select a reference genome and proceed with alignment\n\" >> $log_file && exit 1; fi\n";
			} # CHECKPOINT: fail if there are discordant Kaiju order predictions (only applicable if there are PE and de-facto SE reads for comparison)
		  
		  # (11) DETERMINE HOW BROAD THE SEARCH SPACE FOR MASH WILL BE, I.E. HOW MANY GENOMES WILL BE DOWNLOADED OF A SPECIFIC TAXONOMIC RANK. WE DO THIS BY REQUIRING THAT AN ABOVE-THRESHOLD % OF READS BE CLASSIFIED AT A GIVEN TAXONOMIC RANK. SHOULD THIS NOT BE MET AT THE LOWEST RANK (SPECIES), WE TRY AGAIN AT THE SECOND (GENUS), THIRD (FAMILY), AND FOURTH (ORDER), BEFORE ABORTING. NOTE THAT THE TOP HIT FOR SEARCHING WILL, AT THIS POINT, BE EITHER A SPECIES OR GENUS-LEVEL ID, IRRESPECTIVE OF WHETHER THE RANK IN WHICH WE ARE ACTUALLY CONFIDENT IS A HIGHER ONE, FAMILY OR ORDER. THIS IS BECAUSE SPECIES AND GNEUS ARE THE ONLY RANKS WE CAN *IMMEDIATELY* EXTRACT FROM THE DOWNLOADED NCBI QUERY FILES. WHERE NECESSARY, WE SHALL ASCEND IN RANKS LATER BY REFERENCE TO PRE-COMPUTED TAXONOMIC LOOKUP TABLES, FINDING THE RIGHT FAMILY AND ORDER IDS.
		  print SH "top_hit_for_searching='unknown'\n";
		  print SH "taxonomic_rank_in_which_we_have_highest_confidence='unknown'\n";
		  print SH "if (( \$(echo \"\$pc_top_species_hit_pe > $min_pc_of_reads_for_being_confident_in_species_classification\" | bc -l) )); then top_hit_for_searching=\"\$top_species_hit_pe\" && taxonomic_rank_in_which_we_have_highest_confidence\='species'; elif (( \$(echo \"\$pc_top_genus_hit_pe > $min_pc_of_reads_for_being_confident_in_genus_classification\" | bc -l) )); then top_hit_for_searching=\"\$top_genus_hit_pe\" && taxonomic_rank_in_which_we_have_highest_confidence\='genus'; elif (( \$(echo \"\$pc_top_family_hit_pe > $min_pc_of_reads_for_being_confident_in_family_classification\" | bc -l) )); then top_hit_for_searching=\"\$top_genus_hit_pe\" && taxonomic_rank_in_which_we_have_highest_confidence\='family'; elif (( \$(echo \"\$pc_top_order_hit_pe > $min_pc_of_reads_for_being_confident_in_order_classification\" | bc -l) )); then top_hit_for_searching=\"\$top_genus_hit_pe\" && taxonomic_rank_in_which_we_have_highest_confidence\='order'; fi\n";
		  if ($fail_beyond_order eq 'yes')
			{ print SH "if [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" = 'unknown' ]; then echo \"ERROR: too low a proportion of reads can be classified to a single taxonomic rank; unable to select a reference genome and proceed with alignment\" >> $log_file && exit 1; fi\n"; # CHECKPOINT: fail if too low a % of reads can be assigned to the broadest taxonomic rank
			}
		  elsif ($fail_beyond_order eq 'no')
			{ print SH "if [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" = 'unknown' ]; then taxonomic_rank_in_which_we_have_highest_confidence\='order'; fi\n"; }
		  print SH "GENOME=\$($esearch_path -db genome -query \"\$top_hit_for_searching\"[orgn] | $efetch_path -format docsum | tee \"$out_dir/$run_id/\$top_hit_for_searching.genome.esearch.docsum\")\n";
		  print SH "ACC=`echo \$GENOME | $xtract_path -pattern DocumentSummary -element Assembly_Accession`\n";
		  print SH "ACC=\$(echo \$ACC | awk 'BEGIN {FS=\" \";} { print \$1 }')\n";
		  print SH "RESULT=\$($esearch_path -db assembly -query \"\$ACC\" | $efetch_path -format docsum | tee \"$out_dir/$run_id/\$top_hit_for_searching.assembly.esearch.docsum\")\n";
		  print SH "SPECIESTAXID=`echo \$RESULT | $xtract_path -pattern DocumentSummary -element SpeciesTaxid`\n"; # see ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt for a description of the difference between Taxid and Speciestaxid
		  print SH "rm \"$out_dir/$run_id/\$top_hit_for_searching.genome.esearch.docsum\" \"$out_dir/$run_id/\$top_hit_for_searching.assembly.esearch.docsum\"\n";
		  $ref = "$out_dir/$run_id/ref.fa";
		  $gff = "$out_dir/$run_id/ref.gff";
		  $faa = "$out_dir/$run_id/ref.faa";
		  
		  # (12) USE THE KAIJU BEST SPECIES/GENUS/FAMILY/ORDER HIT TO NARROW DOWN A SEARCH OF "SPECIES SPACE" (OR, MORE BROADLY, "GENUS SPACE"/"FAMILY SPACE"/"ORDER SPACE"). WE DO THIS BY DOWNLOADING THE LATEST COMPLETE GENOMES FOR EACH SPECIES/GENUS/FAMILY/ORDER, CREATING MASH SKETCHES FOR THEM AND THEN PICKING THE SKETCH WITH SHORTEST DISTANCE TO THE READS. TO AVOID RE-DOWNLOADING GENOMES UNNECESSARILY, MASH SKETCHES (WHICH ARE SMALL BINARY REPRESENTATIONS) CAN BE STORED LOCALLY - THESE ALL HAVE UNIQUE NCBI ASSEMBLY ACCESSIONS.
		  # this concept can be illustrated using fig. 3 of the Mash paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x. "Species space" can be considered the higher-order clustering of nodes (species) in a network graph, all of which are assigned the same (species) name. Consequently, a Kaiju 'best species hit' narrows the search to within this region only. Mash will pinpoint the closest species to the reads within this region.
		  # NOTE: see question 15 ("how can I download RefSeq data for all complete bacterial genomes?") of this FAQ: https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/
		  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_mash_genome_prediction\n";
		  my $ref_genome_download_log = "$out_dir/$run_id/logs/$log_number"."_mash_genome_prediction/$run_id.reference_genome_download_log.txt";
		  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. downloading the set of complete genomes corresponding to the best hit taxonomy ID...\n\" >> $log_file\n";
		  print SH "if curl --head --fail --silent ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt >/dev/null; then wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt &>> $log_file; else echo \"ERROR: unable to download ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt; cannot proceed with genome download and mash clustering\n\" >> $log_file && exit 1; fi\n"; # CHECKPOINT: fail if we cannot download assembly_summary.txt (as this would render genome download and mash clustering impossible)
		  print SH "date_last_modified=\"`date -r assembly_summary.txt`\"\n";
		  print SH "echo \"--> downloaded ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt. last modified \$date_last_modified\" >> $log_file\n";
		  print SH "GENUSTAXID=\$(awk -v speciestaxid=\$SPECIESTAXID '\$1==speciestaxid { print \$2 }' $species_id_to_genus_id_lookup)\n";
		  print SH "FAMILYTAXID=\$(awk -v speciestaxid=\$SPECIESTAXID '\$1==speciestaxid { print \$2 }' $species_and_genus_id_to_family_id_lookup)\n";
		  print SH "ORDERTAXID=\$(awk -v speciestaxid=\$SPECIESTAXID '\$1==speciestaxid { print \$2 }' $species_and_genus_id_to_order_id_lookup)\n";
		  print SH "echo \"TOP KAIJU HIT (SPECIES): \$top_species_hit_pe (\$pc_top_species_hit_pe % of PE reads); taxon ID \$SPECIESTAXID\" >> $ref_genome_download_log\n";
		  print SH "echo \"TOP KAIJU HIT (GENUS): \$top_genus_hit_pe (\$pc_top_genus_hit_pe % of PE reads); taxon ID \$GENUSTAXID\" >> $ref_genome_download_log\n";
		  print SH "echo \"TOP KAIJU HIT (FAMILY): \$top_family_hit_pe (\$pc_top_family_hit_pe % of PE reads); taxon ID \$FAMILYTAXID\" >> $ref_genome_download_log\n";
		  print SH "echo \"TOP KAIJU HIT (ORDER): \$top_order_hit_pe (\$pc_top_order_hit_pe % of PE reads); taxon ID \$ORDERTAXID\" >> $ref_genome_download_log\n";
		  print SH "echo \"LOWEST TAXONOMIC RANK IN WHICH WE HAVE HIGHEST CONFIDENCE (i.e. AN ABOVE-THRESHOLD % OF READS)\: \$taxonomic_rank_in_which_we_have_highest_confidence\" >> $ref_genome_download_log\n";
		  print SH "ALL_SPECIES_THIS_GENUS=\$(awk -v genustaxid=\$GENUSTAXID '\$2==genustaxid { print \$1 }' $species_id_to_genus_id_lookup)\n"; # if taxonomic_rank_in_which_we_have_highest_confidence eq 'genus', identify the set of species that also have this ID
		  print SH "ALL_SPECIES_THIS_FAMILY=\$(awk -v familytaxid=\$FAMILYTAXID '\$2==familytaxid { print \$1 }' $species_and_genus_id_to_family_id_lookup)\n"; # if taxonomic_rank_in_which_we_have_highest_confidence eq 'family', identify the set of species that also have this ID
		  print SH "ALL_SPECIES_THIS_ORDER=\$(awk -v ordertaxid=\$ORDERTAXID '\$2==ordertaxid { print \$1 }' $species_and_genus_id_to_order_id_lookup)\n"; # if taxonomic_rank_in_which_we_have_highest_confidence eq 'order', identify the set of species that also have this ID
		  print SH "ALL_SPECIES_THIS_GENUS_ARRAY=(\$ALL_SPECIES_THIS_GENUS)\n";
		  print SH "ALL_SPECIES_THIS_FAMILY_ARRAY=(\$ALL_SPECIES_THIS_FAMILY)\n";
		  print SH "ALL_SPECIES_THIS_ORDER_ARRAY=(\$ALL_SPECIES_THIS_ORDER)\n";
		  print SH "if [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" == 'species' ]; then echo \"NCBI TAXON ID FOR LOWEST TAXONOMIC RANK ASSIGNED AN ABOVE-THRESHOLD % OF READS: \$SPECIESTAXID\" >> $ref_genome_download_log; elif [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" == 'genus' ]; then echo \"NCBI TAXON ID FOR LOWEST TAXONOMIC RANK ASSIGNED AN ABOVE-THRESHOLD % OF READS: \$GENUSTAXID\" >> $ref_genome_download_log; elif [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" == 'family' ]; then echo \"NCBI TAXON ID FOR LOWEST TAXONOMIC RANK ASSIGNED AN ABOVE-THRESHOLD % OF READS: \$FAMILYTAXID\" >> $ref_genome_download_log; elif [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" == 'order' ]; then echo \"NCBI TAXON ID FOR LOWEST TAXONOMIC RANK ASSIGNED AN ABOVE-THRESHOLD % OF READS: \$ORDERTAXID\" >> $ref_genome_download_log; fi\n";
		  # we now download the latest complete genomes for the corresponding taxon ID. "Complete genome" is the highest level of assembly, whereby "all chromosomes are gapless and have no runs of 10 or more ambiguous bases (Ns), there are no unplaced or unlocalized scaffolds, and all the expected chromosomes are present (i.e. the assembly is not noted as having partial genome representation). Plasmids and organelles may or may not be included in the assembly but if present then the sequences are gapless" (ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt).
		  print SH "if [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" == 'species' ]; then awk -v speciestaxid=\$SPECIESTAXID -F \"\\t\" '\$7==speciestaxid && \$12==\"Complete Genome\" && \$11==\"latest\"{print \$20}' $out_dir/$run_id/assembly_summary.txt >> $out_dir/$run_id/ftpdirpaths.txt; elif [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" == 'genus' ]; then for i in \${ALL_SPECIES_THIS_GENUS[@]}; do echo \"identifying species IDs corresponding to genus ID \$GENUSTAXID: \$i...\" >> $log_file && awk -v speciestaxid=\$i -F \"\\t\" '\$7==speciestaxid && \$12==\"Complete Genome\" && \$11==\"latest\"{print \$20}' $out_dir/$run_id/assembly_summary.txt >> $out_dir/$run_id/ftpdirpaths.txt; done; elif [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" == 'family' ]; then for i in \${ALL_SPECIES_THIS_FAMILY[@]}; do echo \"identifying species IDs corresponding to family ID \$FAMILYTAXID: \$i...\" >> $log_file && awk -v speciestaxid=\$i -F \"\\t\" '\$7==speciestaxid && \$12==\"Complete Genome\" && \$11==\"latest\"{print \$20}' $out_dir/$run_id/assembly_summary.txt >> $out_dir/$run_id/ftpdirpaths.txt; done; elif [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" == 'order' ]; then for i in \${ALL_SPECIES_THIS_ORDER[@]}; do echo \"identifying species IDs corresponding to order ID \$ORDERTAXID: \$i...\" >> $log_file && awk -v speciestaxid=\$i -F \"\\t\" '\$7==speciestaxid && \$12==\"Complete Genome\" && \$11==\"latest\"{print \$20}' $out_dir/$run_id/assembly_summary.txt >> $out_dir/$run_id/ftpdirpaths.txt; done; fi\n";
		  # NOTE: the overall intention of the following is to download each complete genome only once (however many iterations of this script are run) and to store a local copy of each genome's mash sketch, as $essentials_dir/sketches/$rootname.msh. The set of appropriate mash genomes - that is, those downloaded and/or those previously download and now locally stored and/or those made from local assemblies - are then combined using 'mash paste'.
		  print SH "while read filename; do name=\$(basename \"\$filename\") && rootname=\${name%_genomic.fna.gz} && if [ -e $essentials_dir/sketches/\$rootname.msh ]; then echo $essentials_dir/sketches/\$rootname.msh >> $out_dir/$run_id/mashsketchpaths.txt; else echo \$filename >> $out_dir/$run_id/ftpdirpaths2.txt; fi; done < $out_dir/$run_id/ftpdirpaths.txt\n"; # for each URL in ftpfilepaths, we first test whether a mash sketch already exists, i.e. whether this genome has previously been downloaded (and sketched). We do this by creating a file of sketch paths that will exist only if this is the case, mashsketchpaths.txt. Otherwise, we create a filtered subset of ftpfilepaths.txt that contains only those genomes not yet downloaded.
		  print SH "if [ -e $out_dir/$run_id/ftpdirpaths2.txt ]; then mv $out_dir/$run_id/ftpdirpaths2.txt $out_dir/$run_id/ftpdirpaths.txt; else truncate -s 0 $out_dir/$run_id/ftpdirpaths.txt; fi\n"; # if ftpdirpaths2.txt exists, this means that we have already downloaded either none of the genomes, or some of them - having pre-made mash sketches for the remainder. If so, we overwrite the original list of all genomes (that is, ftpdirpaths.txt) to contain only the subset of URLs we don't have (if we haven't previously downloaded any of them, all we're doing here is overwriting a list of URLs with an identical list). Alternatively, if we haven't created ftpdirpaths2.txt above, this means we already have a mash sketch of EVERY existing genome. In this case, we empty ftpdirpaths so we don't downloaded them again.
		  print SH "awk 'BEGIN{FS=OFS=\"/\";filesuffix=\"genomic.fna.gz\"}{ftpdir=\$0;asm=\$10;file=asm\"_\"filesuffix;print ftpdir,file}' $out_dir/$run_id/ftpdirpaths.txt > $out_dir/$run_id/ftpfilepaths.txt\n";
		  print SH "rm $out_dir/$run_id/ftpdirpaths.txt\n";
		  print SH "wget -i $out_dir/$run_id/ftpfilepaths.txt --spider -nv -a $out_dir/$run_id/linktestlog.txt 2>&1\n"; # test whether each URL in ftpfilepaths.txt exists
		  print SH "cat $out_dir/$run_id/linktestlog.txt | awk '/fna.gz\" \\[1\\]/{print \$4}' > $out_dir/$run_id/confirmedftpfilepaths.txt\n";
		  print SH "rm $out_dir/$run_id/linktestlog.txt\n";
		  print SH "mv $out_dir/$run_id/confirmedftpfilepaths.txt $out_dir/$run_id/ftpfilepaths.txt\n";
		  print SH "mkdir $out_dir/$run_id/temp_dir_genome_storage\n";
		  print SH "wget -i $out_dir/$run_id/ftpfilepaths.txt -P $out_dir/$run_id/temp_dir_genome_storage &>> $log_file\n"; # unlike the other wget commands, which are tested before running, we have already affirmed that each URL in ftpfilepaths.txt exists (see above)
		  print SH "find $out_dir/$run_id/temp_dir_genome_storage -maxdepth 1 -type f > $out_dir/$run_id/fa_filepaths_for_mash_sketch.txt\n";
		  print SH "while read filename; do name=\$(basename \"\$filename\") && rootname=\${name%_genomic.fna.gz} && $mash_path sketch -p $num_threads -o $essentials_dir/sketches/\$rootname \"\$filename\" &>> $log_file; done < $out_dir/$run_id/fa_filepaths_for_mash_sketch.txt\n"; # make one mash sketch per downloaded genome, if it does not already exist
		  print SH "while read filename; do name=\$(basename \"\$filename\") && rootname=\${name%_genomic.fna.gz} && echo $essentials_dir/sketches/\$rootname.msh >> $out_dir/$run_id/msh_filepaths_for_mash_paste.txt; done < $out_dir/$run_id/fa_filepaths_for_mash_sketch.txt\n"; # combine the mash sketches for each relevant genome
		  print SH "rm $out_dir/$run_id/fa_filepaths_for_mash_sketch.txt\n";
		  print SH "if [ -e $out_dir/$run_id/mashsketchpaths.txt ]; then cat $out_dir/$run_id/mashsketchpaths.txt >> $out_dir/$run_id/msh_filepaths_for_mash_paste.txt; fi\n";
		  print SH "$mash_path paste $out_dir/$run_id/$run_id.fa -l $out_dir/$run_id/msh_filepaths_for_mash_paste.txt &>> $log_file\n";
		  print SH "number_of_complete_genomes=\"`wc -l < $out_dir/$run_id/msh_filepaths_for_mash_paste.txt`\"\n";
		  print SH "rm $out_dir/$run_id/msh_filepaths_for_mash_paste.txt\n";
		  print SH "if [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" == 'species' ]; then echo \"--> assembly_summary.txt parsed to identify \$number_of_complete_genomes complete genomes for taxon ID \$SPECIESTAXID\" >> $log_file; elif [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" == 'genus' ]; then echo \"--> assembly_summary.txt parsed to identify \$number_of_complete_genomes complete genomes for taxon ID \$GENUSTAXID\" >> $log_file; elif [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" == 'family' ]; then echo \"--> assembly_summary.txt parsed to identify \$number_of_complete_genomes complete genomes for taxon ID \$FAMILYTAXID\" >> $log_file; elif [ \"\$taxonomic_rank_in_which_we_have_highest_confidence\" == 'order' ]; then echo \"--> assembly_summary.txt parsed to identify \$number_of_complete_genomes complete genomes for taxon ID \$ORDERTAXID\" >> $log_file; fi\n";
		  print SH "echo \"NUMBER OF COMPLETE GENOMES ASSIGNED THIS TAXON ID, OBTAINED AFTER PARSING ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt (DOWNLOADED \$date_last_modified): \$number_of_complete_genomes\" >> $ref_genome_download_log\n";
		  print SH "echo \"URLS OF DOWNLOADED GENOMES:\" >> $ref_genome_download_log\n";
		  print SH "cat $out_dir/$run_id/ftpfilepaths.txt >> $ref_genome_download_log\n";
		  print SH "rm $out_dir/$run_id/ftpfilepaths.txt\n";
		  print SH "echo \"PATHS TO PRE-DOWNLOADED GENOMES:\" >> $ref_genome_download_log\n";
		  print SH "if [ -e $out_dir/$run_id/mashsketchpaths.txt ]; then cat $out_dir/$run_id/mashsketchpaths.txt >> $ref_genome_download_log; fi\n";
		  print SH "if [ -e $out_dir/$run_id/mashsketchpaths.txt ]; then rm $out_dir/$run_id/mashsketchpaths.txt; fi\n";
		  if ($fq_se ne '')
			{ print SH "cat $fq_1 $fq_2 $fq_se > $out_dir/$run_id/$run_id.concatenated_reads.fq.gz\n"; }
		  elsif ($fq_se eq '')
			{ print SH "cat $fq_1 $fq_2 > $out_dir/$run_id/$run_id.concatenated_reads.fq.gz\n"; }
		  print SH "$mash_path sketch -p $num_threads -o $out_dir/$run_id/$run_id.fq $out_dir/$run_id/$run_id.concatenated_reads.fq.gz &>> $log_file\n";
		  print SH "rm $out_dir/$run_id/$run_id.concatenated_reads.fq.gz\n";
		  print SH "$mash_path dist -p $num_threads $out_dir/$run_id/$run_id.fq.msh $out_dir/$run_id/$run_id.fa.msh > $out_dir/$run_id/$run_id.mash-dist_output.tsv\n";
		  print SH "closest_mash_hit=\$(basename \$(sort -g -k3 $out_dir/$run_id/$run_id.mash-dist_output.tsv | head -1 | cut -f2))\n";
		  print SH "closest_mash_hit_root=\${closest_mash_hit%_genomic.fna.gz}\n";
		  print SH "path_to_closest_mash_hit=\$(grep \$closest_mash_hit_root $out_dir/$run_id/assembly_summary.txt | head -1 | awk '{print \$(NF)}')\n";
		  print SH "FTPPATHG=\$path_to_closest_mash_hit/\$closest_mash_hit_root'_genomic.fna.gz'\n";
		  print SH "FTPPATHGFF=\$path_to_closest_mash_hit/\$closest_mash_hit_root'_genomic.gff.gz'\n";
		  print SH "FTPPATHFAA=\$path_to_closest_mash_hit/\$closest_mash_hit_root'_protein.faa.gz'\n";
		  print SH "if curl --head --fail --silent \$FTPPATHG >/dev/null; then wget \$FTPPATHG -O $out_dir/$run_id/ref.fa.gz &>> $log_file; else echo \"ERROR: unable to download the reference genome at \$FTPPATHG\n\" >> $log_file && exit 1; fi\n";
		  print SH "gunzip $out_dir/$run_id/ref.fa.gz\n";
		  print SH "if curl --head --fail --silent \$FTPPATHGFF >/dev/null; then wget \$FTPPATHGFF -O $out_dir/$run_id/ref.gff.gz &>> $log_file; else echo \"ERROR: unable to download the GFF for the reference genome at \$FTPPATHGFF\n\" >> $log_file && exit 1; fi\n";
		  print SH "gunzip $out_dir/$run_id/ref.gff.gz\n";
		  print SH "if curl --head --fail --silent \$FTPPATHFAA >/dev/null; then wget \$FTPPATHFAA -O $out_dir/$run_id/ref.faa.gz &>> $log_file; else echo \"ERROR: unable to download the FAA for the reference genome at \$FTPPATHFAA\n\" >> $log_file && exit 1; fi\n";
		  print SH "gunzip $out_dir/$run_id/ref.faa.gz\n";
		  print SH "echo \"\nURL FOR CLOSEST MASH HIT (SELECTED AS REFERENCE GENOME): \$FTPPATHG\" >> $ref_genome_download_log\n";
		  # the URL for the closest match hit will be, e.g., "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/925/GCF_000016925.1_ASM1692v1/GCF_000016925.1_ASM1692v1_genomic.fna.gz"
		  # this name (stored above both as $closest_mash_hit_root and, in the "while read filename" loops, as 'rootname') embeds within it unique identifying features. For instance, "GCF_000016925.1_ASM1692v1" comprises the "assembly accession" (GCF_000016925.1) and "assembly name" (ASM1692v1); the former is unique but the latter is not. The assembly accession is "a unique identifier for the set of sequences in this particular version of the genome assembly" and the assembly (asm) name is "the submitter's name for the genome assembly, when one was provided, otherwise a default name, in the form ASM#####v#, provided by NCBI". Both are defined at ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt. That this name is unique means it is appropriate for storing the local mash sketches. Genomes that may later become deprecated (anticipated with periodic updates to assembly_summary.txt) will result in the local storage of 'dead' mash sketches, that may never be referred to again. This is not a concern: the sketches dir may be periodically emptied as it is easily remade.
		  print SH "organism_name=\$(grep \$closest_mash_hit_root $out_dir/$run_id/assembly_summary.txt | awk -F '\\t' '{print \$8}')\n"; # "Organism name: the scientific name of the organism from which the sequences in the genome assembly were derived. This name is taken from the NCBI Taxonomy record for the taxid specified in column 6. Some older taxids were assigned at the strain level and for these the organism name will include the strain. Current practice is only to assign taxids at the species level; for these the organism name will be just the species, however, the strain name will be reported in the infraspecific_name field (column 9)" (ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt)
		  print SH "infraspecific_name=\$(grep \$closest_mash_hit_root $out_dir/$run_id/assembly_summary.txt | awk -F '\\t' '{print \$9}')\n"; # "Infraspecific name: the strain, breed, cultivar or ecotype of the organism from which the sequences in the genome assembly were derived. Data are reported in the form tag=value, e.g. strain=AF16. Strain, breed, cultivar and ecotype are not expected to be used together, however, if they are then they will be reported in a list separated by ", /". Empty if no strain, breed, cultivar or ecotype is specified on the genomic sequence records" (ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt)
		  print SH "echo \"ORGANISM NAME: \$organism_name\" >> $ref_genome_download_log\n";
		  print SH "echo \"INTRASPECIFIC NAME: \$infraspecific_name\" >> $ref_genome_download_log\n";
		  print SH "mv $out_dir/$run_id/$run_id.mash-dist_output.tsv $out_dir/$run_id/logs/$log_number"."_mash_genome_prediction/$run_id.mash-dist_output.tsv\n";
		  print SH "rm $out_dir/$run_id/assembly_summary.txt\n";
		  print SH "rm $out_dir/$run_id/$run_id.fq.msh $out_dir/$run_id/$run_id.fa.msh\n";
		  print SH "rm -r $out_dir/$run_id/temp_dir_genome_storage\n";		  
		  $log_number++;
=cut	  
		  ### TO BE COMPLETED ###
		  # (13) IF THE TOP GENUS HIT IS PREDICTED TO BE MYCOBACTERIUM, THEN RUN MYKROBE TO PREDICT THE SPECIFIC STRAIN (WHICH WILL BE CHOSEN AS A REFERENCE GENOME). FOR THIS PURPOSE, MYKROBE OVERRULES THE OUTPUT OF KAIJU/MASH.
		  ## IMPORTANT: the mykrobe predictor codebase is scheduled for deprecation in December 2018; see https://github.com/iqbal-lab/Mykrobe-predictor
		  if ($run_mykrobe_if_tb eq 'yes')
			{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
			  print SH "elapsed_time=\$SECONDS\n";
			  print SH "[ \"\$top_genus_hit_pe\" == 'Mycobacterium' ]; then echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. running Mykrobe on $run_id...\n\" >> $log_file; fi\n";
			  print SH "[ \"\$top_genus_hit_pe\" == 'Mycobacterium' ]; then $mykrobe_path predict $run_id tb --tmp $out_dir/$run_id/$out_dir/$run_id/temp_mykrobe --skeleton_dir $out_dir/$run_id/$out_dir/$run_id/skeletondir_mykrobe --mccortex31_path $mccortex_path -t $num_threads --seq $fq_1 --seq $fq_2; fi\n";
			}
=cut
		}
	  
	  # (14) INDEX THE REFERENCE GENOME
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. indexing reference genome...\n\" >> $log_file\n";
	  print SH "$bwa_path index $ref 2>> $log_file\n";
	  
	  # (15) ALIGN CLEANED, HUMAN-DEPLETED, READS TO REFERENCE GENOME USING BWA. THE BAM IS POST-PROCESSED WITH PICARD TO COORDINATE-SORT, FIX MATE INFORMATION (IF RELEVANT) AND TO MARK DUPLICATES.
	  # (15i) alignment
	  for(my $x=0;$x<=1;$x++) # this loop aligns, in turn, the set of PE and SE reads (those which, after the above trimming/QC, have lost their mate)
		{ my $pe_or_se = '';
		  if ($x == 0) { $pe_or_se = 'paired-end'; } elsif ($x == 1) { $pe_or_se = 'single-end'; }
		  next if (($pe_or_se eq 'single-end') and ($run_trimmomatic eq 'no'));
		  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. aligning the quality- and adapter-trimmed $pe_or_se reads from $run_id...\n\" >> $log_file\n";
		  if ($pe_or_se eq 'paired-end')
			{ print SH "$bwa_path mem -R '\@RG\\tID:group_$run_id\\tSM:sample_$run_id\\tPL:Illumina\\tLIB:lib_$run_id\\tPU:unit_$run_id' -t $num_threads -M $ref $fq_1 $fq_2 2>> $log_file | $sambamba_path view -S -f bam -t $num_threads -o $out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam /dev/stdin\n"; # IF JUST USING SAMTOOLS RATHER THAN SAMBAMBA: print SH "$bwa_path mem -R '\@RG\\tID:group_$run_id\\tSM:sample_$run_id\\tPL:Illumina\\tLIB:lib_$run_id\\tPU:unit_$run_id' -t $num_threads -M $ref $fq_1 $fq_2 2>> $log_file | $samtools_path view -Shb - > $out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam\n";
			}
		  elsif ($pe_or_se eq 'single-end')
			{ print SH "$bwa_path mem -R '\@RG\\tID:group_$run_id\\tSM:sample_$run_id\\tPL:Illumina\\tLIB:lib_$run_id\\tPU:unit_$run_id' -t $num_threads -M $ref $fq_se 2>> $log_file | $sambamba_path view -S -f bam -t $num_threads -o $out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam /dev/stdin\n"; # IF JUST USING SAMTOOLS RATHER THAN SAMBAMBA: print SH "$bwa_path mem -R '\@RG\\tID:group_$run_id\\tSM:sample_$run_id\\tPL:Illumina\\tLIB:lib_$run_id\\tPU:unit_$run_id' -t $num_threads -M $ref $fq_se 2>> $log_file | $samtools_path view -Shb - > $out_dir/$run_id/$run_id.unsorted.$pe_or_se.bam\n";
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
	  # (15ii) merge the PE and SE BAMs
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
	  # (15iii) obtain summary statistics for the final, merged, BAM
	  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_bam_statistics\n";
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. obtaining summary statistics for the merged BAM of $run_id...\n\" >> $log_file\n";
	  print SH "$samtools_path stats $out_dir/$run_id/$run_id.bam > $out_dir/$run_id/logs/$log_number"."_bam_statistics/$run_id.samtools_stats_output.txt\n";
	  print SH "$samtools_path flagstat $out_dir/$run_id/$run_id.bam > $out_dir/$run_id/logs/$log_number"."_bam_statistics/$run_id.samtools_flagstat_output.txt\n";
	  $log_number++;
	  
	  # (16) CALL VARIANTS
	  # (16i) variant calling using mpileup
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. calling variants from $run_id using mpileup...\n\" >> $log_file\n";
	  if ($output_all_sites eq 'no') # if no, we are only interested in calling the variant sites
		{ print SH "$bcftools_path mpileup -Ou -f $ref $out_dir/$run_id/$run_id.bam 2>> $log_file | $bcftools_path call --threads $num_threads --ploidy 1 -mv -Ov -o $out_dir/$run_id/$run_id.vcf\n"; } # the "v" in "-mv" specifies that the output file only contain variant sites; omit this to output all sites
	  elsif ($output_all_sites eq 'yes')
		{ print SH "$bcftools_path mpileup -Ou -f $ref $out_dir/$run_id/$run_id.bam 2>> $log_file | $bcftools_path call --threads $num_threads --ploidy 1 -m -Ov -o $out_dir/$run_id/$run_id.vcf\n"; }
      # (16ii) variant calling using cortex (which first requires that the reference genome be indexed with Stampy)
#	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
#	  print SH "elapsed_time=\$SECONDS\n";
#	  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. calling variants from $run_id using cortex...\n\" >> $log_file\n";
#	  print SH "$stampy_path -G ref $ref 2>> $log_file\n";
#	  print SH "$stampy_path -g ref -H ref 2>> $log_file\n";
#	  print SH "echo \"$ref\" > $out_dir/$run_id/cortex.in.index_ref.fofn\n";
#	  print SH "echo -e \"sample\t$out_dir/$run_id/cortex.in.index_ref.fofn\t.\t.\" > $out_dir/$run_id/cortex.in.index\n";
	  
	  # (17) CLOCKWORK: USE MINOS TO ADJUDICATE BETWEEN VCFs PRODUCED BY BWA/MPILEUP AND BWA/CORTEX, SO CREATING ONE FINAL VCF
	  ## to be added ## 
	  
	  # (18) REGULARISE VCF SO THAT DIFFERENT REPRESENTATIONS OF THE SAME INDEL OR COMPLEX VARIANT ARE NOT COUNTED AS DIFFERENT VARIANTS (THIS IS USEFUL IF LATER COMPARING VCFs PRODUCED BY DIFFERENT TOOLS)
	  if ($regularise_vcf eq 'yes')
		{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
		  print SH "elapsed_time=\$SECONDS\n";
		  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. regularising VCF for $run_id...\n\" >> $log_file\n";
		  print SH "$vcfallelicprimitives_path $out_dir/$run_id/$run_id.vcf > $out_dir/$run_id/$run_id.regularised.vcf 2>> $log_file\n";
		  print SH "rm --interactive=never $out_dir/$run_id/$run_id.vcf\n";
		  print SH "mv $out_dir/$run_id/$run_id.regularised.vcf $out_dir/$run_id/$run_id.vcf\n";
		}
	  
	  # (19) REMOVE LOW-QUALITY RECORDS FROM THE VCF
	  if ($filter_vcf eq 'yes')
		{ print SH "$vcffilter_path -f \"DP4 > $min_DP4 & MQ > $min_MQ & AN = $max_number_of_alleles_at_site\" $out_dir/$run_id/$run_id.vcf >> $out_dir/$run_id/$run_id.filtered.vcf 2>> $log_file\n";
		  if ($keep_original_vcf eq 'no')
			{ print SH "mv $out_dir/$run_id/$run_id.filtered.vcf $out_dir/$run_id/$run_id.vcf\n"; }
		}
	  
	  # (20) SORT, COMPRESS AND INDEX THE FINAL VCF(s)
	  print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
	  print SH "elapsed_time=\$SECONDS\n";
	  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. sorting, compressing and indexing VCF for $run_id...\n\" >> $log_file\n";
	  print SH "$vcfsort_path $out_dir/$run_id/$run_id.vcf > $out_dir/$run_id/$run_id.sorted.vcf 2>> $log_file\n";
	  print SH "mv $out_dir/$run_id/$run_id.sorted.vcf $out_dir/$run_id/$run_id.vcf\n";
	  print SH "$bgzip_path $out_dir/$run_id/$run_id.vcf\n";
	  print SH "$tabix_path -p vcf $out_dir/$run_id/$run_id.vcf.gz\n";
	  if (($filter_vcf eq 'yes') and ($keep_original_vcf eq 'yes'))
		{ print SH "$vcfsort_path $out_dir/$run_id/$run_id.filtered.vcf > $out_dir/$run_id/$run_id.sorted.filtered.vcf 2>> $log_file\n";
		  print SH "mv $out_dir/$run_id/$run_id.sorted.filtered.vcf $out_dir/$run_id/$run_id.filtered.vcf\n";
		  print SH "$bgzip_path $out_dir/$run_id/$run_id.filtered.vcf\n";
		  print SH "$tabix_path -p vcf $out_dir/$run_id/$run_id.filtered.vcf.gz\n";
		}
	  
	  # (21) ASSEMBLE THE SET OF CLEANED, CLASSIFIED READS, QC THE RESULTING SET OF SCAFFOLDS, AND ANNOTATE
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
		  # (21i) QC the SPAdes assembly using Quast
		  if ($qc_assembly eq 'yes')
			{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
			  print SH "elapsed_time=\$SECONDS\n";
			  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. creating summary statistics for the de novo assembly of $run_id...\n\" >> $log_file\n";
			  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_quast_report_on_de_novo_assembly\n";
			  if ($fq_se ne '')
				{ print SH "python $quast_path -r $ref -g $gff -t $num_threads --circos --gene-finding --rna-finding --conserved-genes-finding --pe1 $fq_1 --pe2 $fq_2 --single $fq_se -o $out_dir/$run_id/logs/$log_number"."_quast_report_on_de_novo_assembly $out_dir/$run_id/$run_id.fa &>> $log_file\n"; }
			  elsif ($fq_se eq '')
				{ print SH "python $quast_path -r $ref -g $gff -t $num_threads --circos --gene-finding --rna-finding --conserved-genes-finding --pe1 $fq_1 --pe2 $fq_2 -o $out_dir/$run_id/logs/$log_number"."_quast_report_on_de_novo_assembly $out_dir/$run_id/$run_id.fa &>> $log_file\n"; }
			  $log_number++;
			}
		  # (21ii) annotate the SPAdes assembly using Prokka. FYI, a tutorial for Prokka and Roary ('how to determine a pangenome from a collection of isolate genome sequences') is available here: https://github.com/microgenomics/tutorials/blob/master/pangenome.md
		  if ($annotate_assembly eq 'yes')
			{ print SH "current_date_time=\"`date \"+%Y-%m-%d %H:%M:%S\"`\"\n";
			  print SH "elapsed_time=\$SECONDS\n";
			  print SH "echo \"-->\$current_date_time. elapsed time: \$((\$elapsed_time / 60)) minutes and \$((\$elapsed_time % 60)) seconds. annotating the de novo assembly of $run_id...\n\" >> $log_file\n";
			  my $prokka_outdir = "$out_dir/$run_id/logs/$log_number"."_prokka_annotation_of_de_novo_assembly";
			  my $species_name_for_prokka = 'species'; my $strain_name_for_prokka = 'strain';
			  print SH "species_name_for_prokka=\$( echo \$top_species_hit_pe | perl -pe 's/(.*?) (.*?)/\$2/' )\n";
			  print SH "strain_name_for_prokka=\$( echo \$infraspecific_name | perl -pe 's/strain\\=//' )\n";
			  print SH "$prokka_path --cpus $num_threads --kingdom Bacteria --outdir $prokka_outdir --prefix $run_id --locustag $run_id --compliant --genus \"\$top_genus_hit_pe\" --species \"\$species_name_for_prokka\" --strain \"\$strain_name_for_prokka\" --proteins $faa --usegenus --rfam $out_dir/$run_id/$run_id.fa &>> $log_file\n";
			  $log_number++;
			  ## TO BE ADDED: the Prokka parameter --gram will search for signal peptides, provided you specify what Gram status you're looking for.
			  ## TO BE ADDED: the Prokka parameter --genus will use genus-specific databases, if they are available. If one is NOT available, in $progs/prokka/db, it is possible to make it: see "adding a genus database" in https://github.com/tseemann/prokka
			}
		  ### TO BE ADDED: create a mash sketch of the scaffolds (so that, if taxonomically appropriate and if it passes QC, this assembly can be added to a set of public assemblies downloaded for mash clustering; essentially, this feeds data from this pipeline back into itself, for iterative improvement)
		}
	  
	  # (22) RUN TETYPER TO SEARCH FOR THE blaKPC TRANSPOSON Tn4401
	  if (($run_tetyper eq 'yes') and ($create_assembly eq 'yes'))
		{ # IMPORTANT: see https://github.com/aesheppard/TETyper, which states that "TETyper was designed with the blaKPC transposon Tn4401 in mind. A Tn4401b reference sequence is provided with TETyper (Tn4401b-1.fasta), as well as example profile definitions for SNVs / deletions with respect to this reference (struct_profiles.txt and snp_profiles.txt)." These files are used on the command line here but may be varied (see REQUIREMENTS at the top of this script)
		  # IMPORTANT: note that TETyper expects blastn on the $PATH but that more recent versions of the BLAST+ toolkit have deprecated this command
		  print SH "mkdir -p $out_dir/$run_id/logs/$log_number"."_tetyper_output\n";
		  my $tetyper_output_prefix = "$out_dir/$run_id/logs/$log_number"."_tetyper_output/$run_id";
		  print SH "$tetyper_path --ref $tetyper_ref --fq1 $fq_1 --fq2 $fq_2 --assembly $out_dir/$run_id/$run_id.fa --threads $num_threads --flank_len $tetyper_flank_length --struct_profiles $tetyper_struct_profiles --snp_profiles $tetyper_snp_profiles --outprefix $tetyper_output_prefix &>> $log_file\n";
		  $log_number++;
		}
	  
	  # (23) IMPLEMENT PHYLOGENOMIC ANALYSES (USING ELEPHANT WALK?)
	  ## to be added ##
	  
	  # (24) TIDY UP
	  # delete the now redundant reads
	  if ($keep_cleaned_fq eq 'no')
		{ if 	($fq_se ne '') { print SH "rm $fq_1 $fq_2 $fq_se\n"; }
		  elsif ($fq_se eq '') { print SH "rm $fq_1 $fq_2\n"; 		 }
		}
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
	  if (($filter_vcf eq 'yes') and ($keep_original_vcf eq 'yes'))
		{ print SH "md5sum $out_dir/$run_id/$run_id.filtered.vcf.gz > $out_dir/$run_id/$run_id.filtered.vcf.gz.md5sum\n";
		  print SH "md5sum $out_dir/$run_id/$run_id.filtered.vcf.gz.tbi > $out_dir/$run_id/$run_id.filtered.vcf.gz.tbi.md5sum\n";
		}
	  
	  # (25) FINALLY, HOW LONG DID IT TAKE TO COMPLETE THIS PIPELINE?
	  print SH "duration=\$SECONDS\n";
	  print SH "echo \"--> total time taken to process $run_id: \$((\$duration / 60)) minutes and \$((\$duration % 60)) seconds\n\" >> $log_file\n";
	  print SH "exit 0\n";
	  close(SH) or die $!;
	  close(LOG) or die $!;
	  if ($where eq 'eddie') { print QSUB "qsub $sh_file\n"; }
	}
if ($where eq 'eddie') { close(QSUB) or die $!; }
exit 1;