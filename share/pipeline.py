import kDNA_annotation as ka
from os import system

config_file = 'config.yaml'

ka.clean_mini_and_maxicircles()
# ka.mRNA_process()
# system('align_maxi')
# system('align_mini')
# ka.hq_gRNAs()
# system(f'meme Work_files/hq_gRNAs.fasta -dna -oc Work_files/Meme -mod zoops -nmotifs 3 -minw 5 -maxw 30 -objfun classic -markov_order 0 -p 4')
# ka.extract_motifs( )
# ka.mO_scoring()
# ka.identify_cassettes()
# ka.identify_CSB_gRNAs()

########### Only run if transcriptomics is available ############
# ka.transcripts()
# ka.find_transcript_p()
# ka.predict_transcript_end_pos()
# ka.predict_expression()
########### Only run if transcriptomics is available ############

# ka.annotate_minicircles()
