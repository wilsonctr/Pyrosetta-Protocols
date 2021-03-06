from Bio import SeqIO
import pandas as pd
import os
import subprocess
from urllib import urlopen
from Bio import AlignIO
from Bio import pairwise2
import urllib, urllib2
import gzip
import requests

rosetta_path = '/share/siegellab/software/Rosetta_group_0618'

def phmmer_local_old(filename):

    os.system('phmmer -A output_alignment %s /share/siegellab/software/pdb_seqres.txt'%filename)
    os.system('esl-reformat fasta output_alignment > hmmer_results.fasta')
    for i in SeqIO.parse('hmmer_results.fasta', 'fasta'):
        os.system('esl-sfetch /share/siegellab/software/pdb_seqres.txt "%s" >> full_fasta.fasta'%i.id[0:6])

def phmmer_local(filename):

    os.system('phmmer -A output_alignment %s /share/siegellab/software/pdb_seqres.txt'%filename)
    os.system('esl-reformat fasta output_alignment > hmmer_results.fasta')


def get_pdb_seq(filename):

    with open('full_fasta.fasta', 'w') as f:
        for i in SeqIO.parse(filename, 'fasta'):
            data = urllib2.urlopen('https://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=%s&compressionType=uncompressed'%i.id[0:4])
            f.write(data.read())
            f.write('\n')
            
        
def get_pdb(pdb_code):
    
    #Get pdb file, substitute pdb code
    pdb_file_url = "https://files.rcsb.org/view/%s.pdb"%pdb_code[0:4]    
    pdb_file_html = urlopen(pdb_file_url).read()  
    with open('%s_A.pdb'%pdb_code[0:4].upper(), 'w') as f:
        
        for line in pdb_file_html.split('\n'):
            if line[0:4] == "ATOM":
                f.write(line)
                f.write('\n')
    
    return pdb_code

def remove_duplucate_seq(fasta_file, query_seq):

    names, sequences = [], []
    for i in SeqIO.parse(fasta_file, 'fasta'): 
        names.append(i.id)
        sequences.append(str(i.seq))
    
    df = pd.DataFrame()
    df['names'] = names
    df['sequences'] = sequences

    query_seq_length = len(query_seq)
    df = df[df['sequences'].map(len) < query_seq_length * 1.2]
    df = df[df['sequences'].map(len) > query_seq_length * 0.5]
    
    df.drop_duplicates(subset=['sequences'], inplace=True)
    
    with open('full_fasta_dupliecates_dropped.fasta', 'w') as f:
        for i,j in zip(df['names'], df['sequences']):
            f.write('>%s\n'%i[0:6].replace(':','_'))
            f.write('%s\n'%j)

def write_slurm_submittion_file(rosetta_path, seq_id):            

    with open('rosetta_cm_%s/submit.sh'%seq_id, 'w') as f:
        f.write('''#!/bin/bash                                                                                                                    
#SBATCH --output=log.txt                                                                                                       
#SBATCH --array=1-100                                                                                                           
#SBATCH --mem-per-cpu=2G                                                                                                       
#SBATCH --time 1440                                                                                                         
#SBATCH --partition production

echo %s/main/source/bin/rosetta_scripts.linuxgccrelease @%s/rosetta_cm_%s/flags -suffix $RANDOM
%s/main/source/bin/rosetta_scripts.linuxgccrelease @%s/rosetta_cm_%s/flags -out:path results_%s -suffix $RANDOM'''%(rosetta_path, os.getcwd(), seq_id, rosetta_path, os.getcwd(), seq_id, seq_id))
        
def rewrite_xml_file(rosetta_cm_folder_path):
    
    with open('rosetta_cm_%s/rosetta_cm.xml'%seq.id, 'r') as f:
        old_xml = f.read()

        with open('rosetta_cm_%s/new_cm.xml'%seq.id, 'w') as f1:
            f1.write(old_xml.split('</Hybridize>')[0].replace('rosetta_cm', 'rosetta_cm_%s'%seq.id))
            f1.write('</Hybridize>\n')
            f1.write('''<FastRelax name="relax" scorefxn="fullatom"></FastRelax>''')
            f1.write(old_xml.split('</Hybridize>')[1].split('</PROTOCOLS>')[0])
            f1.write('<Add mover="relax"/>\n')
            f1.write('</PROTOCOLS>')
            f1.write(old_xml.split('</Hybridize>')[1].split('</PROTOCOLS>')[1])

def rewrite_flags(query_name):
    with open('rosetta_cm_%s/flags'%query_name, 'r') as f:
        flags = f.read().replace('-nstruct 20', '-nstruct 1').replace('rosetta_cm.xml','%s/rosetta_cm_%s/new_cm.xml'%(os.getcwd(), query_name))
        with open('rosetta_cm_%s/flags'%query_name, 'w') as f1:
            f1.write(flags)


for seq in SeqIO.parse('round_1_expressed.fasta', 'fasta'):
    
    #Hmmer search, download full length fasta
    os.system('echo ">%s\n%s\n" > %s_query.fasta'%(seq.id, str(seq.seq), seq.id))
    phmmer_local('query.fasta')
    get_pdb_seq('hmmer_results.fasta')
    print 'got pdb seq'
    remove_duplucate_seq('full_fasta.fasta', str(seq.seq))
    print 'dropped duplicates.'
    os.system('echo ">%s\n%s\n" >> full_fasta.fasta'%(seq.id, str(seq.seq)))
    #Make alignment
    os.system('/share/siegellab/software/Muscle_alignment -in full_fasta_dupliecates_dropped.fasta -out %s_aligned.fasta'%(seq.id))

    #Get pdb files
    pdb_list = []
    for i in SeqIO.parse('hmmer_results.fasta', 'fasta'):
        pdb_list.append(i.id[0:4].upper() + '_A')
    
    command = ('%s/tools/protein_tools/scripts/setup_RosettaCM.py '
               '--fasta %s_query.fasta ' 
               '--alignment %s_aligned.fasta ' 
               '--alignment_format fasta ' 
               '--templates %s.pdb '
               '--rosetta_bin %s/main/source/bin/ ')

    print 'Goint to start threading...'
    os.system(command%(rosetta_path, seq.id, seq.id, '.pdb '.join(pdb_list), rosetta_path))
    os.system('mv -f rosetta_cm rosetta_cm_%s'%seq.id)
    
    write_slurm_submittion_file(rosetta_path, seq.id)
    rewrite_xml_file('rosetta_cm_%s'%seq.id)
    rewrite_flags(seq.id)
    os.system('mkdir results_%s'%seq.id)
    os.system('sbatch rosetta_cm_%s/submit.sh'%seq.id)

    os.system('rm -rf rosetta_cm/')
    os.system('rm hmmer_results.fasta')
    os.system('rm full_fasta.fasta')
    os.system('rm full_fasta_dupliecates_dropped.fasta')
    os.system('rm %s_aligned.fasta'%seq.id)
    
