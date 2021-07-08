#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
##                               SDgenerator
##  Generación de muestras sintéticas (en formato FASTQ) compuestas por una 
##  población definida de DVGs (en tipo y proporciones)
## 
###############################################################################
## Author: Maria Jose Olmo-Uceda
## Version: 1.0
## Email: maolu@alumni.uv.es
## Date: 2021/04/21
###############################################################################


"""
Objetivo: generar un dataset simulado que contenga el conjunto de genomas 
        compuesto por el genoma wt y los dvgs especificados en las proporciones
        indicadas en el input.csv
        Las proporciones (en tanto por 1) tienen en cuenta la longitud total, 
        de manera que si un SD tiene un único evento al 0.5 significará que hay
        un genoma completo DVG por cada genoma completo wt.
Argumentos de entrada:
    -f: csv con cabecera BP,RI,DVG_type,proportion
    -r: genoma de referencia (wt) en formato fasta
    -N: número aproximado total de reads con el que queremos generar el dataset
    -l: longitud de las lecturas. [Default is 100 nt].
    -o: nombre que le daremos al fq generado
    -k: mantener los fastqs individuales generados. Options: y/n 
        [Default is n]
    #-c: genoma contaminante (fasta) si queremos que el fq esté contaminado
Output:
    SD[_name].fq


Usage python3 SDgenerator.py -f input.csv -r reference.fasta -N total_reads
                        [-l length_reads] [-o name_SD] [-k]
"""

import os
import argparse
import time
import datetime

# Third party imports 
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


import functionsSDgenerator as my

def main():
    t0 = time.time()

    ## Lectura de los parámetros introducidos por el usuario
    input_csv, reference, N, len_reads, name, outer_distance, keep_individuals \
     = my.read_arguments()

    ### name outputs
    #### csv with necessary info
    composition_table = name + '_composition.csv'
    sd_name = name + '.fq'

    ## Extract some common variables
    ### reading reference sequence
    for seq_record in SeqIO.parse(reference, "fasta"):
        reference_sequence = seq_record.seq
    ### length of wt
    len_wt = len(reference_sequence)

    #--------------------------------------------------------------------------
    ## Generación de los fasta de cada dvg
    #--------------------------------------------------------------------------
    t1 = time.time()
    # Calculo del número de reads de cada especie
    ## Put csv info into a pd.DataFrame and add the length dvg column
    df = my.add_lengths_to_df(input_csv, len_wt)

    # proportion of the wt genome
    proportion_wt = my.prop_wt(df)

    ## Calculate the N_dvg and write the information estracted in a new csv
    df = my.add_N_dvgs(df, N, len_reads, len_wt, proportion_wt)
    

    ## Generate de basename of each event to facilitate the posterior naming
    ## of the individual files
    df = my.add_basename_dvg_files(df)

    # Save the anotated table
    df.to_csv(composition_table, index=0)

    ## In case doesn't exist, create the directory Outputs/fastas
    if 'Outputs' not in os.listdir():
        os.system('mkdir Outputs')
    
    if 'fastas' not in os.listdir('Outputs'):
        os.system('mkdir Outputs/fastas')
    
    ## Write the fastas with the sequences of each event
    my.write_all_dvg_sequences(reference_sequence, df)

    #--------------------------------------------------------------------------
    ## Simulación de los fastqs de cada dvg
    #--------------------------------------------------------------------------
    t2 = time.time()
    ## In case doesn't exist, create directory Outputs/fastqs
    if 'fastqs' not in os.listdir('Outputs'):
        os.system('mkdir Outputs/fastqs')
    ## In case doesn't exist, create the specific fastq directory
    dir_fastqs = "Outputs/fastqs/" + name + "/"
    if 'dir_fastqs' not in os.listdir("Outputs/fastqs/"):
        os.system('mkdir {}'.format(dir_fastqs))
    else:
        print("The directory {} has been rewritten.".format(dir_fastqs))
    
    ## Generate de individual fastq 
    df.apply(lambda row : my.generate_simulations(name, df, len_reads, 
                                                    outer_distance))

    ## Generate wt simulation
    my.simulation_wt(name, df, reference, N, len_reads, len_wt, proportion_wt, 
                    outer_distance)
    
    #--------------------------------------------------------------------------
    ## Generación de un único fichero fq con las reads de todos los genomas
    #--------------------------------------------------------------------------
    t3 = time.time()
    #my.concatenate_fastq(df, name)
    my.concatenate_fastq(df, N, len_reads, len_wt, proportion_wt, name,
                     as_single_end=True)

    if not keep_individuals:
        os.system("rm -r {}".format(dir_fastqs))
    
    tf = time.time()
    print("-"*79)
    print("The synthetic dataset {} has been generated in the directory '{}'.".\
        format(name + '.fq', "Outputs/"))
    print("Total time: {} s".format (tf-t0))
    print("-"*79)


if __name__=='__main__':
    main()