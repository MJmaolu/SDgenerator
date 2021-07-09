#!/usr/bin/env python
# -*- coding: utf-8 -*-

## functionsSDgenerator.py
###############################################################################
##                        Funciones para SDgenerator
##  Generación de muestras sintéticas (en formato FASTQ) compuestas por una 
##  población definida de DVGs (en tipo y proporciones)
## 
###############################################################################
## Author: Maria Jose Olmo-Uceda
## Version: 1.0
## Email: maolu@alumni.uv.es
## Date: 2021/04/21
###############################################################################

import os
import argparse
import time
import datetime

# Third party imports 
import pandas as pd
import numpy as np



#------------------------------------------------------------------------------
## Parsear argumentos de entrada
#------------------------------------------------------------------------------

def read_arguments():
    ## Lectura de los parámetros introducidos por el usuario
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--table_csv', 
                        help="Complete path to the csv file with the info:\
                            BP,RI,DVG_type,proportion. [Default 'input_table.csv']")
    parser.add_argument('-r', '--reference_genome', 
                        help="Complete path to the reference genome (wt) in\
                        fasta format. [Default SARS_CoV_2_Wuhan1]." )
    parser.add_argument('-N', '--number_reads', 
                        help="Number of reads (aproximate) wanted in the final dataset.\
                        [Default 1000000]")
    parser.add_argument('-l', '--len_reads', 
                        help="Length of the reads. [Default 100]")
    parser.add_argument('-o', '--name',
                        help="Name of the final dataset ({name}.fq). [Default 'date_SD.fq']")
    parser.add_argument('-d', '--outer_distance', 
                        help="Outer distance between the two ends. [Default 300]")
    parser.add_argument('-k', '--keep_individuals', 
                        help="""Keep the individuals fastqs. 
                        \nOptions: y/n [Default is n].""")
    parser.add_argument('-c', '--contaminant_genome',
                        help="Path to the contaminant genome in fasta format.\
                        Optional.")

    args = parser.parse_args()

    if args.table_csv:
        input_csv = args.table_csv
    else:
        input_csv = 'input_table.csv'
    
    if args.reference_genome:
        reference = args.reference_genome
    else:
        reference = 'SARSCoV2WuhanHu1.fasta'
    
    if args.number_reads:
        N = int(args.number_reads)
    else:
        N = 1000000
    
    if args.len_reads:
        len_reads = int(args.len_reads)
    else:
        len_reads = 100
    
    if args.name:
        name = args.name
    else:
        name = 'SD' #mejor fecha en la que se genera
    
    if args.outer_distance:
        outer_distance = int(args.outer_distance)
    else:
        outer_distance = 300

    if args.keep_individuals == 'y':
        keep_individuals = True
    else:
        keep_individuals = False
    
    if args.contaminant_genome:
        contaminant = args.contaminant_genome

    return input_csv, reference, N, len_reads, name, outer_distance, \
            keep_individuals


#------------------------------------------------------------------------------
## Cálculo de la longitud teórica de los DVGs y del número de reads a generar
#------------------------------------------------------------------------------

def first_seqment(len_wt, DVG_type, BP):
    """
    Longitud del primer segmento del DVG (hasta BP incluido) en función del tipo
    de DVG.
    Dos procedimientos diferentes para averiguar la longitud:
        a) Deletion_forward, Insertion_forward, 5cb/sb --> 1:BP
        b) Deletion_reverse, Insertion_reverse, 3cb/sb --> len_wt - (RI - 1)


    Args:
        len_wt   (int)   [in]    Longitud del genoma tomado como referencia
        DVG_type        (str)   [in]    Tipo de DVG
        BP              (int)   [in]    Break point: último nucleótido del 
                                        primer segmento
    Returns:
        len_first_seqment   (int)   [out]   Longitud del 1er segmento del DVG
    """

    a = ['Deletion_forward', 'Insertion_forward', '5cb/sb']
    b = ['Deletion_reverse', 'Insertion_reverse', '3cb/sb']

    if DVG_type in a:
        len_first_seqment = BP
    
    elif DVG_type in b:
        len_first_seqment = len_wt - BP + 1

    else:
        len_first_seqment = 'Error'
    
    return len_first_seqment


def second_seqment(len_wt, DVG_type, RI):
    """
    Longitud del segundo segmento del DVG (desde RI incluido) en función del tipo
    de DVG.
    Dos procedimientos diferentes para averiguar la longitud:
        a) Deletion_forward, Insertion_forward, 5cb/sb --> 1:BP
        b) Deletion_reverse, Insertion_reverse, 3cb/sb --> len_wt - (RI - 1)


    Args:
        len_wt   (int)   [in]    Longitud del genoma tomado como referencia
        DVG_type        (str)   [in]    Tipo de DVG
        BP              (int)   [in]    Break point: último nucleótido del 
                                        primer segmento
    Returns:
        len_second_seqment   (int)   [out]   Longitud del 1er segmento del DVG
    """

    a = ['Deletion_reverse', 'Insertion_reverse', '5cb/sb']
    b = ['Deletion_forward', 'Insertion_forward', '3cb/sb']

    if DVG_type in a:
        len_second_seqment = RI
    
    elif DVG_type in b:
        len_second_seqment = len_wt - RI + 1
    else:
        len_second_seqment = 'Error'
    
    return len_second_seqment

def len_dvg(len_wt, DVG_type, BP, RI):
    """
    Longitud completa del dvg: len_first_seqment + len_second_seqment
    """
    len_first_seqment = first_seqment(len_wt, DVG_type, BP)
    len_second_seqment = second_seqment(len_wt, DVG_type, RI)

    len_dvg = len_first_seqment + len_second_seqment

    return len_dvg

def add_lengths_to_df(input_csv, len_wt):
    from Bio import SeqIO
    import pandas as pd
    """
    Lee el fichero de input y lo convierte en un df al que añade la longitud
    teórica del DVG.

    Args:
        input_csv   (str)   [in]    Path to the input
        reference   (str)   [in]    Path to the reference sequence in fasta format
    Returns:
        df  (pd.DataFrame)  [out]   Table with the theoric length of the DVG
                                    added
    """
    # csv --> df
    df = pd.read_csv(input_csv, header=0)

    df[['length_dvg']] = df.apply(lambda row : len_dvg(len_wt,
                        row['DVG_type'], row['BP'], row['RI']), axis=1)
    
    #write
    df.to_csv('tempLen.csv', index=0)
    
    return df

def prop_wt(df):
    """
    Calcula la proporción en la que se encuentra el genoma wt como diferencia de 
    la proporción de genomas defectivos total
    Args:
        df  (pd.DataFrame)  [in]    Dataframe accesible

    Returns:
        prop_wt (float) [out]   Proporción en la que se encuentra el wt en una
                                región común a todos los genomas del dataset
    """
    
    prop_wt = 1 - sum(df['proportion'])

    return prop_wt

def prop_per_len(prop, length):
    """
    Calcula el producto entre proporción y longitud del genoma
    """
    return prop * length

def depth_common_coord_total(df, N, len_reads, len_wt, proportion_wt):
    """
    Calcula el valor de cobertura media que tendrá el dataset en una región en 
    la que estén presentes todos los genomas
    """
    # product prop * length of each dvg
    dvg_proxlen = df.apply(lambda row : prop_per_len(row['proportion'], 
                row['length_dvg']), axis=1)
    # summatory of the products prop * length of all dvgs
    sum_dvg_proxlen = sum(dvg_proxlen)

    depth_common = (N * len_reads)/((proportion_wt * len_wt) + sum_dvg_proxlen)

    return depth_common

def calc_N_dvg(prop, depth_common, len_dvg, len_reads):
    """
    Calcula el número aproximado de reads (N_dvg) con el que se generarán los
    fq simulados de los genotipos individuales. Redondearemos al alza el valor
    por lo que el N_total_real será aproximado al N introducido como parámetro.

    Args:
        prop    (float) [in]    Proporción que representa el genoma en una 
                                región común
        N   (int)   [in]    Número de reads totales introducidas como parámetro
        len_reads   (int)   [in]    Longitud de las reads, introducido como 
                                argumento
        len_wt  (int)   [in]    Longitud del genoma de referencia

    Returns:
        N_dvg   (int)   [out]   Número de reads con las que generar el fq
                                simulado para que se ajuste a las proporciones
                                introducidas como argumento
    """
    
    N = (prop * depth_common * len_dvg)/len_reads

    return round(N)

def add_N_dvgs(df, N, len_reads, len_wt, proportion_wt):
    """
    Añade al df una columna con el número de reads a generar de cada genoma
    Args:
        df  (pd.DataFrame)  [in]    Dataframe accesible
        N   (int)   [in]    Número aproximado de reads con las que se va a 
                            generar el dataset. Argumento de entrada.
        len_reads   (int)   [in]    Longitud de las reads. Argumento de entrada.
        len_wt  (int)   [in]    Longitud del genoma de referencia o wt
        
    Returns:
        df  (pd.DataFrame)  [out]   Tabla con la columna de N_dvg añadida
    """

    depth_common = depth_common_coord_total(df, N, len_reads, len_wt, proportion_wt)

    df[['N_dvg']] = df.apply(lambda row : calc_N_dvg(row['proportion'],
             depth_common, row['length_dvg'], len_reads), axis=1)
    
    return df


#------------------------------------------------------------------------------
## Generación de los fasta de cada dvg
#------------------------------------------------------------------------------

def complement(seq):
    """
    Genera la secuencia complementaria de seq
    Args:
        seq (str)   [in]    Secuencia de la que queremos obtener su 
                            complementaria
    """
    dict_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    seq = seq.upper()
    complement = []
    for nt in seq:
        complement.append(dict_complement[nt])
    return "".join(complement)

def reverse_complement(seq):
    """
    Genera la secuencia complementaria y reversa
    Args:
        seq (str)   [in]    Secuencia de la que queremos obtener su 
                            complementaria
    Returns:
        reverse_complement  (str)   [out]   Secuencia reversa y complementaria
    """
    seq = seq.upper()

    return complement(seq)[::-1]

def generate_sequence(reference_sequence, BP, RI, DVG_type):
    """
    Genera la secuencia teórica del dvg
    Args:
        reference_sequence  (str)   [in]   Secuencia de referencia que usaremos como
                                    molde.
        BP  (int)   [in]    Posición en base 1 del último nucleótido (incluído)
                            que forma parte del primer segmento del dvg
        RI  (int)   [in]    Posición en base 1 del primer nucleótido (incluído)
                            que forma parte del segundo segmento del dvg
        DVG_type    (str)   [in]    Tipo de DVG:
                            Deletion_forward, Deletion_reverse, Insertion_forward,
                            Insertion_reverse, 5cb/sb y 3cb/sb
                            En función del tipo de DVG se construye la 
                            secuencia.

    Returns:
        dvg_sequence    (str)   [out]   Secuencia completa del dvg
    """
    plus_plus = ['Deletion_forward', 'Insertion_forward']
    minus_minus = ['Deletion_reverse', 'Insertion_reverse']
    plus_minus = ['5cb/sb']
    minus_plus = ['3cb/sb']

    if DVG_type in plus_plus:
        first_seqment = reference_sequence[:BP]
        second_seqment = reference_sequence[RI-1:].lower()

    elif DVG_type in minus_minus:
        first_seqment = reverse_complement(reference_sequence[BP-1:])
        second_seqment = reverse_complement(reference_sequence[:RI]).lower()

    elif DVG_type in plus_minus:
        first_seqment = reference_sequence[:BP]
        second_seqment = reverse_complement(reference_sequence[:RI]).lower()
    
    elif DVG_type in minus_plus:
        first_seqment = reverse_complement(reference_sequence[BP-1:])
        second_seqment = reference_sequence[RI-1:].lower()

    dvg_sequence = first_seqment + second_seqment

    return str(dvg_sequence)

def basename_dvg(BP, RI, DVG_type):
    """
    Genera el nombre base o identificador del dvg. Es para evitar problemas 
    con el 'slash' en el nombre de los ficheros.
    Args:
        BP  (int)   [in]    Posición en base 1 del último nucleótido (incluído)
                            que forma parte del primer segmento del dvg
        RI  (int)   [in]    Posición en base 1 del primer nucleótido (incluído)
                            que forma parte del segundo segmento del dvg
        DVG_type    (str)   [in]    Tipo de DVG:
                            Deletion_forward, Deletion_reverse, Insertion_forward,
                            Insertion_reverse, 5cb/sb y 3cb/sb
                            En función del tipo de DVG se construye la 
                            secuencia.
    Returns:
        basename    (str)   [out]   Nombre base que se le va a dar a los 
                            ficheros individuales generados durante el proceso

    """
    basename = DVG_type + '_' + str(BP) + '_' + str(RI)
    cb = ['5cb/sb', '3cb/sb']
    if DVG_type in cb:
        basename = basename.replace('/sb', '')
    
    return basename

def add_basename_dvg_files(df):
    """
    Añade al df el nombre base que se le dará a los diferentes fichero generados
    (fasta y fastqs)

    Args:
        df  (pd.DataFrame)  [in]    Dataframe accesible

    Returns:
        df  (pd.DataFrame)  [in]    Tabla con los basename de cada dvg anotados
    """
    df[['basename_files']] = df.apply(lambda row : basename_dvg(row['BP'], 
                        row['RI'], row['DVG_type']), axis=1)
    
    return df
    

def write_sequence(reference_sequence, BP, RI, DVG_type):
    """
    Genera un fichero fasta con cada evento. El nombre del fichero se
    corresponde con el identificador del dvg

    Args:
        reference_sequence  (str)   [in]   Secuencia de referencia que usaremos
                            como molde.
        BP  (int)   [in]    Posición en base 1 del último nucleótido (incluído)
                            que forma parte del primer segmento del dvg
        RI  (int)   [in]    Posición en base 1 del primer nucleótido (incluído)
                            que forma parte del segundo segmento del dvg
        DVG_type    (str)   [in]    Tipo de DVG:
                            Deletion_forward, Deletion_reverse, Insertion_forward,
                            Insertion_reverse, 5cb/sb y 3cb/sb
                            En función del tipo de DVG se construye la 
                            secuencia.
    """
    dir_fastas = "Outputs/fastas/"
    basename = basename_dvg(BP, RI, DVG_type)

    fasta_file = dir_fastas + basename + '.fasta'

    header = '>' + basename
    dvg_sequence = generate_sequence(reference_sequence, BP, RI, DVG_type)
    
    # Escribimos la secuencia en un fichero fasta
    with open(fasta_file, 'w') as f:
        f.write(header + '\n' + dvg_sequence)

def write_all_dvg_sequences(reference_sequence, df):
    """
    Escribe los ficheros fasta de todas los dvgs que hay en la tabla
    Args:
        reference_sequence  (str)   [in]   Secuencia de referencia que usaremos
                                    como molde.
        df  (pd.DataFrame)  [in]    Tabla cargada en memoria y accesible para
                                    su acceso
    """
    # Genera la secuencia de cada dvg del df y escribe en un fasta
    df.apply(lambda row : write_sequence(reference_sequence, row['BP'], 
                        row['RI'], row['DVG_type']), axis=1)



#------------------------------------------------------------------------------
## Simulación de las reads con WGSIM
#------------------------------------------------------------------------------

def simulation_illumina_sequencing_wgsim(name, BP, RI, DVG_type, N, len_reads, 
                                        outer_distance):
    """
    Lanza la simulación de la secuenciación (short read simulation) con WGSIM
    con los parámetros introducidos y escribe los fastq (paired ende) de la 
    secuencia dvg.
    En esta versión ponemos todos los argumentos de mutación/error... a 0 
    Nota: Para ver todos los parámetros posibles consultar directamente la 
            ayuda de wgsim

    Args:
        BP  (int)   [in]    Posición en base 1 del último nucleótido (incluído)
                            que forma parte del primer segmento del dvg
        RI  (int)   [in]    Posición en base 1 del primer nucleótido (incluído)
                            que forma parte del segundo segmento del dvg
        DVG_type    (str)   [in]    Tipo de DVG:
                            Deletion_forward, Deletion_reverse, Insertion_forward,
                            Insertion_reverse, 5cb/sb y 3cb/sb
                            En función del tipo de DVG se construye la 
                            secuencia.
        N   (int)   [in]    Número de reads a generar
        len_reads   (int)   [in]    Longitud de las reads a generar
        outer_distance  (int)   Separación externa entre los dos finales de 
                            cada par generado. Realmente no nos interesa porque
                            solo nos quedaremos con una de los pares (los
                            trataremos con single ends)

    """
    dir_fastas = "Outputs/fastas/"
    dir_fastqs = "Outputs/fastqs/" + name + "/"

    basename = basename_dvg(BP, RI, DVG_type)
    fasta_file = dir_fastas + basename + '.fasta'
    
    fastq1 = dir_fastqs + basename + '_N' + str(N) + '_1.fq'
    fastq2 = dir_fastqs + basename + '_N' + str(N) + '_2.fq'
    
    # establecemos en 0 las variables relacionadas con generación de mutación, 
    # error o indels. 
    mut_error0 = "-e0 -r0 -R0 -X0"

    # Lanzamos la simulación
    os.system('wgsim -N{} -1{} -2{} -d{} {} -h {} {} {}'.format(N, len_reads,
    len_reads, outer_distance, mut_error0, fasta_file, fastq1, fastq2))

def generate_simulations(name, df, len_reads, outer_distance):
    """
    Genera las simulaciones de todos los dvgs presentes en la tabla introducida
    como input. El número de reads (N_dvg) vendrá determinado por el cálculo 
    realizado en el paso 1 del programa
    Args:
        df  (pd.DataFrame)  [in]    Tabla anotada con la columna 'N_dvg' y
                                    accesible como pd.DataFrame
        len_reads   (int)   [in]    Longitud de las reads a generar
        outer_distance  (int)   Separación externa entre los dos finales de 
                            cada par generado. Realmente no nos interesa porque
                            solo nos quedaremos con una de los pares (los
                            trataremos con single ends)
    """
    # simulations of each DVG event
    df.apply(lambda row : simulation_illumina_sequencing_wgsim(name, row['BP'], 
                    row['RI'], row['DVG_type'], row['N_dvg'], len_reads, 
                                        outer_distance), axis=1)
    # simulation of the native genome in the correct proportion
    
 
def simulation_wt(name, df, reference, N, len_reads, len_wt, proportion_wt, 
                outer_distance):
    """
    Lanza la simulación del genoma nativo con el número de reads 
    correspondientes. ///Mejorar porque calculo cosas que ya había calculado

    Args:
        df
        reference
        N
        len_reads
        len_wt
        proportion_wt
        outer_distance
    """
    dir_wt = "Outputs/fastqs/" + name + "/"
    depth_common = depth_common_coord_total(df, N, len_reads, len_wt, 
                    proportion_wt)
    N_wt = calc_N_dvg(proportion_wt, depth_common, len_wt, len_reads)



    # establecemos en 0 las variables relacionadas con generación de mutación, 
    # error o indels. 
    mut_error0 = "-e0 -r0 -R0 -X0"

    wt1 = dir_wt + 'wt_N' + str(N_wt) + '_1.fq'
    wt2 = dir_wt + 'wt_N' + str(N_wt) + '_2.fq'
    
    # run simulation with wgsim
    os.system('wgsim -N{} -1{} -2{} -d{} {} -h {} {} {}'.format(N_wt, len_reads,   
    len_reads, outer_distance, mut_error0, reference, wt1, wt2))

                
def list_fastqs(df, name, N_wt, as_single_end=True):
    """
    Extrae la lista de fastqs a concatenar. Por defecto solo usaremos los 
    {basename}_1.fq, es decir, uno de los ficheros generados en la simulación.
    Args:
        df
        as_single_end   (bool)  [in]    Cuando sea True, solo tendremos en 
                                cuenta el fq 1 de cada par. [Default True]
    
    Returns:
        fq_list (list)  [out]   Lista con todos los nombres de los ficheros 
                                fq que van a componer el fichero final
    """
    # si no existe el directorio fastqs del proyecto lo generamos pero si 
    # existe se reescribirá
    dir_fastqs = "Outputs/fastqs/" + name + "/"
    
    # wt fastqs 
    wt1 = dir_fastqs + 'wt_N' + str(N_wt) + '_1.fq'
    wt2 = dir_fastqs + 'wt_N' + str(N_wt) + '_2.fq'
    if as_single_end:
        end = ['_1.fq']
    else:
        end = ['_1.fq', '_2.fq']
    
    basename_list = list(df['basename_files'])
    N_list = list(df['N_dvg'])
    fq_list = []
    # Add {basename}_1.fq to the list
    for basename in basename_list:
        fq = dir_fastqs + basename + '_N' \
            + str(N_list[basename_list.index(basename)]) + end[0]
        fq_list.append(fq)
    fq_list.append(wt1)
    if len(end) == 2:
        for basename in basename_list:
            fq2 = dir_fastqs + basename + '_N' \
                + str(N_list[basename_list.index(basename)]) + end[1]
            fq_list.append(fq2)
        fq_list.append(wt2)
    return fq_list


def concatenate_fastq(df, N, len_reads, len_wt, proportion_wt, name,
                     as_single_end=True):
    """
    Escribe en un único fichero fastq ({name}.fq) todas las reads generadas en 
    la simulación. Este fq es el output final y debe contener aproximadamente
    el número de reads introducidas como parámetro por el usuario en las
    proporciones indicadas en el 'input.csv'
    
    Args:
        df  (pd.DataFrame)  [in]    Tabla con la información necesaria, tiene 
                                que tener el nombre base de cada evento 
                                anotado.
        N_wt (int)  [int]   Number of reads of the wild type genome
    """

    name_final_fq = "Outputs/" + name + '.fq'
    depth_common = depth_common_coord_total(df, N, len_reads, len_wt, 
                    proportion_wt)
    N_wt = calc_N_dvg(proportion_wt, depth_common, len_wt, len_reads)

    # Extracting the list of fq to concatenate
    fq_list = list_fastqs(df, name, N_wt, as_single_end=True)
    # preparing the list to the command line execution
    fq_list2command = " ".join(fq_list)

    # Writing a unique fq with all the individual fastqs:
    os.system('cat {} > {}'.format(fq_list2command, name_final_fq))

