<!-- PROJECT INFO -->
<br />
<p align="center">

  <h3 align="center">SDgenerator v2</h3>
  <p align="center">
    SDgenerator is a python program to create synthetic samples (in FASTQ format) with a controled population of DVGs (proportion & type) 
    
The version 2 consider the case that all the genomes in the sample doesn't share a common genomic coordinate. The way to make the v2 calculation was carried by @jmusan
    
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">External Program</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>

<!-- EXTERNAL PROGRAM -->
## Third party program

`wgsim` (reads simulator)

<https://github.com/lh3/wgsim>
  

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.


### Installation

1. Clone the repo in the directory of your choice
   ```sh
   git clone https://github.com/MJmaolu/SDgenerator.git
   ```

2. Go to the DVGfinder directory
   ```sh
   cd SDgenerator
   ```
      
3. Give execution permission to all the scripts in the SDgenerator directory
   ```sh
   chmod -R +x .
   ```
5. Create a new environment with conda with all the dependencies needed to run SDgenerator
   ```sh
   conda env create -f sdgenerator_env.yaml
   ```
   
6. Activate SDgenerator environment 
   ```sh
   conda activate sdgenerator_env
   ```

<!-- USAGE EXAMPLES -->
## Usage

```python3 SDgenerator.py -t population_csv -N number_total_reads [-l length_reads. Default 100] -o output_basename [-d outer_distance. Default 300]```

<!-- INPUT FORMATS -->
## Input formats

### Population_csv example

BP,RI,DVG_type,proportion      
6194,8961,Deletion_forward,0.1       
23000,120,Deletion_reverse,0.1      
9478,9489,3cb/sb,0.1       
10000,9800,5cb/sb,0.1        
2200,2100,Insertion_forward,0.1     
1000,1200,Insertion_reverse,0.1

The proportion of the wild type genome is calculated as the substraction to 1 of the sum of the DVG proportions


### Supported DVG types 

- Deletion_forward

- Deletion_reverse

- 3cb/sb

- 5cb/sb

- Insertion_forward

- Insertion_reverse

<!-- CONTACT -->
## Contact

María José Olmo-Uceda - maolu@alumni.uv.es

PhD student

*EvolSysVir Group*, I<sup>2</sup>SysBio (CSIC-UV) 

<!-- MARKDOWN LINKS & IMAGES -->
