<!-- PROJECT INFO -->
<br />
<p align="center">

  <h3 align="center">SDgenerator</h3>

  <p align="center">
    SDgenerator is a simple program to simulate a FASTQ, the product of 'sequencing' a sample conformed of a wild type virus and a known population of DVGs. 
    <br />
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
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

<!-- BUILT WITH -->
## Built with
### wgsim for reads simula
 
  

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


<!-- CONTACT -->
## Contact

María José Olmo-Uceda - maolu@alumni.uv.es
PhD student

*EvolSysVir Group*

I^2SysBio (CSIC-UV)

<!-- MARKDOWN LINKS & IMAGES -->
