----------------------------------------------------
## Debris Flow Numerical Simulation — Basilisk C
----------------------------------------------------
## Results

<p align="center">
  <img src="Example : Results.png" width="45%">
</p>

## DESCRIPTION

Numerical simulation of debris flow propagation <br>
Using the Saint-Venant shallow-water equations with Voellmy friction law, on an adaptive quadtree grid (AMR).


This code was developed by students for a final year Master's project<br>
It is an adaptation of an existing Gerris simulation, transposed to Basilisk C<br>
Some parts were improved with the assistance of AI<br>



## REQUIREMENTS

- [Basilisk C] (http://basilisk.fr) installed
- Input rasters in ESRI ASCII format (.asc)
- A 'terrain' folder in your Basilisk working directory
- An empty 'results' folder in your Basilisk working directory


## PROJECT STRUCTURE

```text
project/
├── code.c
├── README.md
├── simulation_log.txt
├── terrain/
│ ├── without_rockfall.asc
│ ├── only_rockfall.asc
│ └── DoD.asc
└── results/
```

## INPUT FILES

Place the following files into the 'terrain' folder in your Basilisk working directory : 

|          File        |                   Description                     |       
|----------------------|---------------------------------------------------|      
|'without_rockfall.asc'| DEM without the rockfall deposit (bed topography) |
|'only_rockfall.asc'   | Release area - Initial flow depth		            |
|'DoD.asc' (optional)  | Difference of DEMs - for post processing          |



## SET UP FOR A NEW PROJECT 
Mandatory step for a new project (otherwise, skip to the next section)

1. Domain origin coordinates (TX, TY)
These must match the bottom-left corner - replace with your xllcorner (x lower left corner). Same for your yllcorner. 

2. Input files with the same names 

3. Parameters 
MU, XI, END, LEVEL : adjust as needed



## SET-UP FOR THIS PROJECT - Debris Flow
1. Folder terrain with the input files
2. Create a 'results' folder



## COMPILATION 

- Standard (single core) :<br>
   bash : **qcc debris_flow.c -O2 -o debris_flow -lm**

- Faster (multi-core) : <br>
1. Check how many cores are available <br>
   bash : **nproc**
2. Compile with the desired number of threads (here 8 out of 12) <br>
   bash : **qcc -O2 -fopenmp debris_flow.c -o debris_flow -lm && OMP_NUM_THREADS=8**



## RUNNING THE SIMULATION 

bash : **./debris_flow**



## OUTPUT 

Results are saved in the 'results/' folder (must be created before running)

1. 'localDepth-mu***-xi***-L**-t***.asc' final flow depth (Name encoded with mu, xi, level and time parameters)
2. 'simulation_log.txt' containing the run metadata



## KEY PARAMETERS SUMMARY 

| Parameter | Default | Description                            |
|-----------|---------|----------------------------------------|
| `MU`      | 0.155   | Coulomb friction coefficient           |
| `XI`      | 100     | Turbulent friction coefficient (m/s²)  |
| `LEVEL`   | 11      | Max refinement level                   |
| `END`     | 600     | Simulation duration (s)                |
| `DRY`     | 1e-4    | Dry cell threshold (m)                 |



## WORKFLOW SIMPLIFIED

```text
main()
  │
  ├── Domain Initialisation
  │     ├── size(), origin()
  │     └── parameters (G, CFL, DT)
  │
  ├── run()
  │     │
  │     ├── event init (t = 0)
  │     │     ├── reads raster (.asc)
  │     │     ├── grid projection
  │     │     ├── slope calculation (dzdx, dzdy)
  │     │     └── initial condition (h, u)
  │     │
  │     ├── Time loop (i++)
  │     │     │
  │     │     ├── Saint-Venant solver 
  │     │     │
  │     │     ├── event voellmy (i++)
  │     │     │     └── add friction
  │     │     │
  │     │     ├── event adapt (i++)	
  │     │     │     └── AMR : refinement / coarsening
  │     │     │
  │     │     └── (optional) intermediate outputs 
  │     │
  │     ├── event diagnostics (t = END)
  │     │     └── localDepth computation
  │     │
  │     └── event end (t = END)
  │           ├── final raster export (.asc)
  │           └── update simulation_log.txt
  │
  └── End of program
```


## SIMULATION TIME

With 8 threads it usually takes : 
- For a LEVEL = 11, END = 600 (s) >> ~ 20 minutes
- For a LEVEL = 10, END = 300 (s) >> ~ 3 minutes



## POST-PROCESSING 

- Visualisation >> QGIS (.asc is already georeferenced)
- Volume calculation >> Python (comparison with the observed DoD volume)


------------------------------------------------------------------------

## AUTHOR

**Paul SIESS** - paul.siess39@gmail.com <br>
Master's student - ENSEGID <br>
Final year project 2025/2026 <br>


Supervised by: **Anne-Laure ARGENTIN** 
