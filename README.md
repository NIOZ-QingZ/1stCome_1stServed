# basic info
Update date: 31/10/2024
authors: Qing Zhan, Jonas Luca Mauch
License: CC-4.0 (reshare with acknowldegement of this GLEON project)

# 1stCome_1stServed
To gain mechanistic understanding about the potential role of the size as well as arrival time of inoculums, we developed a simple model to answer the question: to what extent, a phytoplankton species that arrives early in the system and with a large amount of biomass, can remain their dominance in the community.

# code:
Title: "Lotka-Volterra_1stCome1stServed.R"
Description: Model is based on Lotka-Volterra competition equations. 
Version 1: two species and four parameters including initial biomass, growth rates, competition effects, carrying capacity, species introudcing day.
The code includes:
    section 1: Lotka-Volterra model; 
        notes: code and dependencies for establishing the model
    section 2: sensitivity analysis; 
        section 2.1.: senstivity analysis of individual parameter
        section 2.2. sensitivity analysis of all parameters:
        notes: the sensitivity analysis takes more than 1hr to run. The simulation results are saved in the csv file entitled "sen_tIntro_K_P.csv". If you just want to play with the simulation results and don't do new sensitivy analysis, you can skip this section and jump to the next section about visualization
    section 3: visualization of senstivity analysis
        notes: package plot3Drgl can open a new device, and scatter3Drgl allows user to interact with the plot, meaning zooming, rotating, etc...


# data
Title :"sen_tIntro_K_P.csv"
Description: simulation data resulting in sensitivity analysis. Codes for generating the dataset can be found in the code entitled: "Lotka-Volterra_1stCome1stServed.R". Sensitivity analysis involve three parameteres: initial biomass, specias introducing day, and carrying capacity.

