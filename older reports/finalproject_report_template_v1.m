%% Microgrid Design and Simulation
% _GROUP NAMES_ 
%
% *PROJECT TITLE* 



%% Microgrid Description and Overview
% A few sentences describing your project.  Please include images of your
% project location, including an overhead satellite image.
% Below is an example of how to include an image.  The image file should be
% in a folder named "html" at the same level as this script.
% You should have at least two descriptive images, one of which is an
% overhead satelillite image (2 pts).
% You should have a description of the location, some of its
% characteristics, and why it might benefit from a microgrid (2 pts).
% An additional 1 point for neatness and completeness.
%
% Here is my sdewscription of our microgrid.  Its in this location and so
% forth.
%
% <html>
% <img src="water_kiosk.png" width="75%">
% </html>



%% Load Calculations
% Start with a few sentence description of the load and its characterstics.
% You'll need to determine the daily kWh of your load, and the load factor.
% (The load factor is the ratio of daily average power to daily peak power;
% see the instructions for more discussion). Include all assumptions and
% calculations.
% 2 pts for a good description of how you arrived at your estimated daily
% kWh load.
% 2 pts for a description of how you arrived at an estimate of load factor.
% 1 pt for neatness and completeness.



%% Summary of Homer Output
% Provide a short summary of the design from Homer you are choosing to
% simulate in MATLAB/Simulink.  (Choose a design with solar and wind, if
% possible.)  
% 
% Design the energy storage secondary control, and the deferred load
% control.  The energy storage secondary control should manage the SOC in
% some way.  The deferred load control is optional, but gives you a tool 
% to keep the grid stable under the condition that the wind and solar
% generation exceed the load, and the battery is fully charged so that it
% cannot take any additional power.  (Hints: something like Pessec_ref =
% K*(0.5 - SOC) can work, as it acts to always drive the battery back to
% 0.5 SOC.  For the deferred load, you could create a check for which when
% wind and solar generation exceeds load, and the SOC is full or near full,
% we turn on the pumps.  Remember to avoid hard logic checks.  Try to use
% either smooth (proportional) control actions, or use relays.)
%
% Include a screenshot of the Homer output in your shared
% presentation.



%% Translating Homer's Outputs to MATLAB/Simulink
% Describe and demonstrate how you will translate the outputs of Homer to
% the MATLAB/Simulink parameters:
% 
% * meanLoad
% * peakLoad
% * es.Prated
% * es.Erated
% * wt.scale
% * pv.Prated
% * ht.Prated
% * mg.Pbase (should be set to either peak or average load)



%% Simulation Results
% Simulate your system and plot the following:
% 
% # Pload, Pwt, Ppv, Pes, Pht, Pload2
% # Pes on left axis and SOC on right axis
% # delfpu
%
% Comment on if the design is stable.
% Make three additional comments on anything of interest in the performance
% of your microgrid.



%% Final Project Instructions
%
% You can delete this section before the final publish (see publish instructions at the bottom).  
% 
% Design a microgrid (wind, solar, battery, load, and hydro) using Homer,
% and simulate your microgrid for 24 hours using MATLAB Simulink
%
% You'll first need to open a free account with Homer Quickstart
%
% You'll need the following information to design your microgrid:
%
% * Location of the microgrid and local cost of diesel fuel
% * Type of load (residential, commercial, industrial, community)
% * Amount of load (peak W in a day, and kWh per day)
% * Cost of energy storage
% * Cost of PV 
% * Cost of Wind
% * Cost of diesel fuel
%
% Cost of diesel: you can use the default of $1/liter.
%
% Electrical load:
%
% * Can use 11.3 kWh/sqft/yr for residential loads (EIA,
% https://www.eia.gov/consumption/residential/data/2015/c&e/pdf/ce1.1.pdf)
% * Can use 14.6 kWh/sqft/yr for commercial loads (EIA,
% https://www.eia.gov/consumption/commercial/data/2012/c&e/cfm/pba4.php)
% * Industrial loads are more specific to the site.  (I.e., what does the
% factory do?)  Make a reasonable estimate and justify it.  (One example
% approach would be to find a reference for a similar site and benchmark
% against that.)
% * As a sanity check, consider that a typical residential home will have
% 1500 to 2500 square feet, and an energy consumption of 30-40 kWh per day.
% A typical Fred Meyer will have tens of thousands of square feet and an
% energy consumption of up to thousands of kWh per day.
% * We also need to determine the load factor, which is the ratio of
% average daily power to peak daily power.  As a rough rule of thumb, very
% large loads distributed over very wide areas may have a load factor of
% around 0.85.  Smaller loads have smaller load factors.  An individual
% house might have a load factor of around 0.3.
% * From the load factor and the daily energy consumption (kWh), you can
% estimate the peak load.
% * IMPORTANT: after Homer has designed your grid, click on the "Load" in
% the "Schematic" and adjust the "variability" accordingly to match your
% peak load.  You'll need to re-calculate the design after changing this.
% 
% PV:
%
% * For the price of capacity, use
% https://www.nrel.gov/solar/market-research-analysis/solar-installed-system-cost.html
% * The y-axis is dollars per watt.
% 
% Wind turbine:
% 
% * Assume a 100 kW XANT turbine
% * Assume a capacity cost of $1/W for wind turbines
%
% Battery:
%
% * Use generic 1 kWh lithium ion batteries.
% * Use $140 per battery.
%
% Homer outputs:
%
% * Homer will provide a table of possible designs, ordered by Net Present
% Cost.  The columns are:
% * kW of diesel capacity
% * kW of PV capacity
% * kW of battery power capacity
% * Number of 1 kWh battery cells
%
% Rubric:
%
% * Description (5 pts)
% * Load Calculations (5 pts)
% * Homer Output Summary (5 pts)
% * Simulation Setup (5 pts)
% * Simulation Results and Analysis (5 pts)
% * Presentation (5 pts)
%
% *Note on publishing*
% For the html style image inclusion code to work, you should publish to
% html.  Any images should be a folder named "html" at the same level as
% this script.  When this script is "publish"ed, it will produce an html
% file in the html folder.  You can then print that html file to pdf using
% any browser.