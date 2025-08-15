# LBM_lid_driven_flow_in_polygons
Lattice Boltzmann MATLAB code for 2D lid-driven flow in polygons 

ME7751 Term Project Lid driven flow in a circular cavity using Lattice Boltzmann Method, ME Soohwan Kim
 c7  c3   c6  D2Q9 model
   \  |  /    
 c4 -c1 - c2  
   /  |  \    
 c8  c5   c9      
D2Q9 solver, SRT (Single Relaxation Time) model, BGK (Bhatnagar-Gross-Krook) collision operator
based on Textbook (Mohamad 2011) with clearer distinction between the physical,
nondimensional and numerical/discrete parameters, according to Jonas Latt's 2008 paper on the topic.

add function that monitors the current flow field in real-time while MATLAB is solving for the flowfield
