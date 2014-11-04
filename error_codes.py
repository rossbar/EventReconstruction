SUCCESS = 0                                                                     
ORPHAN = -1     # Orphan = readout only worked on one side of the det           
ACONLY = -2     # No corresponding DC readouts                                  
DCONLY = -3     # No corresponding AC readouts                                  
NOT22_CS = -4   # If you have 4 readouts from 1 detector and there aren't 2 on  
NOT22_CL = -5   # each readout (i.e. 2 on AC, 2 on DC) then you likely have a   
NOT22_BOTH = -6 # big charge-sharing or charge-loss mechanism                   
ENERGY_MATCH_FAILURE = -7       # Energy matching fails for a 1-1 interaction   
THREE_INT_MATCH_FAILURE = -8    # Failed to find a 3-interaction event          
NOT_THREE_INTERACTIONS = -9     # Input not right format for a 3-interaction ev 
