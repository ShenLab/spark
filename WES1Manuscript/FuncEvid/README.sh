# 
# Layer
#
# Reformat layer data
perl scripts/refmt_layer_data.pl

# Plot layer expression
Rscript Layer/HumanLayers.R

#
# Cell types
#
Rscript CellType/allKI.R

