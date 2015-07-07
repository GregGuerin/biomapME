#Example for range.metrics.R

source("range.metrics.R")

data.eg <- structure(list(SPECIES = structure(c(7L, 7L, 7L, 7L, 7L, 7L, 
7L, 7L, 7L, 7L, 32L, 32L, 32L, 32L, 32L, 32L, 32L, 32L, 32L, 
32L), .Label = c("Atriplex acutibractea", "Atriplex acutibractea subsp. acutibractea", 
"Atriplex acutibractea subsp. karoniensis", "Atriplex acutiloba", 
"Atriplex angulata", "Atriplex australasica", "Atriplex cinerea", 
"Atriplex cordifolia", "Atriplex crassipes", "Atriplex crassipes var. crassipes", 
"Atriplex cryptocarpa", "Atriplex eardleyae", "Atriplex eichleri", 
"Atriplex elachophylla", "Atriplex fissivalvis", "Atriplex holocarpa", 
"Atriplex humifusa", "Atriplex incrassata", "Atriplex intermedia", 
"Atriplex kochiana", "Atriplex leptocarpa", "Atriplex limbata", 
"Atriplex lindleyi", "Atriplex lindleyi subsp. conduplicata", 
"Atriplex lindleyi subsp. inflata", "Atriplex lindleyi subsp. lindleyi", 
"Atriplex lindleyi subsp. quadripartita", "Atriplex lobativalvis", 
"Atriplex macropterocarpa", "Atriplex morrisii", "Atriplex nessorhina", 
"Atriplex nummularia", "Atriplex nummularia subsp. nummularia", 
"Atriplex nummularia subsp. omissa", "Atriplex nummularia subsp. spathulata", 
"Atriplex obconica", "Atriplex paludosa", "Atriplex paludosa subsp. cordata", 
"Atriplex paludosa subsp. paludosa", "Atriplex papillata", "Atriplex prostrata", 
"Atriplex pseudocampanulata", "Atriplex pumilio", "Atriplex quadrivalvata", 
"Atriplex quadrivalvata var. quadrivalvata", "Atriplex quadrivalvata var. sessilifolia", 
"Atriplex quinii", "Atriplex semibaccata", "Atriplex spongiosa", 
"Atriplex stipitata", "Atriplex suberecta", "Atriplex turbinata", 
"Atriplex undulata", "Atriplex velutinella", "Atriplex vesicaria", 
"Atriplex vesicaria subsp. calcicola", "Atriplex vesicaria subsp. macrocystidia", 
"Atriplex vesicaria subsp. variabilis"), class = "factor"), LONGITUDE = c(138.50138, 
137.790562, 133.683333, 138.616667, 139.584436, 137.783333, 137.783333, 
138.516667, 138.516667, 137.533333, 131.7, 135.416667, 135.416667, 
138.722222, 138.949722, 135.166667, 135.483056, 136.65, 135.25, 
139.4), LATITUDE = c(-35.077061, -35.820682, -32.133333, -35.55, 
-36.065392, -35.816667, -35.816667, -35.5, -35.5, -33.983333, 
-27.1, -27.5, -27.5, -35.4225, -26.933333, -32.85, -26.466389, 
-29.233333, -26.066667, -35.2)), .Names = c("SPECIES", "LONGITUDE", 
"LATITUDE"), row.names = c(15L, 17L, 19L, 28L, 33L, 39L, 40L, 
45L, 48L, 58L, 1L, 23L, 24L, 75L, 178L, 210L, 225L, 226L, 348L, 
396L), class = "data.frame")

range.metrics(data.eg)

range.metrics(data.eg, weight.type="geo")