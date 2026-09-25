ImmunoFluorescence_data_structure.mat contains a data structure named "Sessions" that contains raw gap ( Hb, Kr, Gt and Kni) and pair rule (Eve) gene expression (protein) intensity profiles measured using immunofluorescence imaging for fly embryos during n.c. 14 for 6 fly lines from our Bicoid-GFP fly line library. 

Each record in the structure contains the following fields:
1) FlyLineNumber: fly line number as listed in column 1 in Table S1 of Liu et al. PNAS 2013.
2) FlyLineName: fly line name as listed in column 2 of Table S1 of Liu et al. PNAS 2013.
3) GeneName: Gene name of the proteins immunostained and imaged in the measurement.
3) Embryos: structure that contains the source data for each embryo from one representative imaging session. 
---Embryos.dorsal: contains a n*1000 matrix with the ith column representing the dorsal profile of gene i listed in Genename. The 1000 points are equally spaced along the AP axis. NaN means that profile intensity was not reliably detected at that position (usually near the edges).
---Embryos.ventral: contains the ventral profile with the same formate as Embryos.dorsal.
---Embryos.InvaginationDepth: contains the furrow canal depth (delta_FC) measured in micrometers.
--Embryos.AP: contains information about the azimuthal orientation of the embryo on the slide (AP=1 if the confocal plane is closer to midsagittal plane, AP=0 if it is closer to the coronal plane)
---Embryos.EggLength: egg length in micrometer.
