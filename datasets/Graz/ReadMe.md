%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Graz Dataset
#### Dataset reference: 
"Leeb, R., Brunner, C., Müller-Putz, G., Schlögl, A., & Pfurtscheller, G. J. G. U. O. T. (2008). BCI Competition 2008–Graz data set B.
Graz University of Technology, Austria, 16, 1-6."

#### Training AT
Link for training data: <a href="http://www.bbci.de/competition/iv/"> www.bbci.de/competition/iv/ </a> -> Download of data sets -> agree submit -> Data sets 2a: ‹4-class motor imagery>
gdf files

#### Test AE
Link for test data: <a href="http://www.bbci.de/competition/iv/"> www.bbci.de/competition/iv/ </a> -> News -> Results of the

#### Labels
BCI Competition IV -> True Labels of Competition's Evaluation Sets -> Data sets 2a: 

Link to BioSig:   -> <a href="https://biosig.sourceforge.net/download.html"> biosig.sourceforge.net/download.html </a>

addpath -> path to ('biosig/t200_FileAccess/');
addpath -> path to ('biosig/t250_ArtifactPreProcessingQualityControl/');

# Specifics:

- 1:    Left Hand
- 2:    Right Hand
- 3:    Foot
- 4:    Tongue
