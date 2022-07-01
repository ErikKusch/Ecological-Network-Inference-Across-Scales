# Ecological Network Inference across Scales
Chapter 3 of my PhD project.

## Abstract
1.	Ecosystem stability has been linked to the topology of ecological networks; an attribute difficult to quantify at macroecological scales due to a limited number of ecological networks inferred using field observations. Several methods for inference of ecological networks have been proposed to overcome this shortfall using distributional and population information, but the consistency among these remains unstudied.
2.	Here, we use two data sets describing distributional and attribute information at local and continental scales for woody species across North America. We subsequently infer ecological networks from these data using four different frameworks and compare the resulting ecological networks.
3.	Inferred network topologies were neither consistent among approaches when comparing networks belonging to the same scale, nor when comparing networks inferred by the same method across geographic scales. We did find some scale consistency for node-level topology metrics across approaches and scales.
4.	Choice of network inference framework, particularly at macroecological scales, significantly impacts inferred 1-1 associations and network-level topology metrics. We argue that the choice of inference framework needs to align the research questions with the framework applicability, network-link interpretation, geographic scale of assessment, and ecological processes governing community assembly at the assessment scale.  

## Research Questions
1. How consistent are the topological attributes of ecological networks derived from different methodologies when assessing the same ecosystem?   
2. How scale-dependent are the topological attributes of ecological networks derived using the same inference approach? 

## Primary Contact and Collaborators
### Primary Contact
Erik Kusch (erik.kusch@bio.au.dk, [OrcID](https://orcid.org/my-orcid?orcid=0000-0002-4984-7646))  
PhD Student  
Department of Biology  
Section for Ecoinformatics & Biodiversity  
Center for Biodiversity Dynamics in a Changing World (BIOCHANGE)  
Aarhus University  

### Collaborators
- Alejandro Ordonez ([OrcID](https://orcid.org/0000-0003-2873-4551))  
- Malyon Bimler ([OrcID](https://orcid.org/0000-0003-0059-2360))  
- James A. Lutz ([OrcID](https://orcid.org/0000-0002-2560-0710))

## Data
This research is focussed on two data sets of woody plants across North America:  

### Yosemite Forest Dynamics Plot (YFDP)  
These data are available upon reuqest to James A. Lutz.

### Forest Inventory Analysis Database (FIA)
These data have been obtained using the rFIA package (https://cran.r-project.org/web/packages/rFIA/index.html) .

### Data Processing
YFDP and FIA data have been matched with climate data using the KrigR package (https://github.com/ErikKusch/KrigR) and subsequently been divided into subsets corresponding to (1) biomes across North America according to the WWF (https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world), (2) state-lines according to naturalearthdata (https://www.naturalearthdata.com/) as well as (3) the region of the yosemite national park according to the National Park Service (https://databasin.org/datasets/7ae51cf6-aac8-4636-b1db-eb206ed012d0/).

## Analyses
We infer ecological networks using four separate ecological network inference methods:
- COOCCUR (https://cran.r-project.org/web/packages/cooccur/index.html)
- NETASSOC (https://cran.r-project.org/web/packages/netassoc/index.html)
- HMSC (https://cran.r-project.org/web/packages/Hmsc/index.html)
- IF-REM (https://www.biorxiv.org/content/10.1101/2022.03.28.486154v1)

We do so for three distinct scales: (1) YFDP-scale (plot scale, (2) FIA data for the Yosemite National Park region (regional scale), and (3) FIA data for all Temperate Conifer Forests across North America (macro-scale).

We subsequently compare inferred ecological networks and their topologies among methods within scales and across scales within methods.

## Key Findings
Our research highlights a lack of consensus in inferred ecological networks across scales and methods. Only for node-level topology metrics did we find some consistency across methods and scales.

We argue that the choice of ecological network inference approach needs to consider: (1) alignment of link-interpretation with research questions, (2) applicability of the network inference method at the geographic scale of study, (3) data availability and information content at the geographic scale of study, and (4) biological mechanisms which are identifiable using network inference at the scale of study. 

## Funding
Aarhus University Research Foundation Strat-up Grant (grant no. AUFF-2018-7-8) 
