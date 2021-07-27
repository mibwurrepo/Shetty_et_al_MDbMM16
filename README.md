Functional analysis of a synthetic gut microbiome reveals niche segregation as a driving force for assembly, stability, co-existence and activity       
----------------------------------------------------------------------------------------------------------------------

Authors: Sudarshan A. Shetty<sup>1#⸸*</sup>, Ioannis Kostopoulos<sup>1#⁋</sup>, Sharon Geerlings<sup>1#</sup>, Hauke Smidt<sup>* 1</sup>, Willem M. de Vos<sup>1,2* </sup>, Clara Belzer<sup>*1</sup>  

<sup>1</sup>Laboratory of Microbiology, Wageningen University & Research   
<sup>2</sup>Human Microbiome Research Program, Faculty of Medicine, University of Helsinki, Helsinki, Finland   
<sup>⸸</sup>Present address: University Medical Center Utrecht, Utrecht, The Netherlands & Rijksinstituut voor Volksgezondheid en Milieu, Bilthoven, The Netherlands  
<sup>⁋</sup>Present address: Danone Nutricia Research, Utrecht, The Netherlands  
<sup>#</sup>Shared first authors; These authors contributed equally.  
<sup>*</sup>Co-Corresponding authors; These authors contributed equally.  
  

Manuscript Correspondence:  
Sudarshan A. Shetty; sudarshanshetty9[@]gmail[.]com  
Clara Belzer; clara.belzer[@]wur[dot]nl   
Hauke Smidt; hauke.smidt[@]wur[dot]nl   
Willem M. de Vos; willem.devos[@]wur[dot]nl    

This repository contains codes for analysis done in the research article by Sudarshan A. Shetty<sup>1#⸸*</sup>, Ioannis Kostopoulos<sup>1#⁋</sup>, Sharon Geerlings<sup>1#</sup>, Clara Belzer<sup>*1</sup>, Hauke Smidt<sup>*1</sup>, Willem M. de Vos<sup>1,2</sup>* (2021).Functional analysis of a synthetic gut microbiome reveals niche segregation as a driving force for assembly, stability, co-existence and activity. [in prep](XXX).   

**We developed `syncomR` package as a research compendium for organizing and sharing files processed files and custom functions**   
For re-running the analysis, first install `syncomR` package from [GitHub](https://github.com/microsud/syncomR)  

```r
install.packages("devtools")

devtools::install_github("microsud/syncomR")

```
The scripts that we used for processing raw reads are in `amp_seq/` (dada2, raw reads to ASVs --> transferred to syncomR) and `rna_seq/` (modified samsa2, raw reads to gene count tables --> transferred to syncomR). The main analysis codes and figures reported in the manuscript are located in `analysis/`, `codes/`, and `data/` folders.  

Correspondence for analysis:  
Author: Sudarshan A. Shetty
Contact: sudarshanshetty9[at]gmail[dot]com  
