# PEDIA analysis in Translate NAMSE study

This repository hosts the code utilized in the PEDIA analysis within the Translate NAMSE study titled 
"Next-generation phenotyping integrated in a national framework for patients with ultra-rare disorders improves genetic diagnostics and yields new molecular findings".
You can find the preprint of the study here (https://www.medrxiv.org/content/10.1101/2023.04.19.23288824v1).

Our analysis leveraged the open-source versions of PEDIA and GestaltMatcher as outlined in the paper titled 
"Enhancing Variant Prioritization in VarFish through On-Premise Computational Facial Analysis" available here
(https://www.mdpi.com/2073-4425/15/3/370). 
Additionally, we benchmarked various HPO and variant analysis approaches including LIRICAL, AMELIE, Exomizer, and Xrare alongside PEDIA.

We provide the codes and instructions necessary for replicating the same analysis with your own cohort.
Please note that in order to ensure patient data protection, access to the training data is available upon request.
For inquiries regarding data access, kindly reach out to Dr. Tzung-Chien Hsieh at thsieh@uni-bonn.de.

## Step-by-step analysis
PEDIA integrates three distinct scores—GestaltMatcher, CADA, and CADD—into a unified PEDIA score, enhancing the prioritization of genes and variants.
Below are reference of each approach:
* GestaltMatcher: Facial image analysis (https://www.nature.com/articles/s41588-021-01010-x)
* CADA: clinical feature analysis (https://academic.oup.com/nargab/article/3/3/lqab078/6363753)
* CADD: exome variant analysis (https://www.nature.com/articles/ng.2892)
* PEDIA: integrating three different scores into PEDIA score (https://www.nature.com/articles/s41436-019-0566-2)

### VarFish: variants analysis platform with PEDIA
We used Varish (https://github.com/ahujameg/varfish-server) to enable the PEDIA analysis.
VarFish is a variant analysis tool that requires integration with the PEDIA middleware to connect with GestaltMatcher and PEDIA services.
The input of the patient are:
* Facial image
* Clinical features in HPO terms
* Exome variants in VCF

![VarFish](https://github.com/igsb/PEDIA-TNAMSE/blob/main/image/Enhancing_applicability_of_PEDIA.svg)

Please install the VarFish to enable the variant analysis.

VarFish requires integration with the PEDIA middleware to access GestaltMatcher and PEDIA services.
Refer to the documentation provided in the PEDIA middleware repository for instructions on connecting VarFish with the middleware 
(https://github.com/igsb/pedia-middleware).

Once you launch Varfish and import a VCF, you can perform GestaltMatcher and PEDIA shown in the screenshot below.
![PEDIA](https://github.com/igsb/PEDIA-TNAMSE/blob/main/image/GM-Varifish.png)

For detailed step-by-step instructions on analyzing a patient using VarFish and its integration with the PEDIA middleware, refer to the results section of the article associated with this study 
(https://www.mdpi.com/2073-4425/15/3/370).

## Results
After performing the PEDIA analysis, the resulting variants are displayed in the table below.
These variants can be further sorted by the PEDIA score in descending order.
In this patient, the disease-causing variant is identified within the MED13L gene, which holds the top-ranking position based on the PEDIA score.
![results](https://github.com/igsb/PEDIA-TNAMSE/blob/main/image/result_table.png)

After performing the PEDIA analysis, you have the option to export the results in Excel sheet format. Additionally, you can visualize the accuracy of your cohort by plotting figures using the provided Jupyter Notebook (pedia_revision.ipynb).

To protect the patient data, please contact us for accessing the data to reproduce the figures shown in results_figure folder. 

## Benchmark with different feature analysis approaches
In the Translate NAMSE (TNAMSE) study, various alternative feature and exome analysis approaches were explored in addition to the original CADA and CADD scores. These alternative approaches include LIRICAL, AMELIE, Exomizer, and Xrare.

Detailed scripts for running these approaches can be found in the 'run_prioritization_tools_on_vcfs' folder of this repository.

Since these tools are external and may require significant resources for installation, we recommend referring to their original documentation for detailed instructions on installation and usage. Please consult the documentation provided by the respective developers or organizations for comprehensive guidance on installing and utilizing LIRICAL, AMELIE, Exomizer, and Xrare.

## References
GestaltMatcher: https://www.nature.com/articles/s41588-021-01010-x

GestaltMatcher-Arc: https://openaccess.thecvf.com/content/WACV2023/html/Hustinx_Improving_Deep_Facial_Phenotyping_for_Ultra-Rare_Disorder_Verification_Using_Model_WACV_2023_paper.html

LIRICAL: https://www.cell.com/ajhg/fulltext/S0002-9297(20)30230-5

AMELIE: https://www.science.org/doi/10.1126/scitranslmed.aau9113

Exomizer: https://www.nature.com/articles/nprot.2015.124

Xrare: https://www.nature.com/articles/s41436-019-0439-8


## Contact
For inquiries regarding data access or any questions, kindly reach out to Dr. Tzung-Chien Hsieh at thsieh@uni-bonn.de.


