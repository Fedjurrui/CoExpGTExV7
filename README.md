## CoExpGTExV7
Co-expression networks for the GTEx project in his seventh version [https://gtexportal.org/home/]. Here you will find annotated networks for 48 tissues of a total of 53. 

Each network is compound of:
1. Dataset corrected with the name of network. Type: RDS
2. The network for the tissue starting with the prefix net, then the name of the tissue. Type: RDS.
3. Enrichment results for the network modules using Gene Ontology, Reactome and KEGG, and the cell-type markers founds on each tissue. Type: TSV.


Recommended that you install first CoExpNets
```{r}
devtools::install_github('juanbot/CoExpNets')
```

And then this package
```{r}
devtools::install_github('Fedjurrui/CoExpGTExV7')
```

