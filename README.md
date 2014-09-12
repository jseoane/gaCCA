gaCCA  Genetic Algorithm feature selection CCA analysis
=====
This software was developed by
Jose A. Seoane, Ian Day, Colin Campbell, Juan Pablo Casas and Tom Gaunt

This software is based in:
-gene-based CCA:
A multivatiate test for Association. Ferreira M and Purcell S. Bioinformatics 2009
A gene based test of association using canonical correlation. Tang C and Ferreira M. Bioinformatics 2012
Original source can be found in Manuel Ferreira homepage:
http://genepi.qimr.edu.au/staff/manuelF/gene/main.html
-genalg R package:
Developed by Egon Willghagen
http://cran.r-project.org/web/packages/genalg/

The entry points are:
gaCCAgenMultipleGene_MultiplePhen_par.R: For multiple gene / multiple phenotype association
gaCCAgenMultipleGene_SinglePhen_par.R: For multiple gene / single phenotype association
gaCCAgenSingleGene_MultiplePhen_par.R For single gene / multiple phenotype association
gaCCAgenSingleGene_SinglePhen_par.R: For single gene / single phenotype association
genealg3.R modifies binary genalg using foreach for parallel execution in population fitness evaluation
genealgTwoCromosom.R: Mofifies the genealg behaviour to support two populations in the genetic algorithm. One for phenotype and othe for genotype.
gene-based.cca3.2.R: Modifies gene-based CCA improving some computation aspects.
