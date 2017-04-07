This directory contains rewrite rules for all enzymes listed in the KEGG
LIGAND database[1] (Release 58.1 June 2011). The structure reflects the Enzyme
Nomenclature of the IUBMB, where E.C. numbers classify enzymes according to
chemistry, substrates, cofactors etc. For example, you would find the rewrite
rule(s) for E.C. 1.2.3.4 in the file 1/2/3/ec:1.2.3.4.

The KEGG database annotates all reactions catalyzed by an enzyme. If there is
more than one, all reactions are given in the GML file as separate rewrite
rules with the database ID as identifier.

The rewrite rules were built automatically from a reaction mapping based on
the Cut Successive Largest algorithm proposed in [2]. The generated atom maps
follow the Principle of Minimal Chemical Distance, which states that the most
likely reaction mechanism is the one requiring the least bond changes between
substrates and products. This may not always be the true mechanism.

[1] Kanehisa M, Goto S, Sato Y, Furumichi M, Tanabe M. KEGG for
integration and interpretation of large-scale molecular data sets. Nucleic
Acids Research. 2011:1-6.

[2] Crabtree JD, Mehta DP. Automated Reaction Mapping. Journal of
Experimental Algorithmics. 2009;13:1.15--1.29. 
