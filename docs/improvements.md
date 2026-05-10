# Homemade assemblies to database in step1

GTF/gff3, use cds value from column3. if the user only wants to use some cds, then make it clear on the readme to only include those in the GTF/gff3 file. 

In o2t step1:
Add a way to combine homemade assemblies to assemblies downloaded from ncbi, each dna file is for species so it would be just to add the seq from that homemade assembly/species to the db folder.


Luego, estaba pensando en que el usuario pueda implementar otra cosa ademas de CDS. 

Opcion para virus con mature peptides: que el script de python que extrae cds del gff se ajuste segun ocpions de extract mtature peptides o only matur epeptides.... Teniendo en cuenta tmb posibilida de agrupacion? mmmm

Quiero ver bien esa logica de agrupacion del script de python, que pasa si no hya protein id? o s i se usa otra cosa mm o transcirtp id mmm

Aguanta, pero los programas de anotacion generales no botan "mature_peptide" no? no lo hacen simplemente, entonces no es necesario añadir nada de eso en el parser del gff3? pero y si sí lo hacen de otro modo, tal vez programas d enaotaicon solo para virus, investigar


Me quedé en (viernes 8 de mayor):

okey, haz esto: 


Pendiente: arreglar minimamente los logs. Que si skipea solo accessiones, diga skipping step 1.3a. Si skipea local assemblies diga skipping step 1.3b. 

Revisar: que se integre bien con todos step2 y step3


Asunto de accession: 
En local assemblies:
Se usa stem de fasta file
En file con accessions:
Se usa la accesion

Sale en columna accession de og_genes.tsv


En paso 3 es obligatorio colocar en metadata columna accesiones. Pero no se verifica contra anda, es solo informativa del id. Y termina llamandose sample_id.

Tener en cuenta que o2tview verifica que ninguna columna se llame source y no recuerdo qué mas verificaba.



Que siempre salga la version actual d e omni2tree, se puede conectar de alguna forma a tag actual?

Colocar autores al final de cada help en o2t? 