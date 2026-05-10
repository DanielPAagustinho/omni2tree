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

Pendiente: que en todo o2t step1.sh se trate NA como valor normal, y tmb se añada NA cuando no se sabe que valor es. Creo quee sto ya es asi, solo verificar que siga siendo asi. Añadir no ta en el readme para alertar al usuario con ese nombre de gen, sobre todo cuando venga descargado de ncbi

Y ademas a location solo pongo location, no añado protein id ni strand., solo start..end


Pendiente: arreglar minimamente los logs. Que si skipea solo accessiones, diga skipping step 1.3a. Si skipea local assemblies diga skipping step 1.3b. 

