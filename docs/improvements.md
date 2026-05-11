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

Eliminar este info:
[2026-05-11 12:33:30] [INFO] Writing step1 resume state (omni2tree_step1_resume.json, complete=false)?

Eliminar este info:
[INFO] Taxon HRSVAUSATXHOU1986702021 is marked for repair; keeping it in the download queue.
[2026-05-11 12:36:35] [INFO] Taxon hRSVBAustraliaVICRCH0562019 is marked for repair; keeping it in the download queue.
[2026-05-11 12:36:35] [INFO] Taxon RSVBhumanUSAWAS253042021 is marked for repair; keeping it in the download queue.
[2026-05-11 12:36:35] 

Que no repita para todos : [INFO] Validated 11 local sequence unit(s) from 11 feature row(s) for local taxon HRSV-A/IC688/12: /home/adrian_chalco/projects/fritz/o2t_local_assemblies/data/local_assemblies/KP663728.1.gff3
[2026-05-11 12:36:34] [INFO] Validated 11 local sequence unit(s) from 11 feature row(s) for local taxon HRSV/A/England/397/2017: /home/adrian_chalco/projects/fritz/o2t_local_assemblies/data/local_assemblies/PP109421.1.gff3
[2026-05-11 12:36:34]  

Que siempre salga la version actual d e omni2tree, se puede conectar de alguna forma a tag actual?

Colocar autores al final de cada help en o2t? 