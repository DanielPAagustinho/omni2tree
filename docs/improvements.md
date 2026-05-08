# Homemade assemblies to database in step1

GTF/gff3, use cds value from column3. if the user only wants to use some cds, then make it clear on the readme to only include those in the GTF/gff3 file. 

In o2t step1:
Add a way to combine homemade assemblies to assemblies downloaded from ncbi, each dna file is for species so it would be just to add the seq from that homemade assembly/species to the db folder.


Luego, estaba pensando en que el usuario pueda implementar otra cosa ademas de CDS. 

Opcion para virus con mature peptides: que el script de python que extrae cds del gff se ajuste segun ocpions de extract mtature peptides o only matur epeptides.... Teniendo en cuenta tmb posibilida de agrupacion? mmmm

Quiero ver bien esa logica de agrupacion del script de python, que pasa si no hya protein id? o s i se usa otra cosa mm o transcirtp id mmm

Aguanta, pero los programas de anotacion generales no botan "mature_peptide" no? no lo hacen simplemente, entonces no es necesario añadir nada de eso en el parser del gff3? pero y si sí lo hacen de otro modo, tal vez programas d enaotaicon solo para virus, investigar


Me quedé en (viernes 8 de mayor):

okey, haz esto: 1. En o2t step1, en la parte d evlaidacion de gff, edita brevmeente el scirpt para no imprimir un mensaje por cada gff validado, sino decir algo como que se está validando, y luego que se valdio y botar algun warn si es necesario. 2. Luego, veo que siempre generas esas imagenes de cds dos veces, cuando el input es local y accessions downloaded, puedes evitar eso minimamente o seria muy pesado? piensa como evitarlo minimamente, teniendo en cuenta que el parametro resume debe mantener su comportamiento de siempre. Funcionara mover la generacion del png al ultimo? primero propon y dime si es factible. 

Pendiente: añadir "frame" en caracteristicas a recuperar en el gff3 columna 9 y tmb en el og generate tsv, si es que estan presentes en lo descargado desde ncbi. Ah no, pero si añado eso no generaria muchos NA en los archviso resumenes de OG? Bueno añadirlo al final entonces

Pendiente: hacer algo cuando el gen se llama "NA" tal cual

Pendiente: Protein y Product son practiamente lo mismo  par amis intereses .. asi qeu fine!