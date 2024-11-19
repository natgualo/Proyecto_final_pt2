# Paso 1: Preparar las paqueterías necesarias, funciones y cargar la base de datos

#Se necesitan paquetes de Biostrings para hacer el manejo de secuencias 

#Cargar los paquetes con library seqinr, Biostrings y si no los tienes entonces instalarlos con 
#install.packages

#funcion para detectar el tipo de secuencia
detectar_tipo_secuencia <- function(secuencia) {
  # Convertir la secuencia a mayúsculas para que no haya problemas al momento de leer las secuecnias
  #si son en minusculas 
  secuencia <- toupper(secuencia)
  
  # Verificar si todos los caracteres están en A, T, C, G, U, incluso con el ARN 
  if (all(strsplit(secuencia, NULL)[[1]] %in% c("A", "T", "C", "G", "U"))) {
    return("ADN")
  } else {
    return("proteina")
  }
}

#funcion para convertir a AA si es secuencia de nucleotidos y quieres compararlos con aminoacidos

traducir_secuencia<- function(secuencia){
  #primero se correra la funcion que detecta el tipo de secuencia y guarda el tipo en una variable
  
  tipo<- detectar_tipo_secuencia(secuencia)
  
  #ahora con un if else checa que si tipo es ADN
  if(tipo=="ADN"){
    
    #primero se vuelve a cargar la secuencia pero como un objeto de DNA string. de la secuencia
    #cargada se extrae la secuencia con [1][[1]](el primer [] es para ver la primer secuencia y
    #el [[]] es para sacar los caracteres), eso se tranforma a tipo caracter con as.character
    #para ser cargado con la funcion DNAStringSet
    
    secuencia_problema<-DNAStringSet(as.character(secuencia[1][[1]]))
    
    # se traduce a proteina y renombramos la variable secuencia_problema(tuvimos que llamar directamente 
    #la funcion translate con Biostrings::translate porque habíamos tenido problemas con
    #el uso de la función)
    
    secuencia_problema<-Biostrings::translate(secuencia_problema)
  }else {
    #si si era de aminoacidos la secuencia, se queda igual
    return(secuencia)
  }
}
library(Biostrings)
library(seqinr)

##################################################################################################

# Cargar la base de datos (archivo FASTA)
base_datos <- readAAStringSet("Data/BacMet_EXP.BAC0000.BIOCIDES.325.fasta.txt")
base_datos

# Paso 2: Cargar nuestra secuencia problema
secuencia_problema <- readAAStringSet("Data/BacMet2_EXP_database.fasta.txt")
secuencia_problema

# Paso 3: Detectar si es secuencia de nucleótidos o aminoácidos y traducir de ser necesario

detectar_tipo_secuencia(secuencia_problema)
#Traducir secuencia si es de nucleotidos 
traducir_secuencia(secuencia_problema)

# Paso 4: hacer el alineamiento de las secuencias y crear dataframe de resultados
#creamos un dataframe vacio donde se van a ir añadiendo las filas de los resultados de cada
#alineamiento

comparassion_results<-c()

#con un ciclo for se llevara a cabo un alineamiento o "comparación" entre cada una de las secuencias
#de la base de datos 
#para cada i secuencia de la base de datos desde la primera hasta la ultima
for (i in 1:length(base_datos)) {
  #se extrae el nombre de la secuencia de la base de datos que se esta analizando. [i] indica que
  #es la i secuencia y el @ranges@NAMES es para extraer del objeto de DNAstring el nombre de la
  #secuencia
  nameseq<-base_datos[i]@ranges@NAMES
  
  #se realiza un alineamiento de la i secuencia de la base de datos contra la secuencia problema
  #con una apertura y extension de gap -10 y -1 y que sea un alineamiento local y el alineaiento
  #se guarda en una variable
  alignment<-pairwiseAlignment(base_datos[i],secuencia_problema,gapOpening = -10, gapExtension = -1, 
                               type = "local")
  
  #de los datos asociados al alineamiento se obtiene el score y se guarda en una variable
  score<-alignment@score
  
  #con la función pid se btiene el porcentaje de identidad
  p_identity<-pid(alignment)
  
  #en una variable llamada temp se concatenan cada una de las variables generadas y con t para
  #que las variables concatenadas se acomoden como una fia
  temp<-t(c(nameseq,score,p_identity))
  
  #por ultimo se añade al dataframe vacio creado antes del ciclo los resultados del alineamiento 
  #contra la i secuencia de la base de datos
  comparassion_results<-rbind(comparassion_results,temp)
}

View(comparassion_results)

#renombramos el nombre de las columnas del dataframe
colnames(comparassion_results)<-c("Protein_with_more_similarity","Score","P_identity")

#lo guardamos como una variable de tipo dataframe
comparassion_results<-as.data.frame(comparassion_results)

#con esta función convertimos cada uno de los elementos de la segunda hasta la ultima columna
#a tipo numerico
for (j in 2:length(colnames(comparassion_results))){
  col_name<-colnames(comparassion_results)[j]
  comparassion_results[[col_name]]<-as.numeric(comparassion_results[[col_name]])
}

View(comparassion_results)

#Paso 6 del dataframe generado buscar las 10 primeras secuencias 
#que tienen un bitscore más alto(significa que tienen mayor match)

#ordena de mayor a menor y selecciona las primeras 10 secuencias
top_indices<- order(-comparassion_results$Score)[1:10]
top_indices

#Filta las filas correspondientes, dando la informacion de las primeras 10
top_sequences <- comparassion_results[top_indices, ]   
top_sequences

#paso 7 guardar las secuecnias en nuestra caretas de resultados 
write.csv(top_sequences,file = "top_sequences.csv",row.names = FALSE)
