#Combinar os dados de DEPS com a tabela de dados normalizados para todos
# os pacientes

# Usando a funcao merge

# D0

D0_merge <- merge(D0.Control, DEPs_2, by="Accession", all=T)

write.table(D0_merge, file = "D0_merge.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

# D7

D7_merge <- merge(D7.Control, DEPs_2, by="Accession", all=T)

write.table(D7_merge, file = "D7_merge.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

# cs_31

CS31_merge <- merge(D30.Control, DEPs_2, by="Accession", all=T)

write.table(CS31_merge, file = "CS31_merge.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
