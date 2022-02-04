library(STRINGdb)
library(igraph)
library(linkcomm)

E <- read.table(file = 'data/E.tsv', sep = '\t', header = TRUE)
M <- read.table(file = 'data/M.tsv', sep = '\t', header = TRUE)
N <- read.table(file = 'data/N.tsv', sep = '\t', header = TRUE)
S <- read.table(file = 'data/Spike.tsv', sep = '\t', header = TRUE)
nsp1 <- read.table(file = 'data/nsp1.tsv', sep = '\t', header = TRUE)
nsp2 <- read.table(file = 'data/nsp2.tsv', sep = '\t', header = TRUE)
nsp3 <- read.table(file = 'data/nsp3.tsv', sep = '\t', header = TRUE)
nsp4 <- read.table(file = 'data/nsp4.tsv', sep = '\t', header = TRUE)
nsp5 <- read.table(file = 'data/nsp5.tsv', sep = '\t', header = TRUE)
nsp6 <- read.table(file = 'data/nsp6.tsv', sep = '\t', header = TRUE)
nsp7 <- read.table(file = 'data/nsp7.tsv', sep = '\t', header = TRUE)
nsp8 <- read.table(file = 'data/nsp8.tsv', sep = '\t', header = TRUE)
nsp9 <- read.table(file = 'data/nsp9.tsv', sep = '\t', header = TRUE)
nsp10 <- read.table(file = 'data/nsp10.tsv', sep = '\t', header = TRUE)
nsp11 <- read.table(file = 'data/nsp11.tsv', sep = '\t', header = TRUE)
nsp12 <- read.table(file = 'data/nsp12.tsv', sep = '\t', header = TRUE)
nsp13 <- read.table(file = 'data/nsp13.tsv', sep = '\t', header = TRUE)
orf3a <- read.table(file = 'data/orf3a.tsv', sep = '\t', header = TRUE)
orf3b <- read.table(file = 'data/orf3b.tsv', sep = '\t', header = TRUE)
orf6 <- read.table(file = 'data/orf6.tsv', sep = '\t', header = TRUE)
orf8 <- read.table(file = 'data/orf8.tsv', sep = '\t', header = TRUE)
orf9b <- read.table(file = 'data/orf9b.tsv', sep = '\t', header = TRUE)
orf9c <- read.table(file = 'data/orf9c.tsv', sep = '\t', header = TRUE)
orf10 <- read.table(file = 'data/orf10.tsv', sep = '\t', header = TRUE)

proteins_relations <- rbind(E, rbind(M, rbind(N, rbind(S, rbind(nsp1, rbind(nsp2, rbind(nsp3, rbind(nsp4, rbind(nsp5, rbind(nsp6, rbind(nsp7, rbind(nsp8, rbind(nsp9, rbind(nsp10, rbind(nsp11, rbind(nsp12, rbind(nsp13, rbind(orf3a, rbind(orf3b, rbind(orf6, rbind(orf8, rbind(orf9b, rbind(orf9c, orf10)))))))))))))))))))))))

# New STRING_db object
string_db <- STRINGdb$new(version="11",
                          species=9606,
                          score_threshold=995,
                          input_directory="")

# Map to STRING
proteins_mapped = string_db$map(proteins_relations, "node1", removeUnmappedRows = TRUE)

# Plot the STRING network 
png(file = "../results/figuraSTRINGdb.png")

string_db$plot_network(proteins_mapped$STRING_id)

dev.off()

nodes_diseases <- read.table(file = 'data/nodes&diseases.tsv', fill=TRUE, header=TRUE, quote="", sep="\t", encoding="UTF-8")
subset_nodes_diseases <- subset(nodes_diseases, disease!="")

subset_mat <- as.matrix(subset_nodes_diseases[,1:2])

png(file = "../results/figuraiGraph.png")

plot(graph.edgelist(subset_mat, directed = FALSE))


dev.off()


#### Link Communities ######

ppi = graph_from_data_frame(human_pp,directed = F)
subset_nodes_diseases_graph = graph_from_data_frame(subset_nodes_diseases,directed = F)

png(file = "../results/figuraLinkcomm.png")

plot(subset_nodes_diseases_graph,
     vertex.color = "tomato",
     vertex.size = degree(ppi)/5,
     vertex.label.color = "black",
     vertex.label.family = "Helvetica",
     vertex.label.cex = 0.6,
     layout = layout.kamada.kawai
)

dev.off()
