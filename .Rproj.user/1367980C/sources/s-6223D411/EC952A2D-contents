#Create nice plots for coral aDNA
library(tidyverse)
#need mapdam tables- this is for quality 20. This will be the benchmark across the board
samplenames <- c("S17463", "S17464", "S17465", "S17466")
readcounts <- c(14564, 3833, 7844, 1586)
age <- c(1099, -4215, 951, -3737)
meta <- data.frame(samplenames, readcounts, age)

misincorp <- list()
for (i in 1:length(samplenames)){
  misincorp[[i]] <- read.table(paste("../aDNA_local/mapDamage_totalanalysis/", samplenames[i], ".default.filtered20.coral.sorted.mapDamage/misincorporation.txt", sep=''), header = T)
}
names(misincorp) <- samplenames

tidy_misincorp <- bind_rows(misincorp, .id = "column_label")
freqs_5p <- tidy_misincorp %>% filter(End=="5p", Std == "+", Pos <=15) %>% group_by(column_label) %>% mutate(ct_freq = C.T/C) %>%
  mutate(ag_freq = A.G/A, ga_freq = G.A/G, tc_freq = T.C/T, gc_freq = G.C/G, ac_freq = A.C/A, 
         at_freq = A.T/A, cg_freq = G.C/C, ca_freq= C.A/C, tg_freq = T.G/T, ta_freq=T.A/T, gt_freq = G.T/G) %>%
  mutate(firstpos = ifelse(Pos == 1, "Y", "N" ))

avg_freqs_5p <- freqs_5p %>% ungroup() %>% mutate(bkg = rowMeans(select(.,ag_freq:gt_freq)))

#all samps 

sample = "S17463"
#one of these is subset, the other is not
S17463_p <- avg_freqs_5p %>% filter(column_label==sample) %>% select(column_label, Pos, ct_freq:bkg)
S17463_p <- avg_freqs_5p %>% group_by(column_label) %>% select(column_label, Pos, ct_freq:bkg)

S17463_p_tidy <- S17463_p %>% pivot_longer(cols = ct_freq:gt_freq, names_to = "mut", values_to = "freq") %>% mutate(ct_mut = ifelse(mut == "ct_freq", "Y", "N")) %>%
  group_by(column_label, ct_mut, Pos) %>% summarise(N = n(), avg.freq = mean(freq), SE.low = avg.freq - (sd(freq)/sqrt(N)), SE.high =avg.freq + (sd(freq)/sqrt(N))) %>%
  mutate(age = ifelse((column_label %in% c("S17463", "S17465")), "Newer", "Older")) 


ggplot(filter(S17463_p_tidy, ct_mut == "Y"), aes(x = Pos, y = avg.freq, linetype = age, colour=column_label)) +
  geom_line() +
  geom_ribbon(aes(ymin=SE.low, ymax=SE.high), lty = 2, alpha = 0.2) +
  theme_minimal()


ggplot(S17463_p_tidy, aes(x=Pos, y = avg.freq)) + geom_line() + geom_ribbon(aes(ymin=SE.low, ymax=SE.high), lty = 2, alpha = 0.2) +
  theme_minimal()

#########=========Read in zoox damage patterns for 14763######

zooxdam63 <- read.table("../aDNA_local/overall.modern_meta/S17465.c.filtered20.sorted.mapDamage/misincorporation.txt", header = T)

freqs_5p <- zooxdam63 %>% filter(End=="5p", Std == "+", Pos <=15)  %>% mutate(ct_freq = C.T/C) %>%
  mutate(ag_freq = A.G/A, ga_freq = G.A/G, tc_freq = T.C/T, gc_freq = G.C/G, ac_freq = A.C/A, 
         at_freq = A.T/A, cg_freq = G.C/C, ca_freq= C.A/C, tg_freq = T.G/T, ta_freq=T.A/T, gt_freq = G.T/G) %>%
  mutate(firstpos = ifelse(Pos == 1, "Y", "N" ))

zoox_freqs_5p <- freqs_5p %>% ungroup() %>% mutate(bkg = rowMeans(select(.,ag_freq:gt_freq))) %>% select(Pos, ct_freq:bkg) %>% pivot_longer(cols = ct_freq:gt_freq, names_to = "mut", values_to = "freq") %>% mutate(ct_mut = ifelse(mut == "ct_freq", "Y", "N")) %>%
  group_by(ct_mut, Pos) %>% summarise(N = n(), avg.freq = mean(freq), SE.low = avg.freq - (sd(freq)/sqrt(N)), SE.high =avg.freq + (sd(freq)/sqrt(N)))
#Test read counts:

metagenome_names <- read.table("filtered20files")
metagenome_counts <- read.table("filtered20counts")
metagenome_df <- data.frame(metagenome_names, metagenome_counts)
metagenome_df$V1 <- as.character(metagenome_df$V1)
colnames(metagenome_df) <- c("names", "counts")
samples <- c("S17463", "S17464")
likely_ancient <- list()
for (i in 1:length(samples)){
  likely_ancient[[i]] <- read.table(paste(samples[i], ".likely_ancient", sep=''), header = F, stringsAsFactors = F)
}
ancient <- bind_rows(likely_ancient, .id = "column_label")

metagenome_df <- metagenome_df %>% mutate(anc = ifelse(names %in% ancient$V1, "Y", "N")) %>% separate(names, c("S", "Meta", "sort", "filter", "type1", "type2")) %>%
  select(S, type1, type2, counts, anc) %>% unite("align", type1:type2, sep = "_") %>% mutate()

metagenome_sub <- metagenome_df %>% filter(align != "NA_NA", S %in% samples) %>% arrange(desc(anc))
new_order <- unique(metagenome_sub$align)
metagenome_sub$align <- factor(metagenome_sub$align, levels=new_order)
metagenome_byanc <- metagenome_sub %>% group_by(S, anc) %>% summarise(sum = sum(counts))
library(ggplot2)
# Barplot

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

bp<- ggplot(metagenome_sub, aes(x="", y=counts, fill=align, alpha=anc))+
  geom_bar(width = 1, stat = "identity", position = "fill") + scale_alpha_discrete(range = c(0.35, 1)) +
  theme_minimal() + #coord_polar("y", start = 0) + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    plot.title=element_text(size=14, face="bold")
  ) + theme(legend.position = "none") + facet_grid(facets=. ~ S)
bp

metagenome_onlyanc <- metagenome_sub %>% filter(anc == "Y")
bp<- ggplot(metagenome_onlyanc, aes(x=S, y=counts, fill=align))+
  geom_bar(width = 1, stat = "identity", position = "fill")
bp



bp<- ggplot(metagenome_byanc, aes(x=S, y=sum, fill=anc))+
  geom_bar(width = 1, stat = "identity", position = "fill")
bp

#Get broad level names and counts to make pie charts
master_names <- read.table("master_names")
master_counts <- read.table("master_counts")
master_df <- data.frame(master_names, master_counts)
colnames(master_df) <- c("names", "counts")

master <- master_df %>% separate(names, c("S", "align", "filter", "bam")) %>% select(S, align, counts) %>% mutate(anc = ifelse(align == "coral" | (align == "meta" & S == "S17463"), "Y", "N")) %>%
  arrange(anc)

new_order <- unique(master$align)
master$align <- factor(master$align, levels=new_order)
just1 <- master %>% filter(S == "S17466")
all <- master %>% filter(anc == "Y")
justanc <- just1 %>% filter(anc == "Y")
notanc <- just1 %>% filter(anc == "N", align %in% c("a", "b", "c", "d"))

master_bar<- ggplot(data = just1, aes(x="", y=counts, fill=align))+
  geom_bar(width = .5, stat = "identity", position = "fill") + scale_alpha_discrete(range = c(0.45, 1)) +
  #geom_bar(data = metagenome_byanc, aes(x="", y=sum)) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    plot.title=element_text(size=14, face="bold")
  ) + facet_grid(facets=. ~ S) + scale_fill_manual(values=unlist(pal_diverse, use.names = FALSE)) 
master_bar



master_pie<- ggplot(data = just1, aes(x="", y=counts, fill=align))+
  geom_bar(width = 1, stat = "identity", position = "fill") + scale_alpha_discrete(range = c(0.35, 1)) +
  theme_minimal() + coord_polar("y", start = 0) + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    plot.title=element_text(size=14, face="bold")
  ) + theme(legend.position = "none") + facet_grid(facets=. ~ S)
master_pie


S17463<- ggplot(data = just1, aes(x="", y=counts, fill=align))+
  geom_bar(width = 1, stat = "identity", position = "fill") + scale_alpha_discrete(range = c(0.45, 1)) +
  #geom_bar(data = metagenome_byanc, aes(x="", y=sum)) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    plot.title=element_text(size=14, face="bold")
  ) + facet_grid(facets=. ~ S) + scale_fill_manual(values=unlist(pal_diverse, use.names = FALSE)) 
S17463

S17463_anc<- ggplot(data = justanc, aes(x="", y=counts, fill=align))+
  geom_bar(width = 1, stat = "identity", position = "fill") + scale_alpha_discrete(range = c(0.45, 1)) +
  #geom_bar(data = metagenome_byanc, aes(x="", y=sum)) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    plot.title=element_text(size=14, face="bold")
  ) + facet_grid(facets=. ~ S) + scale_fill_manual(values=unlist(c(pal_diverse$meta, pal_diverse$coral), use.names = FALSE)) 
S17463_anc

S17463_noanc<- ggplot(data = notanc, aes(x="", y=counts, fill=align))+
  geom_bar(width = 1, stat = "identity", position = "fill") + scale_alpha_discrete(range = c(0.45, 1)) +
  #geom_bar(data = metagenome_byanc, aes(x="", y=sum)) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    plot.title=element_text(size=14, face="bold")
  ) + facet_grid(facets=. ~ S) + scale_fill_manual(values=unlist(pal_diverse[c("zoox1", "zoox2", "zoox3", "zoox4")], use.names = FALSE)) 
S17463_noanc


#From the very top of the code 
ggplot(filter(S17463_p_tidy, ct_mut == "N"), aes(x = Pos, y = avg.freq)) + 
  geom_line() + geom_ribbon(aes(ymin=SE.low, ymax=SE.high), lty = 2, alpha = 0.2) + ylim(0, 0.15) +
  geom_line(filter(S17463_p_tidy, ct_mut == "Y"), mapping = aes(color = unlist(pal_diverse["coral"], use.names = F))) + 
  theme_minimal() + theme(legend.position = "none") 

ggplot(filter(zoox_freqs_5p, ct_mut == "N"), aes(x = Pos, y = avg.freq)) + 
  geom_line() + geom_ribbon(aes(ymin=SE.low, ymax=SE.high), lty = 2, alpha = 0.2) + ylim(0, 0.4) +
  geom_line(filter(zoox_freqs_5p, ct_mut == "Y"), mapping = aes(colour = unlist(pal_diverse["zoox3"], use.names = F))) + 
  theme_minimal() + theme(legend.position = "none") 

pal_diverse <- list(
  "zoox1" = "#16BAC5", #bluegreen
  "zoox2" = "#274690", #darkcornflowerblue
  "zoox3" = "#5FBFF9", #mayablue
  "zoox4" = "#81F499", #lightgreen
  "cca" = "#DAD4EF", #languidlav
  "redalg" = "#7B506F", #lavender
  "meta" = "#BDC696", #sage
  "coral" = "#EF8275" #salmon 
  )

nice_colors <- list(
  "lavender" = "#7B506F",
  "rubyred" = "#A20021"
)


#drabolivegreen for special meta
