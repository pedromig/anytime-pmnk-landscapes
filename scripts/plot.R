#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 4) {
  pls_file <- args[1]
  gsemo_file <- args[2]
  b_ibea_file <- args[3]
  a_ibea_file <- args[4]
} else {
  stop(sprintf("Usage: ./plot.R {PLS} {GSEMO} {B-IBEA} {A-IBEA}"))
}

df_pls <- read.csv(
  pls_file,
  head = T,
  col.names = c("Evaluation", "Hypervolume")
)
df_pls$algorithm <- "PLS"

df_gsemo <- read.csv(
  gsemo_file,
  head = T,
  col.names = c("Evaluation", "Hypervolume")
)
df_gsemo$algorithm <- "GSEMO"

df_b_ibea <- read.csv(
  b_ibea_file,
  head = T,
  col.names = c("Evaluation", "Generation", "Hypervolume")
)
df_b_ibea$algorithm <- "B-IBEA"
df_b_ibea <- subset(df_b_ibea, select = -Generation)

df_a_ibea <- read.csv(
  a_ibea_file,
  head = T,
  col.names = c("Evaluation", "Generation", "Hypervolume")
)
df_a_ibea$algorithm <- "A-IBEA"
df_a_ibea <- subset(df_a_ibea, select = -Generation)

df <- rbind(df_pls, df_gsemo, df_b_ibea, df_a_ibea)
df$evaluation <- as.numeric(df$Evaluation)
df$hypervolume <- as.numeric(df$Hypervolume)
df$algorithm <- as.factor(df$algorithm)

ggplot(df, aes(evaluation, hypervolume, color = algorithm)) +
  geom_step(direction = "hv")
