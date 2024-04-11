---
title: "Sim-dataclean"
author: "Jay Meyer"
date: "2024-02-29"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)

data <- read_csv("~/Documents/python/mrp_params.csv") %>%
  select(from_desc,to_desc,conditions,prob,mean_time_diff)

trans_probs <- data %>%
  select(-mean_time_diff) %>%
  pivot_wider(names_from = to_desc, values_from = prob)

margin_probs <- data %>%
  select(-prob) %>%
  pivot_wider(names_from=to_desc, values_from = mean_time_diff)

write_csv(trans_probs, "~/Documents/python/mrp_trans_probs.csv")
write_csv(margin_probs, "~/Documents/python/mrp_mtds.csv")

```



