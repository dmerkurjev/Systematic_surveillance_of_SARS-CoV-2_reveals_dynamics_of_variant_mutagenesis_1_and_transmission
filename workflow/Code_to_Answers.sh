Answers to questions path, please run full script succesfully first to run those answers.

```{r}
# ---- Answers for q1..q5 ----
# Requires objects 'counts', 'dds', and 'res' created above

# q1: How many sequencing lanes were concatenated to form sample normal sample 2?
# If you later generate a manifest upstream, read it here instead of hardcoding
sample_lanes <- c(ym = 1L, yp = 1L, sm = 1L, sp = 1L)
ans_q1 <- unname(sample_lanes["norm2"])

# q2: What is the library size (total read counts) for normal sample 1 and tumor sample 1?
libsize <- colSums(counts)
ans_q2_1 <- unname(libsize["norm1"])
ans_q2_2 <- unname(libsize["tum1"])
ans_q2=ans_q2_1+ans_q2_2

# q3: How many genes have nonzero counts in sample tumor 1?
ans_q3 <- sum(counts[, "norm1"] > 0)

# q4: How many genes are upregulated \uc0\u8805 2-fold (log2FC \u8805  1) in  the TCGA ccRCC cohort compared to normal human kidney  with p-value < 0.001?
res_df <- as.data.frame(res)
ans_q4 <- sum(res_df$log2FoldChange >= 1 & res_df$padj < 0.001, na.rm = TRUE)

# q5: How many genes are shared (upregulated) between species?
res_df$SYMBOL <- rownames(res_df)
up_rank <- res_df %>%
  dplyr::filter(!is.na(log2FoldChange) & log2FoldChange > 0) %>%
  dplyr::arrange(dplyr::desc(log2FoldChange), padj, pvalue)
ans_q5 <- if (nrow(up_rank) >= 3) up_rank$SYMBOL[3] else NA_character_

answers <- tibble::tibble(
  id = c("q1","q2","q3","q4","q5"),
  answer = c(ans_q1, ans_q2, ans_q3, ans_q4, ans_q5)
)

answers
```
