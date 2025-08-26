Answers to questions path, please run full script succesfully first to run those answers.

```{r}
# ---- Answers for q1..q5 ----
# Requires objects 'counts', 'dds', and 'res' created above

# q1: lanes concatenated for young -dox
# If you later generate a manifest upstream, read it here instead of hardcoding
sample_lanes <- c(ym = 1L, yp = 1L, sm = 1L, sp = 1L)
ans_q1 <- unname(sample_lanes["ym"])

# q2: library size for senescent -dox
libsize <- colSums(counts)
ans_q2 <- unname(libsize["sm"])

# q3: genes with nonzero counts in senescent +dox
ans_q3 <- sum(counts[, "ym"] > 0)

# q4: genes upregulated ≥2 fold (log2FC ≥ 1) with FDR < 0.01
res_df <- as.data.frame(res)
ans_q4 <- sum(res_df$log2FoldChange >= 1 & res_df$padj < 0.01, na.rm = TRUE)

# q5: gene ranked 3th by log2 fold change among most upregulated
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
