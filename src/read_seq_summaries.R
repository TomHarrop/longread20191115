library(data.table)
library(bit64)
library(ggplot2)
library(scales)

ss_lists <- list.files('data/wga',
                       pattern = 'sequencing_summary.txt',
                       recursive = TRUE,
                       full.names = TRUE)

names(ss_lists) <- basename(dirname(ss_lists))
ss_table_list <- lapply(ss_lists, fread)
bc_results <- rbindlist(ss_table_list, idcol = "fc", fill = TRUE)
setorder(bc_results, sequence_length_template)
l50_calcs <- bc_results[passes_filtering == TRUE, .(
    sequence_length_template,
    cs = cumsum(as.numeric(sequence_length_template)),
    tl = sum(sequence_length_template)),
    by = fc]
l50_dt <- l50_calcs[cs / tl > 0.5,
                    .(xintercept = min(sequence_length_template)),
                    by = fc]

ggplot(bc_results, 
       aes(x = sequence_length_template,
           y = mean_qscore_template)) +
    theme_grey(base_size = 8,
               base_family = "Lato") +
    facet_wrap(~ fc) +
    ylab("Mean Q score") + xlab("Read length") +
    scale_fill_viridis_c(trans = log_trans(),
                         breaks = log_breaks(),
                         guide = guide_colourbar(
                             title = "Number\nof reads"
                         )) +
    scale_x_continuous(trans = log_trans(base = 4),
                       breaks = trans_breaks(function(x) log(x, 4),
                                             function(x) 4^x)) +
    geom_hex() +
    geom_vline(mapping = aes(xintercept = xintercept),
               data = l50_dt,
               linetype = 2)

rm(ss_table_list)
gc(TRUE)
saveRDS(bc_results, "output/bc_results.Rds")

sample_8 <- fread("data/sample_8/basecalled/sequencing_summary.txt")
setorder(sample_8, sequence_length_template)
sample_8[passes_filtering == TRUE, sum(sequence_length_template)]
s8_l50 <- sample_8[passes_filtering == TRUE][
    cumsum(as.numeric(sequence_length_template)) / sum(sequence_length_template) > 0.5,
    min(sequence_length_template)]

ggplot(sample_8, 
       aes(x = sequence_length_template,
           y = mean_qscore_template)) +
    theme_grey(base_size = 8,
               base_family = "Lato") +
    ylab("Mean Q score") + xlab("Read length") +
    scale_fill_viridis_c(trans = log_trans(),
                         breaks = log_breaks(),
                         guide = guide_colourbar(
                             title = "Number\nof reads"
                         )) +
    scale_x_continuous(trans = log_trans(base = 4),
                       breaks = trans_breaks(function(x) log(x, 4),
                                             function(x) 4^x)) +
    geom_hex() +
    geom_vline(xintercept = s8_l50,
               linetype = 2)
