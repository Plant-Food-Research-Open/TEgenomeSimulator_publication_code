# Convert denovoTE-eval output gff to the format required for the *_sim_genome functions 
convert_denovote_gff <- function(denovote_gff){
  newgff <- denovote_gff |>
    dplyr::mutate(ID = gsub(";.*", "", V9)) |>
    dplyr::mutate(Class = paste0("Classification=", gsub(".*#", "", V9))) |>
    dplyr::mutate(Class = gsub(";.*", "", Class)) |>
    dplyr::mutate(Class = gsub("_n.*", "", Class)) |>
    dplyr::mutate(Class = gsub("_p.*", "", Class)) |>
    dplyr::mutate(Int = gsub(".*fragment=", "", V9)) |>
    dplyr::mutate(Int = gsub(";.*", "", Int)) |>
    dplyr::mutate(Int = gsub(";.*", "", Int)) |>
    dplyr::mutate(Int = as.numeric(Int)/100) |>
    dplyr::mutate(Int = case_when(is.na(Int) ~ 1, TRUE ~ Int)) |>
    dplyr::mutate(Idn = gsub(".*identity=", "", V9)) |>
    dplyr::mutate(Idn = gsub(";.*", "", Idn)) |>
    dplyr::mutate(Idn = as.numeric(Idn)/100) |>
    dplyr::mutate(V9 = paste0(ID, ";Integrity=", Int, ";Identity=", Idn, ";", Class)) |> 
    dplyr::mutate(V1 = "chr1", V2 = "denovoTEeval", V3 = gsub("Classification=", "", Class)) |>
    dplyr::mutate(V3 = gsub("\\/.*", "", V3)) |>
    dplyr::select(-c(ID, Class, Int, Idn))
}

# This function summarise numbers of annotated TE loci in the simulated genome grouped by super family
te_by_superfamily_sim_genome <- function(gff){
  break_gff <- gff |> 
    dplyr::mutate(superfamily = gsub(".*Classification=", "", V9)) |>
    dplyr::mutate(superfamily = gsub(";.*", "", superfamily)) |>
    dplyr::mutate(superfamily = gsub("auto", "", superfamily)) |>
    dplyr::mutate(superfamily = gsub("nona", "", superfamily)) |>
    dplyr::mutate(superfamily = gsub("DNA\\/Helitron", "RC\\/Helitron", superfamily)) |>
    dplyr::mutate(superfamily = gsub("DNA\\/CMC.*", "DNA\\/CACTA", superfamily)) |>
    dplyr::mutate(superfamily = gsub("DNA\\/MULE.*", "DNA\\/MULE", superfamily)) |>
    dplyr::mutate(superfamily = gsub("DNA\\/MuDR.*", "DNA\\/MULE", superfamily)) |>
    dplyr::mutate(superfamily = gsub("DNA\\/PIF-Harbinger", "DNA\\/Harbinger", superfamily)) |>
    dplyr::mutate(superfamily = gsub("DNA\\/TcMar.*", "DNA\\/Mariner", superfamily)) |>
    dplyr::mutate(superfamily = gsub("DNA\\/hAT.*", "DNA\\/hAT", superfamily)) |>
    dplyr::mutate(superfamily = gsub("LTR\\/Gypsy", "LTR\\/Ty3", superfamily)) |>
    dplyr::mutate(superfamily = gsub("DNA$", "DNA\\/unknown", superfamily)) |>
    dplyr::mutate(TE_family = gsub(".*ID=", "", V9)) |>
    dplyr::mutate(TE_family = gsub("\\#.*", "", TE_family)) |>
    dplyr::mutate(TE_family = gsub("_NC_003070.*", "", TE_family))

  a_sum <- break_gff |> group_by(superfamily) |> summarise(count = n())
  b_sum <- break_gff |> group_by(superfamily) |> dplyr::mutate(length = V5 - V4 + 1) |> summarise(bp = sum(length))
  c_sum <- break_gff |> group_by(superfamily) |> distinct(superfamily, TE_family) |> summarise(count = n())
  
  data_sum <- a_sum
  data_sum$bp <- b_sum$bp
  data_sum$family_count <- c_sum$count
  colnames(data_sum)[1] <- "TE"

  return(data_sum)
}



# This function breakdown the information in column 9 of the gff file
te_breakdown_sim_genome <- function(gff){
  break_gff <- gff |> 
    dplyr::mutate(sfam = gsub(".*Classification=", "", V9)) |>
    dplyr::mutate(sfam = gsub(";.*", "", sfam)) |>
    dplyr::mutate(sfam = gsub("auto", "", sfam)) |>
    dplyr::mutate(sfam = gsub("nona", "", sfam)) |>
    dplyr::mutate(sfam = gsub("DNA\\/Helitron", "RC\\/Helitron", sfam)) |>
    dplyr::mutate(sfam = gsub("DNA\\/CMC.*", "DNA\\/CACTA", sfam)) |>
    dplyr::mutate(sfam = gsub("DNA\\/MULE.*", "DNA\\/MULE", sfam)) |>
    dplyr::mutate(sfam = gsub("DNA\\/MuDR.*", "DNA\\/MULE", sfam)) |>
    dplyr::mutate(sfam = gsub("DNA\\/PIF-Harbinger", "DNA\\/Harbinger", sfam)) |>
    dplyr::mutate(sfam = gsub("DNA\\/TcMar.*", "DNA\\/Mariner", sfam)) |>
    dplyr::mutate(sfam = gsub("DNA\\/hAT.*", "DNA\\/hAT", sfam)) |>
    dplyr::mutate(sfam = gsub("LTR\\/Gypsy", "LTR\\/Ty3", sfam)) |>
    dplyr::mutate(sfam = gsub("DNA$", "DNA\\/unknown", sfam)) |>
    dplyr::mutate(idn = gsub(".*Identity=", "", V9)) |>
    dplyr::mutate(idn = gsub(";.*", "", idn)) |>
    dplyr::mutate(itg = gsub(".*Integrity=", "", V9)) |>
    dplyr::mutate(itg = gsub(";.*", "", itg)) |>
    dplyr::mutate(idn = as.numeric(idn)) |>
    dplyr::mutate(itg = as.numeric(itg)) |>
    dplyr::mutate(div = 1 - as.numeric(idn)) |>
    dplyr::mutate(family = gsub(".*ID=", "", V9)) |>
    dplyr::mutate(family = gsub("\\#.*", "", family))

  nest <- break_gff[grepl("Nest_in", break_gff$V9),]
  notnest <- break_gff[!grepl("Nest_in", break_gff$V9),]
  nest <- nest |> dplyr::mutate(Nested_in = gsub(".*Nest_in=", "", V9)) |> dplyr::mutate(Nested_in = gsub(";.*", "", Nested_in))
  notnest <- notnest |> dplyr::mutate(Nested_in = "NA")

  break_gff <- rbind(nest, notnest)

  cut <- break_gff[grepl("Cut_by", break_gff$V9),]
  notcut <- break_gff[!grepl("Cut_by", break_gff$V9),]
  cut <- cut |> dplyr::mutate(Cut_at = gsub(".*Cut_at=", "", V9)) |>
    dplyr::mutate(Cut_at = gsub(";.*", "", Cut_at)) |>
    dplyr::mutate(Cut_by = gsub(".*Cut_by=", "", V9)) |>
    dplyr::mutate(Cut_by = gsub(";.*", "", Cut_by))
  notcut <- notcut |> dplyr::mutate(Cut_at = "NA", Cut_by = "NA")

  break_gff <- rbind(cut, notcut)

  return(break_gff)
}



# this function generate a barplot to show TE loci grouped by superfamily
plot_te_loci_barplot <- function(data_sum, genomename, ylim=8e+05){
    data.melted <- reshape2::melt(data_sum, id = "TE_rename")
    colnames(data.melted) <- c("TE", "Genome", "Count")

    # set up color pallet
    blue <- rev(brewer.pal(n = 9, "YlGnBu"))
    bgreen <- rev(brewer.pal(n=9, "PuBuGn"))
    purple <- rev(brewer.pal(n = 9, "Purples"))
    pink <- rev(brewer.pal(n = 9, "RdPu"))
    orange <- rev(brewer.pal(n = 9, "YlOrBr"))
    mycolor <- c(orange[5:7], purple[2], blue[2:6], pink[5], bgreen[2])

    # generate plot
    plot <- ggplot(data.melted, aes(x = Genome, y = Count, fill = TE)) +
    geom_bar(stat = "identity", color="white", position = "stack") +
    scale_fill_manual(values=mycolor) +
    scale_y_continuous(limits = c(0, ylim),
                       breaks = seq(0, 8e+05, by = 2e+05)) +
    scale_x_discrete(labels= genomename) +
    ggtitle("Annotated TE loci") + ylab ("Count of TE loci") +
    theme(title = element_text(size = 16),
          # plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = rel(1.3)),
          axis.text.x = element_text(angle = 45, hjust = 1, size = rel(2)),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = rel(1)),
          # legend.title = element_blank(),
          legend.text = element_text(size = rel(1)),
          legend.key.size = unit(0.5, "cm"),
          # legend.position = c(0.3, 0.9),
          # legend.background = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"))

    return(plot)
}

# this function generate a piechart to show the proportion of the genome occupied by TE and non-TE
plot_sim_genome_bp_piechart <- function(total_genome_bp, total_te_bp, gname){
    data_sum <- data.frame(group = c('non-TE', 'TE'),
                           bp = c(total_genome_bp - total_te_bp, total_te_bp))
    data_sum <- data_sum |> dplyr::mutate(cumulative = cumsum(bp),
                                    percentage = round(100*(bp/total_genome_bp), 2),
                                    midpoint = (cumulative - bp/2), label = paste0(percentage, "%"))

    # set up color pallet
    mycolor <- c('grey70', 'black')

    # generate plot
    plot <- ggplot(data_sum, aes(x = "", y = bp, fill = group)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 1, color = 'white', linewidth = 2, alpha = 0.8) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values=mycolor) +
      annotate(geom = "text", x = 1, y = data_sum |> filter(percentage > 1) |> pull(midpoint), hjust = 0.5,
               label = data_sum |> filter(percentage > 1) |> pull(label) , size = rel(4), color = "white", fontface = "bold") +
      ggtitle(paste0("Genome: ", gname, "\n", "Genome size: ", format(total_genome_bp, big.mark=",",scientific=FALSE), " bp")) +
      theme(title = element_text(size = 8),
            #plot.title = element_text(hjust = 0.5),
            axis.text = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(size = rel(1)),
            legend.key.size = unit(0.5, "cm"),
            # legend.position = c(0.3, 0.9),
            legend.background = element_blank(),
            panel.background = element_rect(fill = "white", color = "white"))

    return(plot)
}


# this function generate a piechart to show TE loci of the genome grouped by superfamily
plot_te_loci_piechart <- function(data.sum, genomename){
  total = sum(data.sum$count)
  data.sum$total <- rep(total, nrow(data.sum))

  # order_te_superfamilies()/te_superfamily_colors() give the same fixed, locale-
  # independent order and colour used by plot_te_superfamily_bars_meansd() - no digit
  # prefix baked into the label; TE is releveled and the rows re-arranged to match so
  # the manually-computed cumulative/midpoint stay in sync with the actual stack order
  te_order <- order_te_superfamilies(data.sum$TE)
  te_color_map <- te_superfamily_colors(data.sum$TE)

  data.sum <- data.sum |>
    dplyr::mutate(TE = factor(as.character(TE), levels = te_order)) |>
    dplyr::arrange(TE) |>
    dplyr::mutate(cumulative = cumsum(count),
                  percentage = round(100*(count/total), 1),
                  midpoint = (cumulative - count/2),
                  label = paste0(percentage, "%"))

  # generate plot
  plot <- ggplot(data.sum, aes(x = "", y = count, fill = TE)) +
    geom_bar(stat = "identity", color = 'white', position = position_stack(reverse = TRUE), width = 1, linewidth = 0.5, alpha = 0.8) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = te_color_map) +
    annotate(geom = "text", x = 1.2, y = data.sum |> filter(percentage > 5) |> pull(midpoint), hjust = 0.5, vjust = 0.7,
             label = data.sum |> filter(percentage > 5) |> pull(label) , size = rel(4.5), color = "white", fontface = "bold") +
    ggtitle(paste0("Genome: ", genomename,
                   "\nTotal TE loci: ", format(total, big.mark=",", scientific=FALSE))) +
    theme(title = element_text(size = 8),
          axis.text = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = rel(1)),
          legend.key.size = unit(0.5, "cm"),
          legend.background = element_blank(),
          panel.background = element_rect(fill = "white", color = "white"))

  return(plot)
}


# this function generate a piechart to show simulated TE sequence proportion grouped by superfamily
plot_te_bp_piechart <- function(data.sum, genomename){
    total = sum(data.sum$bp)
    data.sum$total <- rep(total, nrow(data.sum))

    # see plot_te_loci_piechart() - same fixed, locale-independent order/colour source,
    # no digit prefix on the label
    te_order <- order_te_superfamilies(data.sum$TE)
    te_color_map <- te_superfamily_colors(data.sum$TE)

    data.sum <- data.sum |>
      dplyr::mutate(TE = factor(as.character(TE), levels = te_order)) |>
      dplyr::arrange(TE) |>
      dplyr::mutate(cumulative = cumsum(bp), percentage = round(100*(bp/total), 1), midpoint = (cumulative - bp/2), label = paste0(percentage, "%"))

    # generate plot
    plot <- ggplot(data.sum, aes(x = "", y = bp, fill = TE)) +
    geom_bar(stat = "identity", color = 'white', position = position_stack(reverse = TRUE), width = 1, linewidth = 0.5, alpha = 0.8) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = te_color_map) +
    annotate(geom = "text", x = 1.2, y = data.sum |> filter(percentage > 5) |> pull(midpoint), hjust = 0.5,
             label = data.sum |> filter(percentage > 5) |> pull(label) , size = rel(4.5), color = "white", fontface = "bold") +
    ggtitle(paste0("Genome: ", genomename,
                   "\ntotal TE bases: ", format(total, big.mark=",", scientific=FALSE), " bp")) +
    theme(title = element_text(size = 8),
          axis.text = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = rel(1)),
          legend.key.size = unit(0.5, "cm"),
          legend.background = element_blank(),
          panel.background = element_rect(fill = "white", color = "white"))

  return(plot)
}




# Plot TE loci divergence of the sim genome; y-axis = TE loci count
plot_sim_divergence_by_loci <- function(gff_breakdown, ymax, ybreak, genomename){
  # same fixed, locale-independent order/colour as the pie charts and bar charts
  te_color_map <- te_superfamily_colors(gff_breakdown$sfam)
  gff_breakdown <- gff_breakdown |> dplyr::mutate(sfam = factor(as.character(sfam), levels = order_te_superfamilies(sfam)))

  plot <- ggplot(gff_breakdown, aes(as.numeric(div), fill = sfam)) +
    #call geom_histogram with position="dodge" to offset the bars and manual binwidth of 2
    geom_histogram(position = "stack", binwidth = 0.02, color = "white", alpha = 0.8) +
    scale_fill_manual(values = te_color_map) +
    ggtitle(paste0("Genome: ", genomename)) +
    xlab("Divergence (p-dist)") + ylab("TE loci") +
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, by = ybreak)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    coord_cartesian(xlim = c(0, 0.6)) +  # <-- changed from scale_x_continuous(limits=...
    theme(title = element_text(size = 8),
          axis.text.y = element_text(size = rel(1.8)),
          axis.text.x = element_text(size = rel(1.8)),
          axis.title.x = element_text(size = rel(1.8)),
          axis.title.y = element_text(size = rel(1.8)),
          legend.title = element_blank(),
          legend.text = element_text(size = rel(1)),
          legend.key.size = unit(0.4, "cm"),
          legend.background = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"))
  return(plot)
}


# Plot TE loci integrity of the sim genome; y-axis = TE loci count
plot_sim_integrity_by_loci <- function(gff_breakdown, ymax, ybreak, genomename){
  # same fixed, locale-independent order/colour as the pie charts and bar charts
  te_color_map <- te_superfamily_colors(gff_breakdown$sfam)
  gff_breakdown <- gff_breakdown |> dplyr::mutate(sfam = factor(as.character(sfam), levels = order_te_superfamilies(sfam)))

  plot <- ggplot(gff_breakdown, aes(as.numeric(itg), fill = sfam)) +
    #call geom_histogram with position="dodge" to offset the bars and manual binwidth of 2
    geom_histogram(position = "stack", binwidth = 0.02, color = "white", alpha = 0.8) +
    scale_fill_manual(values = te_color_map) +
    ggtitle(paste0("Genome: ", genomename)) +
    xlab("Integrity") + ylab("TE loci") +
    scale_y_continuous(breaks = seq(0, ymax, by = ybreak)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, ymax)) +  # <-- changed from scale_x_continuous(limits=... and scale_y_continuous(limits=..
    theme(title = element_text(size = 8),
          axis.text.y = element_text(size = rel(1.8)),
          axis.text.x = element_text(size = rel(1.8)),
          axis.title.x = element_text(size = rel(1.8)),
          axis.title.y = element_text(size = rel(1.8)),
          legend.title = element_blank(),
          legend.text = element_text(size = rel(1)),
          legend.key.size = unit(0.4, "cm"),
          legend.background = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"))
  return(plot)
}

# Plot TE loci divergence of the original genome; y-axis = TE loci count
plot_ori_divergence_by_loci <- function(processed_rm_out, ymax, ybreak, genomename){
  # set up color pallet
  mycolor <- c(rep(c('#CC9999', '#CC6666', '#993344', '#996699', '#996633',
                     '#666699', '#336699', '#0099AA', '#99CCCC', '#006655', 
                     '#669966', '#669922', '#99CC99', '#999966', '#CCCC99',
                     '#CCCC77', '#CCCC11', '#337799', '#336666', '#004499', 
                     '#333366', '#3399CC', '#99CC66', '#CC9900', '#CC6633',
                     '#663333', '#CC3333', '#993399', '#CC99CC', '#CCCCCC'),2))
  
  plot <- ggplot(processed_rm_out, aes(as.numeric(div), fill = superfamily)) +
    #call geom_histogram with position="dodge" to offset the bars and manual binwidth of 2
    geom_histogram(position = "stack", binwidth = 0.02, color = "white", alpha = 0.8) +
    scale_fill_manual(values=mycolor) +
    ggtitle(paste0("Genome: ", genomename)) +
    xlab("Divergence") + ylab("TE loci") +
    scale_y_continuous(breaks = seq(0, ymax, by = ybreak)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    coord_cartesian(xlim = c(0, 0.6), ylim = c(0, ymax)) +  # <-- changed from scale_x_continuous(limits=... and scale_y_continuous(limits=..
    theme(title = element_text(size = 8),
          axis.text.y = element_text(size = rel(1.8)),
          axis.text.x = element_text(size = rel(1.8)),
          axis.title.x = element_text(size = rel(1.8)),
          axis.title.y = element_text(size = rel(1.8)),
          legend.title = element_blank(),
          legend.text = element_text(size = rel(1)),
          legend.key.size = unit(0.4, "cm"),
          legend.background = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"))
  return(plot)
}



# Plot TE loci integrity of the original genome; y-axis = TE loci count
plot_ori_integrity_by_loci <- function(processed_rm_out, ymax, ybreak, genomename){
  # set up color pallet
  mycolor <- c(rep(c('#CC9999', '#CC6666', '#993344', '#996699', '#996633',
                     '#666699', '#336699', '#0099AA', '#99CCCC', '#006655', 
                     '#669966', '#669922', '#99CC99', '#999966', '#CCCC99',
                     '#CCCC77', '#CCCC11', '#337799', '#336666', '#004499', 
                     '#333366', '#3399CC', '#99CC66', '#CC9900', '#CC6633',
                     '#663333', '#CC3333', '#993399', '#CC99CC', '#CCCCCC'),2))
  
  plot <- ggplot(processed_rm_out, aes(as.numeric(itg), fill = superfamily)) +
    #call geom_histogram with position="dodge" to offset the bars and manual binwidth of 2
    geom_histogram(position = "stack", binwidth = 0.02, color = "white", alpha = 0.8) +
    scale_fill_manual(values=mycolor) +
    ggtitle(paste0("Genome: ", genomename)) +
    xlab("Integrity") + ylab("TE loci") +
    scale_y_continuous(breaks = seq(0, ymax, by = ybreak)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, ymax)) +  # <-- changed from scale_x_continuous(limits=... and scale_y_continuous(limits=..
    theme(title = element_text(size = 8),
          axis.text.y = element_text(size = rel(1.8)),
          axis.text.x = element_text(size = rel(1.8)),
          axis.title.x = element_text(size = rel(1.8)),
          axis.title.y = element_text(size = rel(1.8)),
          legend.title = element_blank(),
          legend.text = element_text(size = rel(1)),
          legend.key.size = unit(0.4, "cm"),
          legend.background = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"))
  return(plot)
}


# ---------------------------------------------------------------------------
# Replicate-analysis helpers (item28/29: mean +/- SD across duplicated runs)
# ---------------------------------------------------------------------------

# Make sure a fasta has a .fai index (samtools), creating one if missing.
ensure_fai <- function(fasta_path){
  fai_path <- paste0(fasta_path, ".fai")
  if (!file.exists(fai_path)) {
    system2("samtools", c("faidx", fasta_path))
  }
  return(fai_path)
}

# Dodged TE vs non-TE bar plot with mean +/- SD error bars across replicates.
# df_summary needs columns: simulator, feature (TE/nonTE), mean_bases, sd_bases, mean_perc
# (sd_bases/sd_perc may be NA for a non-replicated genome, e.g. the original reference)
plot_te_nonte_bars_meansd <- function(df_summary, title_text){
  plot <- ggplot(df_summary, aes(x = simulator, y = mean_bases, fill = feature)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "white") +
    geom_errorbar(aes(ymin = mean_bases - sd_bases, ymax = mean_bases + sd_bases),
                  position = position_dodge(width = 0.8), width = 0.25, na.rm = TRUE) +
    geom_text(aes(y = mean_bases + ifelse(is.na(sd_bases), 0, sd_bases),
                  label = sprintf("%.1f%%", mean_perc)),
              position = position_dodge(width = 0.8), vjust = -0.5, size = rel(3)) +
    scale_fill_manual(values = c("TE" = "black", "nonTE" = "grey50")) +
    scale_y_continuous(labels = label_number(scale = 1e-6), expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Simulator", y = "Million bases", fill = "Feature") +
    ggtitle(title_text) +
    theme(title = element_text(size = 8),
              axis.text.y = element_text(size = rel(1.3)),
              axis.text.x = element_text(size = rel(1.3)),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size = rel(1.4)),
              legend.title = element_blank(),
              legend.text = element_text(size = rel(1)),
              legend.key.size = unit(0.5, "cm"),
              legend.background = element_blank(),
              panel.background = element_rect(fill = "white", color = "black"))
  return(plot)
}

# Bin a per-locus numeric column (e.g. integrity) into fixed-width bins per replicate,
# express each bin as a proportion of that replicate's total loci (so reps with slightly
# different total loci counts are normalised), then summarise the mean and SD of the
# proportion across replicates for each bin.
summarise_binned_proportion <- function(df, value_col, binwidth = 0.02, range = c(0, 1)){
  breaks <- seq(range[1], range[2], by = binwidth)
  df <- df |> dplyr::filter(.data[[value_col]] >= range[1], .data[[value_col]] <= range[2])
  # round before cutting: values like div = 1 - idn land as e.g. 0.09999999999999998
  # due to floating-point arithmetic, which is indistinguishable from the intended 0.1
  # at any precision this data actually uses, but sits fractionally on the wrong side
  # of a bin edge that coincides exactly with it - misclassifying every such locus into
  # the bin below and creating an artificial dip/spike pair at that boundary
  rounded_val <- round(df[[value_col]], 8)
  df <- df |> dplyr::mutate(.bin = cut(rounded_val, breaks = breaks, include.lowest = TRUE, right = FALSE))
  bin_mid <- breaks[-length(breaks)] + binwidth / 2
  names(bin_mid) <- levels(df$.bin)
  totals <- df |> dplyr::group_by(rep) |> dplyr::summarise(total = dplyr::n(), .groups = "drop")
  per_rep <- df |>
    dplyr::group_by(rep, .bin) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
    tidyr::complete(rep, .bin, fill = list(count = 0)) |>
    dplyr::left_join(totals, by = "rep") |>
    dplyr::mutate(prop = count / total)
  per_rep |>
    dplyr::group_by(.bin) |>
    dplyr::summarise(mean_prop = mean(prop), sd_prop = sd(prop),
                      mean_count = mean(count), sd_count = sd(count),
                      n_rep = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(bin_mid = bin_mid[as.character(.bin)]) |>
    dplyr::arrange(bin_mid)
}

# Extract per-bin mean +/- SD directly from an already-built, per-replicate faceted
# histogram (e.g. the output of plot_sim_integrity_by_loci_reps()/
# plot_sim_divergence_by_loci_reps()), via ggplot_build(), instead of recomputing bins
# independently with cut(). This guarantees the line/ribbon summary uses EXACTLY the
# same bins - position, width, and closed side - as the histogram it's paired with,
# since both come from the same stat_bin computation; no separate binning logic to keep
# in sync (and no risk of the two disagreeing the way cut()'s left-closed default and
# geom_histogram's right-closed default did).
summarise_bins_from_histogram <- function(plot){
  build_data <- ggplot2::ggplot_build(plot)$data[[1]]

  # sum across fill groups (superfamily) within each replicate panel and bin, to get
  # the stacked/combined total per (rep, bin) - matches the histogram's bar height
  per_rep <- build_data |>
    dplyr::group_by(PANEL, xmin, xmax) |>
    dplyr::summarise(count = sum(count), .groups = "drop")

  # per-replicate total loci (across all bins), for a proportion column alongside count
  totals <- per_rep |> dplyr::group_by(PANEL) |> dplyr::summarise(total = sum(count), .groups = "drop")

  per_rep <- per_rep |>
    dplyr::left_join(totals, by = "PANEL") |>
    dplyr::mutate(prop = count / total)

  per_rep |>
    dplyr::group_by(xmin, xmax) |>
    dplyr::summarise(mean_count = mean(count), sd_count = sd(count),
                      mean_prop = mean(prop), sd_prop = sd(prop),
                      n_rep = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(bin_mid = (xmin + xmax) / 2) |>
    dplyr::arrange(bin_mid)
}

# Line + ribbon plot of mean value per bin, +/- SD across replicates (all superfamilies
# pooled). Used for the main-figure E/F panels. Set use_count = TRUE to plot the raw
# number of TE loci per bin instead of the proportion of each replicate's total loci.
plot_binned_mean_sd <- function(bin_summary, genomename, xlab_text = "Integrity", ymax = NULL, use_count = FALSE, xlim = c(0, 1)){
  # label reflects the actual number of replicates behind this genome's bins - for a
  # single, non-replicated genome (e.g. Ori) n_rep is 1 and sd is undefined (NA), so
  # "mean +/- SD, n=5 reps" would be wrong here
  n_rep_val <- suppressWarnings(max(bin_summary$n_rep, na.rm = TRUE))
  replicated <- is.finite(n_rep_val) && n_rep_val > 1
  stat_label <- if (replicated) sprintf("mean ± SD, n=%d reps", n_rep_val) else "single genome, not replicated"

  if (use_count) {
    bin_summary <- bin_summary |> dplyr::rename(mean_val = mean_count, sd_val = sd_count)
    ylab_text <- paste0("Number of TE loci\n(", stat_label, ")")
  } else {
    bin_summary <- bin_summary |> dplyr::rename(mean_val = mean_prop, sd_val = sd_prop)
    ylab_text <- paste0("Proportion of TE loci\n(", stat_label, ")")
  }

  plot <- ggplot(bin_summary, aes(x = bin_mid, y = mean_val))
  # skip the ribbon layer entirely when there's no real SD to show (e.g. a single,
  # non-replicated genome) - avoids a "removed N rows" warning for a layer that would
  # never draw anything anyway
  if (replicated) {
    plot <- plot +
      geom_ribbon(aes(ymin = pmax(mean_val - sd_val, 0), ymax = mean_val + sd_val),
                  fill = "red", alpha = 0.4)
  }
  plot <- plot +
    geom_line(color = "black", linewidth = 0.4) +
    geom_point(color = "black", size = 1) +
    ggtitle(paste0("Genome: ", genomename)) +
    xlab(xlab_text) + ylab(ylab_text) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    coord_cartesian(xlim = xlim, ylim = c(0, ymax)) +
    theme(title = element_text(size = 8),
          axis.text.y = element_text(size = rel(1.8)),
          axis.text.x = element_text(size = rel(1.8)),
          axis.title.x = element_text(size = rel(1.8)),
          axis.title.y = element_text(size = rel(1.3)),
          panel.background = element_rect(fill = "white", color = "black"))
  return(plot)
}

# Line + ribbon plot of mean value per bin, +/- SD across replicates, for MULTIPLE
# groups (e.g. simulation scenarios) overlaid on one plot and distinguished by colour -
# the multi-group analogue of plot_binned_mean_sd(). bin_summary must have a column
# named group_col identifying which group each row belongs to: build it by calling
# summarise_bins_from_histogram() once per group, tagging each result with
# mutate(<group_col> = "...") before dplyr::bind_rows()-ing them together.
plot_binned_mean_sd_bygroup <- function(bin_summary, group_col, group_colors, xlab_text = "Divergence (p-dist)", ymax = NULL, use_count = TRUE, xlim = c(0, 1), alpha_ribbon = 0.25){
  if (use_count) {
    bin_summary <- bin_summary |> dplyr::rename(mean_val = mean_count, sd_val = sd_count)
    ylab_text <- "Number of TE loci\n(mean ± SD across replicates)"
  } else {
    bin_summary <- bin_summary |> dplyr::rename(mean_val = mean_prop, sd_val = sd_prop)
    ylab_text <- "Proportion of TE loci\n(mean ± SD across replicates)"
  }
  bin_summary <- bin_summary |> dplyr::mutate(.group = .data[[group_col]])

  plot <- ggplot(bin_summary, aes(x = bin_mid, y = mean_val, color = .group, fill = .group)) +
    geom_ribbon(aes(ymin = pmax(mean_val - sd_val, 0), ymax = mean_val + sd_val), color = NA, alpha = alpha_ribbon, na.rm = TRUE) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.25) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    xlab(xlab_text) + ylab(ylab_text) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
    coord_cartesian(xlim = xlim, ylim = c(0, ymax)) +
    theme(title = element_text(size = 10),
          axis.text.y = element_text(size = rel(1.3)),
          axis.text.x = element_text(size = rel(1.3)),
          axis.title.x = element_text(size = rel(1.4)),
          axis.title.y = element_text(size = rel(1.2)),
          legend.title = element_blank(),
          legend.text = element_text(size = rel(0.9)),
          legend.key.size = unit(0.4, "cm"),
          legend.background = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"))
  return(plot)
}

# Same visual style as plot_sim_integrity_by_loci(), but faceted by replicate
# (Supplementary use: keeps the original per-locus, stacked-by-superfamily view,
# one small-multiple panel per replicate).
plot_sim_integrity_by_loci_reps <- function(gff_breakdown, ymax, ybreak, genomename){
  # same fixed, locale-independent order/colour as the pie charts and bar charts
  te_color_map <- te_superfamily_colors(gff_breakdown$sfam)
  gff_breakdown <- gff_breakdown |> dplyr::mutate(sfam = factor(as.character(sfam), levels = order_te_superfamilies(sfam)))

  plot <- ggplot(gff_breakdown, aes(round(as.numeric(itg), 8), fill = sfam)) +
    # boundary = 0 matches summarise_binned_proportion()'s bin edges (0, 0.02, 0.04, ...)
    # exactly, instead of geom_histogram's own default (edges offset by half a binwidth);
    # round() avoids floating-point noise (e.g. 0.09999999999999998) misclassifying a
    # locus into the bin below when the true value sits exactly on a bin edge
    geom_histogram(position = "stack", binwidth = 0.02, boundary = 0, closed = "left", color = "white", alpha = 0.8) +
    scale_fill_manual(values = te_color_map) +
    facet_wrap(~rep, ncol = 5, labeller = labeller(rep = function(x) paste0("rep", x))) +
    ggtitle(paste0("Genome: ", genomename, " (per replicate)")) +
    xlab("Integrity") + ylab("TE loci") +
    scale_y_continuous(breaks = seq(0, ymax, by = ybreak)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, ymax)) +
    theme(title = element_text(size = 8),
          axis.text.y = element_text(size = rel(1.3)),
          axis.text.x = element_text(size = rel(1.3)),
          axis.title.x = element_text(size = rel(1.4)),
          axis.title.y = element_text(size = rel(1.4)),
          legend.title = element_blank(),
          legend.text = element_text(size = rel(0.8)),
          legend.key.size = unit(0.3, "cm"),
          legend.background = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"))
  return(plot)
}

# Same visual style as plot_sim_divergence_by_loci(), but faceted by replicate
# (Supplementary use: keeps the original per-locus, stacked-by-superfamily view,
# one small-multiple panel per replicate).
plot_sim_divergence_by_loci_reps <- function(gff_breakdown, ymax, ybreak, genomename){
  # same fixed, locale-independent order/colour as the pie charts and bar charts
  te_color_map <- te_superfamily_colors(gff_breakdown$sfam)
  gff_breakdown <- gff_breakdown |> dplyr::mutate(sfam = factor(as.character(sfam), levels = order_te_superfamilies(sfam)))

  plot <- ggplot(gff_breakdown, aes(round(as.numeric(div), 8), fill = sfam)) +
    # see plot_sim_integrity_by_loci_reps() - boundary = 0 + round() keep this on the
    # same bin grid as summarise_binned_proportion() and avoid floating-point noise
    geom_histogram(position = "stack", binwidth = 0.02, boundary = 0, closed = "left", color = "white", alpha = 0.8) +
    scale_fill_manual(values = te_color_map) +
    facet_wrap(~rep, ncol = 5, labeller = labeller(rep = function(x) paste0("rep", x))) +
    ggtitle(paste0("Genome: ", genomename, " (per replicate)")) +
    xlab("Divergence (p-dist)") + ylab("TE loci") +
    scale_y_continuous(breaks = seq(0, ymax, by = ybreak)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    coord_cartesian(xlim = c(0, 0.6), ylim = c(0, ymax)) +
    theme(title = element_text(size = 8),
          axis.text.y = element_text(size = rel(1.3)),
          axis.text.x = element_text(size = rel(1.3)),
          axis.title.x = element_text(size = rel(1.4)),
          axis.title.y = element_text(size = rel(1.4)),
          legend.title = element_blank(),
          legend.text = element_text(size = rel(0.8)),
          legend.key.size = unit(0.3, "cm"),
          legend.background = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"))
  return(plot)
}

# Per-superfamily TE composition (loci count, bp, family count), mean +/- SD across a
# list of per-replicate GFFs for one simulator. A superfamily missing from a given
# replicate is filled with zero before averaging (not dropped), since absence there
# means "zero copies that run", not missing data.
te_superfamily_meansd <- function(tegff_list, simulator_name){
  per_rep <- dplyr::bind_rows(lapply(seq_along(tegff_list), function(i){
    s <- te_by_superfamily_sim_genome(tegff_list[[i]])
    s$rep <- i
    s
  }))
  per_rep <- per_rep |>
    tidyr::complete(rep, TE, fill = list(count = 0, bp = 0, family_count = 0))
  per_rep |>
    dplyr::group_by(TE) |>
    dplyr::summarise(mean_count = mean(count), sd_count = sd(count),
                      mean_bp = mean(bp), sd_bp = sd(bp),
                      mean_family = mean(family_count), sd_family = sd(family_count),
                      .groups = "drop") |>
    dplyr::mutate(simulator = simulator_name)
}

# Named TE-superfamily -> colour map, keyed by alphabetical order of the superfamily
# name. plot_te_loci_piechart()/plot_te_bp_piechart() take their category order from
# te_by_superfamily_sim_genome() (a group_by(superfamily) |> summarise(), which dplyr
# returns sorted alphabetically) and assign the same 30-colour palette positionally in
# that order, so reproducing "alphabetical rank -> palette position" here keeps a given
# superfamily on the same colour in both the pie charts and any bar chart using this map.
# Canonical TE superfamily display order. This is the order dplyr's
# group_by(superfamily) |> summarise() inside te_by_superfamily_sim_genome() produces
# (dplyr groups characters using a locale-independent, C-locale/byte-order sort) - it is
# NOT necessarily what base R's sort() gives, since that respects the session's collation
# locale (e.g. under en_US.UTF-8, sort() puts "DNA/hAT" before "DNA/Mariner"/"DNA/MULE",
# which disagrees with dplyr and with the pie charts). Hardcoding it here keeps pie
# charts, bar charts and their colours in sync regardless of the machine's locale.
TE_SUPERFAMILY_ORDER <- c(
  "DNA/CACTA", "DNA/Harbinger", "DNA/hAT", "DNA/Mariner", "DNA/MULE",
  "LINE/L1", "LINE/unknown", "LTR/Copia", "LTR/Solo", "LTR/Ty3", "LTR/unknown",
  "MITE/Stow", "MITE/Tourist", "RC/Helitron", "SINE/tRNA", "SINE/unknown"
)

# Order a set of TE superfamily names: known categories first (in TE_SUPERFAMILY_ORDER),
# any category not in that fixed list appended afterwards via a locale-independent
# (radix/C-locale) sort, so an unexpected new category never errors out or gets dropped.
order_te_superfamilies <- function(te_names){
  present <- unique(as.character(te_names))
  known <- TE_SUPERFAMILY_ORDER[TE_SUPERFAMILY_ORDER %in% present]
  extra <- unique(setdiff(present, TE_SUPERFAMILY_ORDER))
  extra <- extra[order(tolower(extra), method = "radix")]
  c(known, extra)
}

te_superfamily_colors <- function(te_names){
  mycolor <- c(rep(c('#CC9999', '#CC6666', '#993344', '#996699', '#996633',
                     '#666699', '#336699', '#0099AA', '#99CCCC', '#006655',
                     '#669966', '#669922', '#99CC99', '#999966', '#CCCC99',
                     '#CCCC77', '#CCCC11', '#337799', '#336666', '#004499',
                     '#333366', '#3399CC', '#99CC66', '#CC9900', '#CC6633',
                     '#663333', '#CC3333', '#993399', '#CC99CC', '#CCCCCC'),2))
  te_sorted <- order_te_superfamilies(te_names)
  stats::setNames(mycolor[seq_along(te_sorted)], te_sorted)
}

# Bar plot of mean +/- SD across replicates for one value column (loci count, bp, or
# family count), one panel per simulator (facet_wrap(~simulator)) and coloured by TE
# superfamily - the mean+-SD analogue of plot_te_loci_piechart()/plot_te_bp_piechart().
# Reuses their TE2 trick (a zero-padded rank prefix on the label, e.g. "01_DNA/MULE")
# to pin down category order explicitly rather than relying on default factor/alpha
# ordering; here the rank comes from sort_by instead of row order, and the colour for
# each TE2 label is still looked up by the underlying superfamily name (alphabetically,
# via te_superfamily_colors()) so bars stay colour-matched to the pie charts regardless
# of what order sort_by puts them in.
plot_te_superfamily_bars_meansd <- function(df_summary, value_col, sd_col, title_text, ylab_text, scale_million = FALSE, sort_by = NULL){
  te_color_by_name <- te_superfamily_colors(df_summary$TE)

  if (!is.null(sort_by)) {
    te_order <- df_summary |>
      dplyr::group_by(TE) |>
      dplyr::summarise(.sort_val = sum(.data[[sort_by]], na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(.sort_val)) |>
      dplyr::pull(TE)
  } else {
    te_order <- order_te_superfamilies(df_summary$TE)
  }
  te_order <- as.character(te_order)

  # order via factor levels only - no digit prefix baked into the visible label (unlike
  # the pie charts' TE2), since a bar chart already shows the order positionally on the
  # x-axis/facet and a "01_" prefix on every tick label would just be visual clutter
  df_summary <- df_summary |>
    dplyr::mutate(mean_val = .data[[value_col]], sd_val = .data[[sd_col]]) |>
    dplyr::mutate(TE = factor(as.character(TE), levels = te_order))

  plot <- ggplot(df_summary, aes(x = TE, y = mean_val, fill = TE)) +
    geom_col(width = 0.7, color = "white") +
    geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), width = 0.3, na.rm = TRUE) +
    scale_fill_manual(values = te_color_by_name) +
    facet_wrap(~simulator, ncol = 2) +
    labs(x = "TE superfamily", y = ylab_text, fill = "TE superfamily") +
    ggtitle(title_text) +
    theme(title = element_text(size = 8),
          axis.text.y = element_text(size = rel(1.3)),
          axis.text.x = element_text(size = rel(1.2), angle = 60, hjust = 1, vjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = rel(1.4)),
          legend.position = "none",
          strip.text = element_text(size = rel(1.1)),
          panel.background = element_rect(fill = "white", color = "black"))
  if (scale_million) plot <- plot + scale_y_continuous(labels = label_number(scale = 1e-6))
  return(plot)
}



# Compute Pearson and Spearman correlation between TE loci divergence and integrity
# for one genome/dataset. `df` must contain numeric columns `div` and `itg`.
# `genome` and `dataset` are labels carried through to identify the row in the
# combined summary table (e.g. genome = "tair10", dataset = "Real (original)").
compute_corr_identity_integrity <- function(df, genome, dataset){
  df <- df |> dplyr::filter(is.finite(div), is.finite(itg))
  ct_p <- cor.test(df$div, df$itg, method = "pearson")
  ct_s <- cor.test(df$div, df$itg, method = "spearman", exact = FALSE)

  data.frame(
    genome = genome,
    dataset = dataset,
    n_loci = nrow(df),
    pearson_r = unname(ct_p$estimate),
    pearson_ci_low = ct_p$conf.int[1],
    pearson_ci_high = ct_p$conf.int[2],
    pearson_p = ct_p$p.value,
    spearman_rho = unname(ct_s$estimate),
    spearman_p = ct_s$p.value
  )
}

# Family-level robustness check for compute_corr_identity_integrity(): TE loci
# within a family are not independent observations (shared ancestry), which
# inflates the per-locus test's apparent significance - the same
# non-independence issue raised for the Shannon-entropy/Wilcoxon analysis
# elsewhere in this revision. This collapses each family to its median
# divergence and median integrity first (one data point per family, however
# many loci it has), then runs the same correlation test on those medians.
# `df` must contain `div`, `itg` and `family` columns.
compute_corr_identity_integrity_family <- function(df, genome, dataset){
  if (!"family" %in% names(df)) {
    stop("compute_corr_identity_integrity_family(): df has no 'family' column for genome '",
         genome, "', dataset '", dataset, "'. If this is a mode-2 breakdown data frame, ",
         "rebuild it with te_breakdown_sim_genome() after sourcing the current functions.R - ",
         "a copy created before the 'family' column was added won't have it.")
  }
  df_fam <- df |>
    dplyr::filter(is.finite(div), is.finite(itg), !is.na(family), family != "") |>
    dplyr::group_by(family) |>
    dplyr::summarise(div = median(div), itg = median(itg), .groups = "drop")
  out <- compute_corr_identity_integrity(df_fam, genome, dataset)
  names(out)[names(out) == "n_loci"] <- "n_families"
  out
}

# Divergence vs. integrity for ONE genome, real vs. mode-2 simulated side by
# side. The real and simulated loci counts for the same species are always
# close in magnitude (they are a matched pair), so both panels can safely
# share one raw-count log10 fill scale - the legend shows genuine, comparable
# TE loci counts.
plot_corr_identity_integrity_one_genome <- function(df_genome, genome){
  plot <- ggplot(df_genome, aes(x = div, y = itg)) +
    # binwidth matches the 0.02 bins used by the divergence/integrity histograms
    # (plot_ori_divergence_by_loci() etc.) elsewhere in this script
    geom_bin2d(binwidth = c(0.02, 0.02)) +
    scale_fill_viridis_c(trans = "log10", name = "TE loci") +
    geom_smooth(method = "lm", formula = y ~ x, color = "#CC3333", linewidth = 0.6, se = TRUE) +
    facet_wrap(~dataset, nrow = 1) +
    labs(x = "Divergence (1 - identity)", y = "Integrity", title = paste0("Genome: ", genome)) +
    coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 1)) +
    theme(title = element_text(size = 9),
          axis.text.y = element_text(size = rel(1)),
          axis.text.x = element_text(size = rel(1)),
          axis.title.x = element_text(size = rel(1.1)),
          axis.title.y = element_text(size = rel(1.1)),
          strip.text = element_text(size = rel(1)),
          panel.background = element_rect(fill = "white", color = "black"))
  return(plot)
}

# Divergence vs. integrity for ONE mode-1 identity-SD scenario (5 reps
# pooled), a single panel - used to visualise the "by construction ~0"
# baseline scenario by scenario. Same binwidth as the histograms and the
# per-genome real-vs-mode-2 plot above, for visual consistency.
plot_corr_identity_integrity_scenario <- function(df_scenario, scenario_label){
  plot <- ggplot(df_scenario, aes(x = div, y = itg)) +
    geom_bin2d(binwidth = c(0.02, 0.02)) +
    scale_fill_viridis_c(trans = "log10", name = "TE loci") +
    geom_smooth(method = "lm", formula = y ~ x, color = "#CC3333", linewidth = 0.6, se = TRUE) +
    labs(x = "Divergence (1 - identity)", y = "Integrity", title = scenario_label) +
    coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 1)) +
    theme(title = element_text(size = 9),
          axis.text.y = element_text(size = rel(1)),
          axis.text.x = element_text(size = rel(1)),
          axis.title.x = element_text(size = rel(1.1)),
          axis.title.y = element_text(size = rel(1.1)),
          panel.background = element_rect(fill = "white", color = "black"))
  return(plot)
}

