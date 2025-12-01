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
    dplyr::mutate(div = 1 - as.numeric(idn))

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

  data.sum <- data.sum |> mutate(TE2 = sprintf("%02d_%s", row_number(), TE)) |> 
    dplyr::arrange(TE2) |> dplyr::mutate(cumulative = cumsum(count), 
                                        percentage = round(100*(count/total), 1), 
                                        midpoint = (cumulative - count/2), 
                                        label = paste0(percentage, "%")
                                        )

  # set up color pallet
  mycolor <- c(rep(c('#CC9999', '#CC6666', '#993344', '#996699', '#996633',
               '#666699', '#336699', '#0099AA', '#99CCCC', '#006655', 
               '#669966', '#669922', '#99CC99', '#999966', '#CCCC99',
               '#CCCC77', '#CCCC11', '#337799', '#336666', '#004499', 
               '#333366', '#3399CC', '#99CC66', '#CC9900', '#CC6633',
               '#663333', '#CC3333', '#993399', '#CC99CC', '#CCCCCC'),2))

  # generate plot
  plot <- ggplot(data.sum, aes(x = "", y = count, fill = TE2)) +
    geom_bar(stat = "identity", color = 'white', position = position_stack(reverse = TRUE), width = 1, linewidth = 0.5, alpha = 0.8) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values=mycolor) +
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
    data.sum <- data.sum |> mutate(TE2 = sprintf("%02d_%s", row_number(), TE)) |> 
      dplyr::arrange(TE2) |> dplyr::mutate(cumulative = cumsum(bp), percentage = round(100*(bp/total), 1), midpoint = (cumulative - bp/2), label = paste0(percentage, "%"))

    # set up color pallet
    mycolor <- c(rep(c('#CC9999', '#CC6666', '#993344', '#996699', '#996633',
                       '#666699', '#336699', '#0099AA', '#99CCCC', '#006655', 
                       '#669966', '#669922', '#99CC99', '#999966', '#CCCC99',
                       '#CCCC77', '#CCCC11', '#337799', '#336666', '#004499', 
                       '#333366', '#3399CC', '#99CC66', '#CC9900', '#CC6633',
                       '#663333', '#CC3333', '#993399', '#CC99CC', '#CCCCCC'),2))
    # generate plot
    plot <- ggplot(data.sum[,-1], aes(x = "", y = bp, fill = TE2)) +
    geom_bar(stat = "identity", color = 'white', position = position_stack(reverse = TRUE), width = 1, linewidth = 0.5, alpha = 0.8) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values=mycolor) +
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
  # set up color pallet
  mycolor <- c(rep(c('#CC9999', '#CC6666', '#993344', '#996699', '#996633',
                     '#666699', '#336699', '#0099AA', '#99CCCC', '#006655', 
                     '#669966', '#669922', '#99CC99', '#999966', '#CCCC99',
                     '#CCCC77', '#CCCC11', '#337799', '#336666', '#004499', 
                     '#333366', '#3399CC', '#99CC66', '#CC9900', '#CC6633',
                     '#663333', '#CC3333', '#993399', '#CC99CC', '#CCCCCC'),2))
  
  plot <- ggplot(gff_breakdown, aes(as.numeric(div), fill = sfam)) +
    #call geom_histogram with position="dodge" to offset the bars and manual binwidth of 2
    geom_histogram(position = "stack", binwidth = 0.02, color = "white", alpha = 0.8) +
    scale_fill_manual(values=mycolor) +
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
  # set up color pallet
  mycolor <- c(rep(c('#CC9999', '#CC6666', '#993344', '#996699', '#996633',
                     '#666699', '#336699', '#0099AA', '#99CCCC', '#006655', 
                     '#669966', '#669922', '#99CC99', '#999966', '#CCCC99',
                     '#CCCC77', '#CCCC11', '#337799', '#336666', '#004499', 
                     '#333366', '#3399CC', '#99CC66', '#CC9900', '#CC6633',
                     '#663333', '#CC3333', '#993399', '#CC99CC', '#CCCCCC'),2))

  plot <- ggplot(gff_breakdown, aes(as.numeric(itg), fill = sfam)) +
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