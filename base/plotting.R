# plot histograms of gene and pathway distributions
generate_histogram_from_list <- function(input_list){
  df <- data.frame(values=input_list)
  p <-  ggplot(df, aes(x=values)) +
    geom_histogram(color='black', fill='#ff5252') + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(), text=element_text(family='Trebuchet MS', face='bold'), axis.title=element_text(size=12), axis.text=element_text(size=13), plot.title=element_text(size = 15, hjust=0.5)) + ylab('Frequency')
  
  return(p)
}