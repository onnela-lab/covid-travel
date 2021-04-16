#' This function plots our data in both histogram and boxplot.
#' @param values numeric vector
#' @import cowplot
#' @examples
#' \dontrun{#' x = 1:9
#' hisboxfunc(x) }
#' @export

hisboxfunc = function(values){
data = as.data.frame(values)
colnames(data) = "values"
plt1 = ggplot(data, aes(x=values))+geom_histogram(color="black", fill="blue")+ ggtitle("Histogram")


plt2 = ggplot(data, aes(y=values))+
  geom_boxplot(fill = "green")+  ggtitle("Boxplot")

plt = ggdraw() +draw_plot(plt1, x = 0, y = 0, width = .5, height = 1) +
  draw_plot(plt2, x = .5, y = 0, width = .5, height = 1)
return(plt)
}


