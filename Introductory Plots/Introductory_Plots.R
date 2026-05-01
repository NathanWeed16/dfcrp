library(ggplot2)
library(ggforce)

# Code to visualize the full lunar highlands data
pdf("Full_Image_Craters.pdf", width = 8, height = 4)
ggplot(crater_data) + geom_circle(aes(x0 = crater_data[,1], y0 = crater_data[,2], 
                                      r = crater_data[,3]/2), color = "gray60", 
                                  fill = NA, size = 0.5) + 
  coord_fixed() + 
  theme_void()
dev.off()


# Code to visualize crater identification and specification variation between experts
interest <- which(crater_data[,1] > 2850 & crater_data[,1] < 3100 & crater_data[,2] > -1700 & crater_data[,2] < -1500)
subset <- crater_data[interest, ]

x_pad <- diff(range(subset[,1])) * 0.2
y_pad <- diff(range(subset[,2])) * 0.5
x_center <- mean(range(subset[,1]))

for (i in 1:11){

  pdf(paste0(LETTERS[i], ".pdf"), width = 4, height = 4)
  p <- ggplot(subset) +
    geom_circle(
      aes(
        x0 = subset[,1],
        y0 = subset[,2],
        r = subset[,3] / 2,
        color = (subset[,4] == i)
      ),
      fill = NA,
      size = 1.5
    ) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray80")) +
    coord_fixed(
      xlim = range(subset[,1]) + c(-x_pad, x_pad),
      ylim = range(subset[,2]) + c(-y_pad, y_pad)
    ) +
    theme_void() +
    theme(legend.position = "none") +
    annotate(
      "text",
      x = x_center,
      y = Inf,
      label = LETTERS[i],
      fontface = "bold",
      vjust = 1.5,
      size = 6
    )
    print(p)
    dev.off()
    
}
