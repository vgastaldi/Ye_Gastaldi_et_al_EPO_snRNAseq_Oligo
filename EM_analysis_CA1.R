#### Data analysis - Xuan's EM data ####
## Loading required info
library(ggplot2)
library(car)
library(openxlsx)
library(dplyr)
#devtools::install_github("coolbutuseless/ggpattern")
library(ggpattern)
library(coin)
set.seed(0)

#### Comparison of myelinated axons per area - CA1 ####
as.data.frame(read.xlsx("CNPcreEPORfl_CA1_myelinated_average_formatted.xlsx")) -> myelinated_area

## Normality
shapiro.test(myelinated_area[which(myelinated_area$group == "Control"),3]) # 0.3871
shapiro.test(myelinated_area[which(myelinated_area$group == "KO"),3]) # 0.0147
# Not-normally distributed

## Variance
leveneTest(myelinated_area$average_myelinated ~ myelinated_area$group) # 0.5104
# Homogeneous variance

# Statistical test
wilcox.test(myelinated_area$average_myelinated ~ myelinated_area$group)
# 0.3725

#### Figure 1 ####
df <- myelinated_area

df_ko <- df[df$group == "KO", ]
df_control <- df[df$group == "Control", ]

summary_df_control <- df_control %>%
  summarise(sem = sd(average_myelinated, na.rm = TRUE) / sqrt(n()),
            average_myelinated = mean(average_myelinated, na.rm = TRUE))
summary_df_control$group <- "Control"

summary_df_ko <- df_ko %>%
  summarise(sem = sd(average_myelinated, na.rm = TRUE) / sqrt(n()),
            average_myelinated = mean(average_myelinated, na.rm = TRUE))
summary_df_ko$group <- "KO"

# Calculate the y position of the comparison bars
y_pos <- max(summary_df_ko$average_myelinated,summary_df_control$average_myelinated)+max(summary_df_ko$sem,summary_df_control$sem) * 2.3
y_pos2 <- max(summary_df_ko$average_myelinated,summary_df_control$average_myelinated)+max(summary_df_ko$sem,summary_df_control$sem) * 2

# Set the spacing between the comparison bars
spacing <- max(summary_df_ko$average_myelinated,summary_df_control$average_myelinated)+max(summary_df_ko$sem,summary_df_control$sem) * 0.04

# Create the bar plot
ggplot(df, aes(x = group, y = average_myelinated)) +
  geom_bar(data = summary_df_control, aes(fill = group), stat = "identity") +
  geom_bar_pattern(data = summary_df_ko, aes(fill = group, pattern = group), 
                   stat = "identity", 
                   pattern = "stripe",
                   pattern_fill = "blue",
                   pattern_angle = 315,
                   pattern_density = 0.1,
                   pattern_spacing = 0.02,
                   pattern_key_scale_factor = 1) +
  geom_point(position = position_jitter(width = 0.1)) +
  geom_errorbar(data = summary_df_control, aes(ymin = average_myelinated - sem, ymax = average_myelinated + sem), width = 0.2, position = position_dodge(0.9)) +
  geom_errorbar(data = summary_df_ko, aes(ymin = average_myelinated - sem, ymax = average_myelinated + sem), width = 0.2, position = position_dodge(0.9)) +
  theme_minimal()+ 
scale_fill_manual(values = c("Control" = "#D4D4D4", "KO" = "#D4D4D4")) +
  annotate("text", x = 1.5, y = y_pos, label = "0.372", parse = TRUE, size = 5) +
  annotate("segment", x = 1, xend = 2, y = y_pos2, yend = y_pos2) +
  theme(axis.title.x=element_blank(), # Suppress the X axis label
        axis.title.y=element_text(size=14,face="bold"), # Increase size of Y axis label and make it bold
        axis.text=element_text(size=12), # Increase size of axis text
        axis.text.x = element_text(face = "bold"), # Make the group names bold
        legend.position="none", # Suppress the legend
        panel.grid.major.x = element_blank()) + # Remove background vertical lines
  ylab("average of myelinated axons/µm^2") # Set the Y axis label to the value of measurement_plotted
ggsave("myelinated_axons_average_CA1.png", last_plot(), width = 10, height = 10, dpi = 600)