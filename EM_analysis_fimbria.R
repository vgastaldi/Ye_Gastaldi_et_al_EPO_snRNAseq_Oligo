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

#### Comparison of myelinated axons per area ####
as.data.frame(read.xlsx("CNPcreEPORfl_myelinated_formatted.xlsx")) -> myelinated_area

## Normality
shapiro.test(myelinated_area[which(myelinated_area$group == "Control"),3]) # 0.2786
shapiro.test(myelinated_area[which(myelinated_area$group == "KO"),3]) # 0.2166
# Normally distributed

## Variance
leveneTest(myelinated_area$average_myelinated ~ myelinated_area$group) # 0.116
# Homogeneous variance

# Statistical test
t.test(myelinated_area$average_myelinated ~ myelinated_area$group,var.equal = T)
# 0.8865

#### G-Ratio analysis ####
as.data.frame(read.xlsx("CNPcreEPORfl_fimbria_gratio_format.xlsx")) -> gratio

gratio %>%
  group_by(animal,group) %>%
  summarise(average_in_out_diameter = mean(`inside/outside.diameter`, na.rm = TRUE)) -> gratio_average
as.data.frame(gratio_average) -> gratio_average

## Normality
shapiro.test(gratio_average[which(gratio_average$group == "Control"),3]) # 0.5019
shapiro.test(gratio_average[which(gratio_average$group == "KO"),3]) # 0.3264
# Normally distributed

## Variance
leveneTest(gratio_average$average_in_out_diameter ~ gratio_average$group) # 0.455
# Homogeneous variance

# Statistical test
t.test(gratio_average$average_in_out_diameter ~ gratio_average$group,var.equal = T)
# 0.3573


#### Testing for difference in the size of the axons ####
gratio %>%
  group_by(animal,group) %>%
  summarise(average_Inside.diameter = mean(Inside.diameter, na.rm = TRUE)) -> gratio_average_inside
as.data.frame(gratio_average_inside) -> gratio_average_inside

## Normality
shapiro.test(gratio_average_inside[which(gratio_average_inside$group == "Control"),3]) # 0.9457
shapiro.test(gratio_average_inside[which(gratio_average_inside$group == "KO"),3]) # 0.4441
# Normally distributed

## Variance
leveneTest(gratio_average_inside$average_Inside.diameter ~ gratio_average_inside$group) # 0.2951
# Homogeneous variance

# Statistical test
t.test(gratio_average_inside$average_Inside.diameter ~ gratio_average_inside$group,var.equal = T)
# 0.2765

# Using the two-sided Kolmogorov-Smirnov test for testing if they come from the same continuous distribution
ks.test(gratio_average_inside[which(gratio_average_inside$group == "Control"),"average_Inside.diameter"],gratio_average_inside[which(gratio_average_inside$group == "KO"),"average_Inside.diameter"])

# Using the whole table...
## Normality
shapiro.test(gratio[which(gratio$group == "Control"),7]) # not normal
shapiro.test(gratio[which(gratio$group == "KO"),7]) # not normal
# Normally distributed

## Variance
leveneTest(gratio[,7] ~ gratio$group) # 0.2366
# Homogeneous variance

# Statistical test
wilcox.test(gratio[,7] ~ gratio$group)
# 0.2887

ks.test(gratio[which(gratio$group == "Control"),"Inside.diameter"],gratio[which(gratio$group == "KO"),"Inside.diameter"],exact = F)

# Dealing with ties
control_diameters <- gratio[gratio$group == "Control", "Inside.diameter"]
ko_diameters <- gratio[gratio$group == "KO", "Inside.diameter"]

# Combine the diameters into a single vector and create a group factor
all_diameters <- c(control_diameters, ko_diameters)
group_factor <- factor(rep(c("Control", "KO"), c(length(control_diameters), length(ko_diameters))))

# Perform the one-way test
oneway_test(all_diameters ~ group_factor,ties.method = "mid-ranks") #  p-value = 3.21e-08

# Fit an ANOVA model
model <- aov(all_diameters ~ group_factor)

# Get the summary of the model
model_summary <- summary(model)

# Extract the sum of squares for the effect (SS_effect)
SS_effect <- model_summary[[1]]["group_factor", "Sum Sq"]

# Extract the total sum of squares (SS_total)
SS_total <- sum(model_summary[[1]][, "Sum Sq"])

# Calculate eta-squared (η²)
eta_squared <- SS_effect / SS_total

# Print the eta-squared value
print(eta_squared)

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
y_pos <- max(summary_df_ko$average_myelinated,summary_df_control$average_myelinated)+max(summary_df_ko$sem,summary_df_control$sem) * 2.2
y_pos2 <- max(summary_df_ko$average_myelinated,summary_df_control$average_myelinated)+max(summary_df_ko$sem,summary_df_control$sem) * 1.95

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
  annotate("text", x = 1.5, y = y_pos, label = "0.886", parse = TRUE, size = 5) +
  annotate("segment", x = 1, xend = 2, y = y_pos2, yend = y_pos2) +
  theme(axis.title.x=element_blank(), # Suppress the X axis label
        axis.title.y=element_text(size=14,face="bold"), # Increase size of Y axis label and make it bold
        axis.text=element_text(size=12), # Increase size of axis text
        axis.text.x = element_text(face = "bold"), # Make the group names bold
        legend.position="none", # Suppress the legend
        panel.grid.major.x = element_blank()) + # Remove background vertical lines
  ylab("average of myelinated axons/µm^2") # Set the Y axis label to the value of measurement_plotted
ggsave("myelinated_axons_average.png", last_plot(), width = 10, height = 10, dpi = 600)

#### Figure 2 ####
df <- gratio_average

df_ko <- df[df$group == "KO", ]
df_control <- df[df$group == "Control", ]

summary_df_control <- df_control %>%
  summarise(sem = sd(average_in_out_diameter, na.rm = TRUE) / sqrt(n()),
            average_in_out_diameter = mean(average_in_out_diameter, na.rm = TRUE))
summary_df_control$group <- "Control"

summary_df_ko <- df_ko %>%
  summarise(sem = sd(average_in_out_diameter, na.rm = TRUE) / sqrt(n()),
            average_in_out_diameter = mean(average_in_out_diameter, na.rm = TRUE))
summary_df_ko$group <- "KO"

# Calculate the y position of the comparison bars
y_pos <- max(summary_df_ko$average_in_out_diameter,summary_df_control$average_in_out_diameter)+max(summary_df_ko$sem,summary_df_control$sem) * 6
y_pos2 <- max(summary_df_ko$average_in_out_diameter,summary_df_control$average_in_out_diameter)+max(summary_df_ko$sem,summary_df_control$sem) * 4.5

# Create the bar plot
ggplot(df, aes(x = group, y = average_in_out_diameter)) +
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
  geom_errorbar(data = summary_df_control, aes(ymin = average_in_out_diameter - sem, ymax = average_in_out_diameter + sem), width = 0.2, position = position_dodge(0.9)) +
  geom_errorbar(data = summary_df_ko, aes(ymin = average_in_out_diameter - sem, ymax = average_in_out_diameter + sem), width = 0.2, position = position_dodge(0.9)) +
  theme_minimal()+ 
  scale_fill_manual(values = c("Control" = "#D4D4D4", "KO" = "#D4D4D4")) +
  annotate("text", x = 1.5, y = y_pos, label = "0.357", parse = TRUE, size = 5) +
  annotate("segment", x = 1, xend = 2, y = y_pos2, yend = y_pos2) +
  theme(axis.title.x=element_blank(), # Suppress the X axis label
        axis.title.y=element_text(size=14,face="bold"), # Increase size of Y axis label and make it bold
        axis.text=element_text(size=12), # Increase size of axis text
        axis.text.x = element_text(face = "bold"), # Make the group names bold
        legend.position="none", # Suppress the legend
        panel.grid.major.x = element_blank()) + # Remove background vertical lines
  ylab("g-ratio") # Set the Y axis label to the value of measurement_plotted
ggsave("gratio.png", last_plot(), width = 10, height = 10, dpi = 600)

#### Figure 3 ####
ggplot(gratio, aes(x = Inside.diameter, y = `inside/outside.diameter`, color = group)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(face = "bold", size = 16), # Increase size of X axis text
    axis.text.y = element_text(face = "bold", size = 16), # Increase size of Y axis text
    legend.text = element_text(size = 16), # Increase size of legend text
    legend.title = element_text(size = 16, face = "bold") # Increase size and bold the legend title
  ) +
  ylab("g-ratio") + xlab("Diameter [µm]")
ggsave("gratio vs diameter.png", last_plot(), width = 10, height = 10, dpi = 600)

write.table(gratio,"fimbria_Figure5j_scatter.txt",sep = "\t",col.names = T,row.names = F,quote = F)

## Version 2 ##
ggplot(gratio, aes(x = Inside.diameter, y = `inside/outside.diameter`, color = group)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(face = "bold", size = 16),
    axis.text.y = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 16),
    axis.line = element_line(colour = "black", size = 1),
    legend.title = element_text(size = 16, face = "bold"),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(.2, "cm")
  ) +
  ylab("g-ratio") + xlab("Diameter [µm]")  +
  scale_x_continuous(breaks = c(0.4, 0.8, 1.2, 1.6))

ggsave("fimbria gratio vs diameter black line with text.png", last_plot(), width = 10, height = 10, dpi = 600)

ggplot(gratio, aes(x = Inside.diameter, y = `inside/outside.diameter`, color = group)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.text = element_text(size = 16),
    axis.line = element_line(colour = "black", size = 1),
    legend.title = element_text(size = 16, face = "bold"),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(.2, "cm")
  ) +
  ylab("") + xlab("")  +
  scale_x_continuous(breaks = c(0.4, 0.8, 1.2, 1.6))

ggsave("fimbria gratio vs diameter black line no text.png", last_plot(), width = 10, height = 10, dpi = 600)

#### Figure 4 ####
gratio -> gratio_hist
round(gratio_hist$Inside.diameter,digits = 1) -> gratio_hist$Inside.diameter

write.table(gratio_hist,"fimbria_Figure5j_distribution.txt",sep = "\t",col.names = T,row.names = F,quote = F)

# Create the histogram
p <- ggplot(gratio_hist, aes(x = Inside.diameter, fill = group)) +
  geom_histogram(aes(y = after_stat(count / tapply(count, group, sum)[group] * 100)), binwidth = 0.1, alpha = 0.5, position = "dodge") +
  scale_x_continuous(breaks = seq(min(gratio_hist$Inside.diameter), max(gratio_hist$Inside.diameter), by = 0.1)) +
  ylab("Percentage") +
  theme_minimal() +
  theme(axis.title.x=element_text(size=14,face="bold"),
        axis.title.y=element_text(size=14,face="bold"), # Increase size of Y axis label and make it bold
        axis.text=element_text(size=12), # Increase size of axis text
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(face = "bold"))+ # Make the group names bold
  ylab("axons [%]") + xlab("axons diameter [µm]")

# Add white vertical lines at each bin boundary
bin_boundaries <- seq(min(gratio_hist$Inside.diameter), max(gratio_hist$Inside.diameter), by = 0.1)
for (boundary in bin_boundaries) {
  p <- p + geom_vline(xintercept = boundary, color = "white", linewidth = 0.5)
}

print(p)
ggsave("histogram_axon_diameter.png", p, width = 10, height = 10, dpi = 600)

p <- ggplot(gratio_hist, aes(x = Inside.diameter, fill = group)) +
  geom_histogram(binwidth = 0.1, alpha = 0.5, position = "dodge") + # Removed after_stat function
  scale_x_continuous(breaks = seq(min(gratio_hist$Inside.diameter), max(gratio_hist$Inside.diameter), by = 0.1)) +
  ylab("Count") + # Changed label to "Count"
  theme_minimal() +
  theme(axis.title.x=element_text(size=20,face="bold"),
        axis.title.y=element_text(size=20,face="bold"), # Increase size of Y axis label and make it bold
        axis.text=element_text(size=18), # Increase size of axis text
        legend.text = element_text(size = 16), # Increase size of legend text
        legend.title = element_text(size = 16, face = "bold"), # Increase size and bold the legend title
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(face = "bold"))+ # Make the group names bold
  ylab("Axons Count") + xlab("Axons Diameter [µm]") # Updated Y axis label

# Add white vertical lines at each bin boundary
bin_boundaries <- seq(min(gratio_hist$Inside.diameter), max(gratio_hist$Inside.diameter), by = 0.1)
for (boundary in bin_boundaries) {
  p <- p + geom_vline(xintercept = boundary, color = "white", linewidth = 0.5)
}

print(p)
ggsave("histogram_axon_diameter_count.png", p, width = 10, height = 10, dpi = 600)

## Version 2 ##
# Create the histogram
p <- ggplot(gratio_hist, aes(x = Inside.diameter, fill = group)) +
  geom_histogram(binwidth = 0.1, alpha = 0.5, position = "dodge") +
  scale_x_continuous(breaks = seq(min(gratio_hist$Inside.diameter), max(gratio_hist$Inside.diameter), by = 0.1), expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Axons Count") +
  xlab("Axons Diameter [µm]") +
  theme_minimal() +
  theme(axis.title.x=element_text(size=20,face="bold"),
        axis.title.y=element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(face = "bold"),
        axis.ticks.x = element_line(color = "black", size = 0.5),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(.2, "cm"))

# Add white vertical lines at each bin boundary
bin_boundaries <- seq(min(gratio_hist$Inside.diameter), max(gratio_hist$Inside.diameter), by = 0.1)
for (boundary in bin_boundaries) {
  p <- p + geom_vline(xintercept = boundary, color = "white", linewidth = 0.5)
}

print(p)
ggsave("fimbria histogram_axon_diameter_count black line with text.png", p, width = 10, height = 10, dpi = 600)

p <- ggplot(gratio_hist, aes(x = Inside.diameter, fill = group)) +
  geom_histogram(binwidth = 0.1, alpha = 0.5, position = "dodge") +
  scale_x_continuous(breaks = seq(min(gratio_hist$Inside.diameter), max(gratio_hist$Inside.diameter), by = 0.1), expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(axis.title.x=element_text(size=20,face="bold"),
        axis.title.y=element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(color = "black", size = 0.5),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(.2, "cm"))

# Add white vertical lines at each bin boundary
bin_boundaries <- seq(min(gratio_hist$Inside.diameter), max(gratio_hist$Inside.diameter), by = 0.1)
for (boundary in bin_boundaries) {
  p <- p + geom_vline(xintercept = boundary, color = "white", linewidth = 0.5)
}
print(p)

ggsave("fimbria histogram_axon_diameter_count black line no text.png", p, width = 10, height = 10, dpi = 600)