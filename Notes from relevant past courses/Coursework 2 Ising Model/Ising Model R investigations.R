library(ggplot2)
library(dplyr)
library(kableExtra)

# Testing section
data <- read.csv("testing.csv")
head(data)

data_modified <- data %>%
  filter(IV == 0) %>%
  group_by(Sweep) %>%
  summarise(avg_mag = mean(Magnetisation),
            avg_eng = mean(Energy))

data_2 <- data %>%
  group_by(Sweep, IV) %>%
  summarise(avg_mag = mean(Magnetisation),
            avg_eng = mean(Energy)) %>%
  ungroup() %>%
  group_by(IV) %>%
  summarise(min_mag = min(-avg_mag),
            min_eng = min(avg_eng))

ggplot(data_modified, aes(x = Sweep, y = -avg_mag)) +
  geom_line() +
  geom_point()

ggplot(data_modified, aes(x = Sweep, y = avg_eng)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0,4))

func <- function(iv_val) {
  tmp <- data %>%
      filter(IV == iv_val) %>%
      group_by(Sweep) %>%
      summarise(avg_mag = mean(Magnetisation),
              avg_eng = mean(Energy))
  ggplot(tmp, aes(x = Sweep, y = avg_eng)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(limits = c(0,4))
}

func(25)

ggplot(data_2, aes(x = IV, y = min_eng)) + 
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0,4))

ggplot(data_2, aes(x = IV, y = min_mag)) + 
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(-1,1))

# Question 1

data <- read.csv("IMQ1.csv")

# Shows the convergence to a steady state for the initial beta value
data_modified <- data %>%
  filter(Temperature == 2.25) %>%
  group_by(Sweep) %>%
  summarise(avg_mag = mean(Magnetisation),
            avg_eng = mean(Energy))

ggplot(data_modified, aes(x = Sweep, y = -avg_mag)) +
  geom_line() 


ggplot(data_modified, aes(x = Sweep, y = avg_eng)) +
  geom_line() 

tmp <- data %>%
  group_by(Sweep, Temperature) %>%
  summarise(avg_mag = mean(Magnetisation),
            avg_eng = mean(Energy)) %>%
  filter(Sweep != 0)

ggplot(tmp, aes(x = log(Sweep), y = avg_mag, color = factor(Temperature))) +
  geom_line() +
  scale_color_manual(values = c("#6f0000", "#8a101c", "#f62257",
                                "#f42794", "#d74fe4", "#a011ff",
                                "#8c63ff", "#0044ff", "#498afb")) +
  labs(color = "Temperature", x = "Ln(Sweep number)", y = "Average Magnetism") +
  theme_gray(base_size = 18)

ggplot(tmp, aes(x = log(Sweep), y = avg_eng, color = factor(Temperature))) +
  geom_line() +
  scale_color_manual(values = c("#6f0000", "#8a101c", "#f62257",
                                "#f42794", "#d74fe4", "#a011ff",
                                "#8c63ff", "#0044ff", "#498afb")) +
  labs(color = "Temperature", x = "Ln(Sweep number)", y = "Average Energy") +
  theme_gray(base_size = 18)


# Question 2

dataq2 <- read.csv("IMQ2.csv")

tmp <- dataq2 %>%
  group_by(Temperature) %>%
  summarise(avg_mag = round(mean(Magnetisation), 5),
            avg_eng = round(mean(Energy), 5),
            avg_abs = round(mean(abs(Magnetisation)), 5),
            var_mag = round(var(Magnetisation), 5),
            var_eng = round(var(Energy), 5)) %>%
  select(Temperature, avg_eng, var_eng, avg_mag, avg_abs, var_mag)
tmp$Temperature <- c("1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "4.0")

kable(tmp, format = "latex", align = "cc", escape = FALSE, booktabs = TRUE,
      caption = "Insert Caption",
      col.names = c("$T$", "$\\langle E \\rangle$", "Var($E$)", "$\\langle M \\rangle$", "$\\langle |M| \\rangle$", "Var($M$)")) %>%
  kable_styling(latex_options = "hold_position", full_width = FALSE)

# Question 3

dataq3 <- read.csv("IMQ3.csv")
datagiven <- read.csv("exactisingdata.csv")

tmp <- dataq3 %>%
  group_by(Temperature) %>%
  summarise(avg_mag = mean(Magnetisation),
            avg_eng = (mean(Energy/2)),
            avg_abs = mean(abs(Magnetisation)),
            var_mag = var(Magnetisation),
            var_eng = var(Energy/2)) %>%
  # rename(T = Temperature) %>%
  # mutate(var_eng = var_eng * 5000)
    mutate(avg_eng = avg_eng * 2)



tmp2 <- datagiven %>%
  rename(Temperature = T, avg_mag = m, avg_eng = e, var_eng = c) %>%
  mutate(avg_mag = -1 * avg_mag) 

ggplot(tmp, aes(x = Temperature, y = avg_mag)) +
  geom_line(color = "blue") +
  geom_line(data = tmp2, color = "red") +
  labs(x = "Temperature", y = "Magnetism") +
  theme_gray(base_size = 18)

ggplot(tmp, aes(x = Temperature, y = avg_eng)) +
  geom_line(color = "blue") +
  geom_line(data = tmp2, color = "red") +
  labs(x = "Temperature", y = "Energy") +
  theme_gray(base_size = 18)

ggplot(tmp, aes(x = Temperature, y = var_eng)) +
  geom_line(color = "blue") +
  geom_line(data = tmp2, color = "red") +
  labs(x = "Temperature", y = "Heat capacity") +
  theme_gray(base_size = 18)

# test <- cross_join(tmp, datagiven)

# ggplot(test, aes(x = var_eng, y = c, color = T)) +
#   geom_point() +
#   geom_smooth(se = FALSE, method = lm)


# summary(lm(data = test, c ~ var_eng))
#  1.775e+03 - 6.47e+01

# Question 4

dataq4 <- read.csv("IMQ4.csv")
dataq5 <- read.csv("IMQ5.csv")

tmp <- dataq4 %>%
  summarise(avg_cor = mean(Distance))

tmp2 <- dataq5 %>%
  filter(Distance == 1) %>%
  summarise(avg_cor = mean(Correlation))


ggplot(tmp, aes(x = Temperature, y = avg_mag)) +
  geom_line(color = "blue") +
  geom_line(data = tmp2, color = "red") +
  labs(x = "Temperature", y = "Magnetism") +
  theme_gray(base_size = 18)

ggplot(tmp, aes(x = Temperature, y = avg_eng)) +
  geom_line(color = "blue") +
  geom_line(data = tmp2, color = "red") +
  labs(x = "Temperature", y = "Energy") +
  theme_gray(base_size = 18)

ggplot(tmp, aes(x = Temperature, y = var_eng)) +
  geom_line(color = "blue") +
  geom_line(data = tmp2, color = "red") +
  labs(x = "Temperature", y = "Heat capacity") +
  theme_gray(base_size = 18)

test <- left_join(tmp, datagiven)

ggplot(test, aes(x = var_eng, y = c, color = T)) +
  geom_point() +
  geom_smooth(se = FALSE, method = lm)


summary(lm(data = test, c ~ var_eng))

# Question 5

# Question 6

data6 <- read.csv("IMQ6.csv")

# Shows the convergence to a steady state for the initial beta value
data_modified <- data6 %>%
  group_by(Sweep) %>%
  summarise(avg_mag = mean(Magnetisation),
            avg_eng = mean(Energy))

ggplot(data_modified, aes(x = Sweep, y = -avg_mag)) +
  geom_line() 


ggplot(data_modified, aes(x = Sweep, y = avg_eng)) +
  geom_line() 

tmp <- data %>%
  group_by(Sweep, Temperature) %>%
  summarise(avg_mag = mean(Magnetisation),
            avg_eng = mean(Energy)) %>%
  filter(Sweep != 0)

ggplot(tmp, aes(x = log(Sweep), y = avg_mag, color = factor(Temperature))) +
  geom_line() +
  scale_color_manual(values = c("#6f0000", "#8a101c", "#f62257",
                                "#f42794", "#d74fe4", "#a011ff",
                                "#8c63ff", "#0044ff", "#498afb")) +
  labs(color = "Temperature", x = "Ln(Sweep number)", y = "Average Energy") +
  theme_gray(base_size = 18)

# End of R file

