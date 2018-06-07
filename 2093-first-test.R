# Script to measure baseline noise in Spectroradiometer S/N 2092

# install required packages
if (!'googledrive' %in% installed.packages() ) install.packages('googledrive')
if (!'stringr' %in% installed.packages() ) install.packages('stringr')
if (!'spectrolab' %in% installed.packages() ) install.packages('spectrolab')
if (!'tidyverse' %in% installed.packages() ) install.packages('tidyverse')


# load packages
library(googledrive)
library(stringr)
library(spectrolab)
library(tidyverse)


# function to replace commas by dots on individual metadata lines
comma_to_dot <- function(x) {
  tmp <- unlist(strsplit(x, ", ") )
  tmp2 <- str_replace(tmp, ',', '.')
  paste(tmp2, collapse = ', ')
}

#### 2093 wd
folder2093 <- paste0(getwd(), '/2093')
dir.create(folder2093)
setwd(folder2093)

# get list of .sig files in the test working directory (Google Drive)
files <- drive_ls(path = '2018-Girard-MSc-UdeM/spectra/2018-05-24-JBMCBTests-2093')

# create directory to store raw spectra
raw_folder <- paste0(getwd(), '/spectra_raw')
dir.create(raw_folder)
for (i in 1:nrow(files)) drive_download(files[i, ], path = paste0(getwd(), '/spectra_raw', '/', files$name[i]), overwrite = T)

temp_folder <- paste0(getwd(), '/spectra_temp')
dir.create(temp_folder)

# convert commas for all .sig files (were saved using the wrong decimal separator...)
for (i in 1:nrow(files)) {
  file1 <- file(paste0(getwd(), '/spectra_raw', '/', files$name[i]), open = 'r')
  file2 <- file(paste0(temp_folder, '/', files$name[i]), open = 'w')
  part1a <- readLines(file1, n = 3)
  part1b <-  comma_to_dot(readLines(file1, n = 1) ) # integration
  part1c <- readLines(file1, n = 9)
  part1d <-  comma_to_dot(readLines(file1, n = 1) ) # temp
  part1e <-  comma_to_dot(readLines(file1, n = 1) ) # battery
  part1f <- readLines(file1, n = 8)
  part1g <-  comma_to_dot(readLines(file1, n = 1) ) # factors
  readLines(file1, n = 2) # remove the two inclinometer lines
  part1h <- readLines(file1, n = 1)
  part1 <- c(part1a, part1b, part1c, part1d, part1e, part1f, part1g, part1h)
  part2 <- readLines(file1)
  part2 <- gsub(',', '.', part2)
  
  close(file1)
  #file1 <- file(files$name[1], open = 'r+')
  writeLines(part1, file2)
  writeLines(part2, file2)
  #writeLines('', file1)
  close(file2)
}

# read .sig files
spectra <- read_spectra(temp_folder, format = 'sig')
plot(spectra)

### store spectra in data frame for plotting
spectra.df <- as.data.frame(spectra)

# add "minute" factor
spectra.df$min <- factor(paste(substr(as.character(spectra.df$sample_name), 10, 13), 'min' ) )
spectra.df$min.num <- as.numeric(substr(as.character(spectra.df$sample_name), 10, 13))

# convert to long format
spectra.df.long <- gather(spectra.df, wvl, refl, `337.3`:`2513.1`) %>%
  mutate(wvl = as.numeric(wvl)) %>%
  arrange(sample_name, min, wvl)


# baseline 0000
spectra.0 <- filter(spectra.df.long, min == '0000 min' )
head(spectra.0) 

# plot all
spectra_all <- ggplot(spectra.df.long, aes(x = wvl, y = refl)) +
  geom_line(aes(group = min)) +
  facet_wrap(~min, ncol = 5) +
  ylab('Relative reflectance') +
  xlab('Wavelength (nm)') +
  scale_y_continuous(limits = c(0.85, 1.15)) +
  theme_bw()
ggsave(spectra_all, file = 'spectra_all.pdf', height = 65, width = 12, limitsize = F)


# re-calculate ref at 15 min
spectra_15 <- filter(spectra.df.long, min == '0015 min' )
head(spectra_15)

# divide all spectra taken on or after 15 min by ref at 15 min
spectra_from15 <- spectra.df.long %>%
  filter(min.num >= 15) %>%
  group_by(sample_name) %>%
  mutate(refl.corr = refl / spectra_15$refl)
head(spectra_from15)     


# plot all
spectra_plot_from15 <- ggplot(spectra_from15, aes(x = wvl, y = refl.corr)) +
  geom_line(aes(group = min)) +
  facet_wrap(~min, ncol = 3) +
  ylab('Relative reflectance') +
  xlab('Wavelength (nm)') +
  scale_y_continuous(limits = c(0.85, 1.15))
ggsave(spectra_plot_from15, file = 'spectra_plot_from15.pdf', height = 55, width = 6, limitsize = F)


# custom function to divide
div_min <- function(x) {
  ref <- min(x$min.num)
  x.sub <- subset(x, min.num == ref)
  refl.sub <- x.sub$refl
  x <- mutate(x, refl.corr = refl / refl.sub)
  return(x)
}

# do it for every 5 min
spectra_every5 <- spectra.df.long %>%
  mutate(group = ceiling((min.num + 1 ) / 5) )  %>%
  group_by(group) %>%
  do(div_min(.))
head(spectra_every5) 


# plot all - every 5 min
spectra_plot_every5 <- ggplot(spectra_every5, aes(x = wvl, y = refl.corr)) +
  geom_line(aes(group = min)) +
  geom_hline(yintercept = 0.98, linetype = 'dotted') +
  geom_hline(yintercept = 1.02, linetype = 'dotted') +
  geom_vline(xintercept = 2300, linetype = 'dotted') +
  facet_wrap(~min, ncol = 5) +
  ylab('Relative reflectance') +
  xlab('Wavelength (nm)') +
  scale_y_continuous(limits = c(0.85, 1.15)) +
  theme_bw()
ggsave(spectra_plot_every5, file = 'spectra_plot_every5.pdf', height = 65, width = 12, limitsize = F)

    
# calculate mean + sd
spectra_every5_noref <- spectra_every5 %>%
  filter(!min.num %in% (seq(0, 130, 5) ) )

spectra_5_mean_sd <- spectra_every5_noref %>%
  group_by(group, wvl) %>%
  summarise(refl.corr2 = mean(refl.corr, na.rm = T), refl.sd = sd(refl.corr, na.rm = T))

# plot all - every 5 min
spectra_plot_every5_meansd <- ggplot(spectra_5_mean_sd, aes(x = wvl ) ) +
 # geom_line(aes(y = refl.corr2))) +
 geom_ribbon(aes(ymin = 1 -1.96 * refl.sd, ymax = 1 + 1.96 * refl.sd )) +
  facet_wrap(~group, ncol = 5) +
  ylab('Relative reflectance') +
  xlab('Wavelength (nm)') +
  scale_y_continuous(limits = c(0.85, 1.15))
ggsave(spectra_plot_every5_meansd, file = 'spectra_plot_every5_meansd.pdf', height = 15, width = 12, limitsize = F)

# using geom_smooth
spectra_plot_every5_smooth <- ggplot(spectra_every5_noref, aes(x = wvl, y = refl.corr ) ) +
   geom_smooth() +
  facet_wrap(~group, ncol = 5) +
  ylab('Relative reflectance') +
  xlab('Wavelength (nm)') +
  scale_y_continuous(limits = c(0.85, 1.15))
ggsave(spectra_plot_every5_smooth, file = 'spectra_plot_every5_smooth.pdf', height = 15, width = 12, limitsize = F)
