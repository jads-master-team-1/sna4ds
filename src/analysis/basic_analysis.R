# Basic Data Analysis

## Load data
mep <- read.csv("./data/mep.csv", header = T, sep = ",")

## Summary statistics
summary(mep)

## Word frequency

### Political parties frequency table
freq_political_parties <- sort(table(unlist(strsplit(mep$group_abbv, " "))),
                               decreasing = TRUE)

pie(freq_political_parties)

### Country frequency table
freq_country <- sort(table(unlist(strsplit(mep$country, " "))),
                     decreasing = TRUE)

pie(freq_country)

### Voting frequency table
freq_voting <- sort(table(unlist(strsplit(mep$Paragraph.1..amendment.5, " "))),
                    decreasing = TRUE)

pie(freq_voting)

#### TODO: what to do with absence (from paper)

## MEP Names

#### TODO: some weird characters in the MEP names detected
