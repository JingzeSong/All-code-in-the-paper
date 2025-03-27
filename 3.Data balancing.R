setwd("Work path")

df <- read.csv("mRNA_count.csv", 
               header = TRUE, 
               row.names = 1,
               stringsAsFactors = FALSE)

df_part1 <- df[, 1:99]

# Seeds can be specified if a fixed random result is required, which results in a slightly different calculation of the minimum classification error rate and has less impact on final method validation
set.seed(123)  
random_cols <- sample(100:ncol(df), 99)
df_random <- df[, random_cols]

df_final <- cbind(df_part1, df_random)

write.csv(df_final, "mRNA_count_new.csv", row.names = TRUE)
