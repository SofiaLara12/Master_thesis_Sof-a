library(readr)
sexdetresults <- read_csv("sexdets_results.csv", col_types = cols(.default = "c"))

library(dplyr)
library(stringr)

sexdets <- sexdetresults %>%
  mutate(
    Genus = word(Species, 1),
    Epithet = word(Species, 2)
  )

sexdets <- sexdets %>%
  mutate(
    Genus = word(Species, 1),
    Epithet = word(Species, 2)
  ) %>%
  relocate(Genus, Epithet, .before = 1)

sexdets <- sexdets %>%
  mutate(Sex = case_match(Sex,
                          "U" ~ "Unknown",
                          "F" ~ "Female",
                          "M" ~ "Male",
                          "M*" ~ "Male",
                          .default = Sex
  ))


sexdets_known <- sexdets %>%
  filter(!Sex %in% c("Unknown", "unknown"))

sexdets_known <- sexdets_known[-(1:9), ]

sexdetsstrict <- sexdets_known2 %>%
  filter(PCR_2.1 == PCR_2.2)

## caret 
library(caret)
sexdetscaret <- sexdetsstrict %>%
  mutate(
    # Si la PCR dio Yes -> Male; si dio No -> Female
    PCR_Result = case_when(
      PCR_2.1 == "Yes" ~ "Male",
      PCR_2.1 == "No"  ~ "Female",
      TRUE            ~ NA_character_
    )
  )

sexdetscaret$PCR_Result <- factor(sexdetscaret$PCR_Result, levels = c("Male", "Female"))
sexdetscaret$Sex        <- factor(sexdetscaret$Sex,        levels = c("Male", "Female"))

cm <- confusionMatrix(
  data = sexdetscaret$PCR_Result,
  reference = sexdetscaret$Sex,
  positive = "Male"
)

print(cm)

library(ggplot2)

cm_tbl <- as.data.frame(cm$table)

cm_tbl <- cm_tbl %>%
  mutate(
    Type = case_when(
      Prediction == "Male"   & Reference == "Male"   ~ paste0("TP=", Freq),
      Prediction == "Female" & Reference == "Female" ~ paste0("TN=", Freq),
      Prediction == "Male"   & Reference == "Female" ~ paste0("FP=", Freq),
      Prediction == "Female" & Reference == "Male"   ~ paste0("FN=", Freq)
    ),
    Category = ifelse(Prediction == Reference, "Correct", "Incorrect")
  )

p <- ggplot(cm_tbl, aes(x = Prediction, y = Reference, fill = Category)) +
  geom_tile(color = "gray80", linewidth = 1) +
  geom_text(aes(label = Type), fontface = "bold", size = 5) +
  scale_fill_manual(values = c("Correct" = "#E69F00", "Incorrect" = "#56B4E9")) +
  labs(
    x = "PCR test prediction",
    y = "Identified sex"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12, face = "plain"),
    axis.text = element_text(size = 11, color = "black"),
    panel.grid = element_blank()
  )

ggsave("confusion_matrix_plot.png", plot = p, width = 6, height = 5, dpi = 300)

sexdets_reprod <- sexdets %>%
  filter(rowSums(select(., starts_with("PCR_")) == "/", na.rm = TRUE) < 4)

library(dplyr)

sexdets_reprod <- sexdets %>%
  filter(
    rowSums(
      across(starts_with("PCR_"), ~ .x %in% c("Yes", "No", "YES", "NO")), 
      na.rm = TRUE
    ) >= 2
  )

library(dplyr)

sexdets_patrones <- sexdets_reprod %>%
  mutate(across(starts_with("PCR_"), ~ ifelse(.x == "/", NA, .x))) %>%
  rowwise() %>%
  mutate(
    n_Yes = sum(c_across(starts_with("PCR_")) == "Yes", na.rm = TRUE),
    n_No  = sum(c_across(starts_with("PCR_")) == "No",  na.rm = TRUE),
    n_Total = n_Yes + n_No,
        Patron_Concordancia = case_when(
      (n_Yes > 0 & n_No == 0) | (n_Yes == 0 & n_No > 0) ~ "100% Consistente",
      n_Yes == n_No ~ "50 / 50",
      (n_Yes == 2 & n_No == 1) | (n_Yes == 1 & n_No == 2) ~ "2 iguales y 1 diferente",
      (n_Yes == 3 & n_No == 1) | (n_Yes == 1 & n_No == 3) ~ "3 iguales y 1 diferente",
      TRUE ~ "Otro"
    )
  ) %>%
  ungroup()

resumen_patrones <- sexdets_patrones %>%
  count(Patron_Concordancia, name = "Num_Muestras") %>%
  mutate(Porcentaje = round((Num_Muestras / sum(Num_Muestras)) * 100, 2))

print(resumen_patrones)

library(dplyr)

sexsdetsfinal <- sexdets_patrones %>%
  mutate(
    SEXRESULT = case_when(
      n_Yes > 0 & n_No == 0 ~ "Male",
      n_No > 0 & n_Yes == 0 ~ "Female",
      n_Yes == n_No ~ "Indeterminate",
      n_Yes == 3 & n_No == 1 ~ "Male",
      n_No == 3 & n_Yes == 1 ~ "Female",
      n_Yes == 2 & n_No == 1 ~ "Male",
      n_No == 2 & n_Yes == 1 ~ "Female",
      TRUE ~ "Indeterminate"
    )
  )

table(sexdetsfinalDUTCH$SEXRESULT, useNA = "ifany")

generos_femaleD <- sexdetsfinalDUTCH %>%
  filter(SEXRESULT == "Female") %>%
  distinct(Genus) %>%
  pull(Genus)

generos_maleD <- sexdetsfinalDUTCH %>%
  filter(SEXRESULT == "Male") %>%
  distinct(Genus) %>%
  pull(Genus)

generos_indeterminateD <- sexdetsfinalDUTCH %>%
  filter(SEXRESULT == "Indeterminate") %>%
  distinct(Genus) %>%
  pull(Genus)
