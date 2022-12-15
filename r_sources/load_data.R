library(ggplot2)
library(plotly)
library(dplyr)
library(rstan)
library(httpgd)
library(gridExtra)
library(outbreaks)
library(bayesplot)
library(data.table)
library(knitr)
library(kableExtra)
library(reticulate)
library("styler")
load_data <- function(path = "data", location = "cdmx")
{
    sub_path_1 <- getwd()
    tree_path <-
        paste(
            sub_path_1,
            path,
            sep = "/"
        )
    file_name <- "cdmx_prevalence_data.csv"
    data_path <- paste(tree_path, file_name, sep = "/")
    covid19_data <- fread(data_path,
                                    select = c("FECHA_SINTOMAS",
                                              "i_s",
                                              "cumulative_i_s"))
    covid19_data <- data.frame(covid19_data)
#
    reference_date <- as.Date('2020-03-10')
    final_date_sample <- as.Date('2020-03-30')
    data_star_dynamics <- covid19_data %>%
        filter(as.Date(FECHA_SINTOMAS) >= reference_date &
                   as.Date(FECHA_SINTOMAS) <= final_date_sample)
    head(data_star_dynamics)
    data_plot <- ggplot(data = data_star_dynamics,
                        aes(x = FECHA_SINTOMAS, cumulative_i_s)) +
                    geom_bar(stat="identity", width = 0.05) +
                    geom_point() +
                    theme(axis.text.x = element_text(angle = 90)) +
                    ggtitle("CDMX")
    #
    data_plotly <- ggplotly(data_plot)
    file_name_pdf <- "cdmx_input_data.pdf"
    file_name_html <- "cdmx_input_data.html"
    file_name_png <- "cdmx_input_data.png"
    path <- "plots"
    plot_path_pdf <-
        paste(sub_path_1, path, file_name_pdf, sep = "/")
    plot_path_html <-
        paste(sub_path_1, path, file_name_html, sep = "/")
    plot_path_png <-
        paste(sub_path_1, path, file_name_png, sep = "/")
    ggsave(plot_path_pdf)

    htmlwidgets::saveWidget(as_widget(data_plotly), plot_path_html)
    if (!require("processx")) {
      install.packages("processx")
    }
    # orca(fig, "figure_01.png")
    save_image(data_plotly, plot_path_png, width = 1417, height = 875.7726)

    onset <- data_star_dynamics %>%
        select(FECHA_SINTOMAS)
    cum_cases <- data_star_dynamics %>%
        select(cumulative_i_s)

    cases <- data_star_dynamics %>%
        select(i_s)
    # cum_cases <- unlist(cum_cases, use.names = FALSE)
    names(cum_cases)[1] <-"cum_cases"
    data <- list(onset, cum_cases, data_star_dynamics)
    return(data)
}
