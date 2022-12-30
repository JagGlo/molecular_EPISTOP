load(("~/EPISTOP/models/seizures_normalized_c3.RData"))
source(here::here("~/interim_report/models/3_read_chunks.R"))
library(plot3D)

pbapply::pblapply(1:as.numeric(count(models_to_eval)), function(i) {
  browser()
  vars <-
    models_to_eval[[i, "features"]] %>%
    str_split(pattern = "\\s\\+\\s") %>%
    unlist()

  fields <-
    vars %>%
    match(names(raw_values))

  data <-
    raw_values[, c(fields, ncol(raw_values)-2, ncol(raw_values)-1, ncol(raw_values))] %>%
    filter(complete.cases(.)) %>% 
    mutate(seizures = if_else(seizures_any == TRUE, "seizures", "no seizures") %>% as.factor() %>% fct_relevel("seizures", "no seizures")) %>%
    select(-seizures_any)
  
  simpler_model <- glm(formula = data$seizures_any ~ pull(data[1]) + pull(data[2]) + pull(data[3]),
                       family = "binomial")      
  
  plot3D::scatter3D(data[1], data[2], data[3], col.var = data$seizures)
  
  # data %>% 
  #   ggplot(aes(x =as.formula(vars[1]), y = as.formula(vars[1]), color = seizures)) +
  #   geom_point() + 
  #   theme_minimal()

  p <-
    plotly::plot_ly(data,
      x = as.formula(paste0("~", vars[1])),
      y = as.formula(paste0("~", vars[2])),
      z = as.formula(paste0("~", vars[3])),
      color = ~seizures,
      colors = c("#0C4B8E", "#BF382A"),
      type = "scatter3d",
      mode = "markers",
      text = ~paste0("Patient: ", code, ", ", sample_code)
    )
  # 
  # library(plotly)
  # library(magrittr)
  # 
  # plot_ly(data = data) %>%
  #   add_trace(
  #     x = as.formula(paste0("~", vars[1])),
  #     y = as.formula(paste0("~", vars[2])),
  #     z = ~seizures_any,
  #     color = ~seizures_any,
  #     colors = c("#0C4B8E", "#BF382A"),
  #     type = "scatter3d") %>% 
  #   add_trace(z = simpler_model$fitted.values, x = as.formula(paste0("~", vars[1])), y = as.formula(paste0("~", vars[2])), type = "mesh3d", 
  #             name = "Fitted values") %>%
  #   layout(scene = list(xaxis = list(title = 'sales'), yaxis = list(title = 'customer_rate', nticks = 5),
  #                       camera = list(eye = list(x = -0.5, y = 2, z = 0)),
  #                       zaxis = list(title = 'promoted'), aspectmode='cube')) 

  htmlwidgets::saveWidget(
    widget = p,
    paste0(
      "ra/mnt/epistop/plots/",
      str_pad(i, width = 4, side = "left", pad = "0"),
      ".html"
    ),
    selfcontained = FALSE,
    libdir = "/mnt/epistop/plots/static"
  )
})

dataset %>% select(seizures_any, everything())  %>% summarise_all(funs(sum(!is.na(.)))) %>% gather(key, n) %>% 
  ggplot(aes(n)) + geom_bar() 


          dataset %>% group_by(seizures_any) %>%  select(everything()) %>%  # replace to your needs
                summarise_all(funs(sum(!is.na(.)))) %>% gather(key, value, -seizures_any) %>% 
            ggplot(aes(value)) + geom_bar(aes(fill = seizures_any)) + facet_wrap(~ seizures_any, scales = "free")

